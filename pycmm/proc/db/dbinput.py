import ntpath
import itertools
import re
import os.path
from collections import OrderedDict
from vcf.model import _Call as _VcfCall
from vcf.model import _Record as _VcfRecord
from vcf import RESERVED_INFO
from vcf import Reader as VcfReader
from pycmm.utils import is_number
from pycmm.template import pyCMMBase
from pycmm.cmmlib.readerlib import Reader

AVDB_CHR_IDX = 0
AVDB_START_IDX = 1
AVDB_END_IDX = 2
AVDB_REF_IDX = 3
AVDB_ALT_IDX = 4
AVDB_ANNOS_IDX = 5

DATA_TYPE_AVDB_CMM_AF = "DATA_TYPE_AVDB_CMM_AF"
DATA_TYPE_AVDB_PUBLIC_AF = "DATA_TYPE_AVDB_PUBLIC_AF"

CMMGT_WILDTYPE = 'wt'
CMMGT_HOMOZYGOTE = 'hom'
CMMGT_HETEROZYGOTE = 'het'
CMMGT_OTHER = 'oth'
CMM_DUMMY = 'dummy'

class AVDBParser(pyCMMBase):
    """ To parse a record in Annovar AVDB file """

    def __init__(self, rec, *args, **kwargs):
        self.__items = rec.strip().split("\t")
        if len(self.__items) < 6:
            raise IOError("Invalid avdb file")
        super(AVDBParser, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self):
        raw_str = OrderedDict()
        raw_str["Chrom"] = self.chrom
        raw_str["Start"] = self.start
        raw_str["End"] = self.end
        raw_str["Ref"] = self.ref
        raw_str["Alt"] = self.alt
        raw_str["Annotations"] = self.annotations
        return raw_str

    @property
    def chrom(self):
        return self.__items[AVDB_CHR_IDX]

    @property
    def start(self):
        return int(self.__items[AVDB_START_IDX])

    @property
    def end(self):
        return int(self.__items[AVDB_END_IDX])

    @property
    def ref(self):
        return self.__items[AVDB_REF_IDX]

    @property
    def alt(self):
        return self.__items[AVDB_ALT_IDX]

    @property
    def annotations(self):
        return self.__items[AVDB_ANNOS_IDX:]

    def as_tuple(self, data_type, gf_idx=None):
        res = []
        res.append(self.chrom)
        res.append(self.start)
        res.append(self.end)
        res.append(self.ref)
        res.append(self.alt)
#        self.dbg(data_type)
        if data_type == DATA_TYPE_AVDB_CMM_AF:
            res += map(lambda x: int(x),
                       self.annotations[:gf_idx])
            res += [float(self.annotations[gf_idx])]
            af = self.annotations[gf_idx+1]
            if is_number(af):
                af = float(af)
            res += [af]
            res += [float(self.annotations[gf_idx+2])]
        elif data_type == DATA_TYPE_AVDB_PUBLIC_AF:
            res += map(lambda x: float(x), self.annotations)
        else:
            res += self.annotations
        return tuple(res)

class AVDBReader(Reader):
    """ To read and parse file plink map file """

    def __init__(self, 
                 data_type=None,
                 *args,
                 **kwargs
                 ):
        kwargs['parser'] = AVDBParser
        self.__data_type = data_type
        file_name = kwargs['file_name']
        super(AVDBReader, self).__init__(*args, **kwargs)
        self.__header_cols = self.__parse_header_cols(file_name)
        self.__gf_idx = 0
        if data_type == DATA_TYPE_AVDB_CMM_AF:
            for col_name in self.header_cols:
                if col_name.endswith("GF"):
                    self.__gf_idx = self.header_cols.index(col_name)
                    self.__gf_idx -= len(self.pkeys)
                    self.__gf_idx -= 1

    def __parse_header_cols(self, file_name):
        header_cols = []
        header_cols.append('Chr')
        header_cols.append('Start')
        header_cols.append('End')
        header_cols.append('Ref')
        header_cols.append('Alt')
        if self.header is None:
            header_cols.append(os.path.splitext(ntpath.basename(file_name).strip('hg19_'))[0]) 
        else:
            split_header_cols = self.header.strip().split("\t")
            for col_idx in xrange(5, len(split_header_cols)):
                col_name = split_header_cols[col_idx]
                col_name = col_name.replace(":", "_")
                col_name = col_name.replace("+", "")
                col_name = col_name.replace("-", "")
                header_cols.append(col_name)
        return header_cols

    @property
    def header_cols(self):
        return self.__header_cols

    @property
    def avdb_cols(self):
        return self.header_cols[5:]

    @property
    def pkeys(self):
        pkeys = []
        pkeys.append('Chr')
        pkeys.append('Start')
        pkeys.append('Ref')
        pkeys.append('Alt')
        return pkeys

#    def next(self):
#        return super(AVDBReader, self).next().as_tuple()
#
    @property
    def record_tuples(self):
        for item in super(AVDBReader, self).__iter__():
            yield item.as_tuple(self.__data_type, self.__gf_idx)
#        return super(AVDBReader, self).next().as_tuple()

try:
    import cparse
except ImportError:
    cparse = None

class _TAVcfCall(_VcfCall, pyCMMBase):
    """
    An encapsulated version of vcf._Call from pyVCF package to
      - add extra genotype translation, "cmm_gts", to handle record
        with more than one alternate alleles
      - add extra genotype translation, "actual_gts", to determine the
        actual genotype based on "cmm_gts" and alleles frequency
      - add an indicator, "mutated", to identify if each genotype is
        actually mutated
    """

    def __init__(self, site, sample, data):
        _VcfCall.__init__(self,
                          site=site,
                          sample=sample,
                          data=data,
                          )
        pyCMMBase.__init__(self)
        self.__cmm_gts = None

    def __cal_cmm_gts(self):
        """
        to calculate type of genotype given allele index
          - 0 -> REF
          - 1 -> first ALT
          - and so on

        ** NOTE ** the calculation here based on GT value alone
        """
        raw_GT = self.data.GT
        cmm_gts = [CMM_DUMMY]
        for allele_idx in xrange(1, len(self.site.alleles)):
            if (raw_GT == ".") or (raw_GT == "./."):
                cmm_gts.append(".")
                continue
            raw_gts = raw_GT.split("/")
            # 0/0 will be translated as "wild type" no matter what allele index is
            if (raw_gts[0] == "0") and (raw_gts[1] == "0"):
                cmm_gts.append(CMMGT_WILDTYPE)
                continue
            # The gt of which the allele index doesn't match will be considered
            # as "other", ex. allele index = 1 but the gt is 2/3
            if (raw_gts[0] != str(allele_idx)) and (raw_gts[1] != str(allele_idx)):
                cmm_gts.append(CMMGT_OTHER)
                continue
            # The rest are in the cases one or both of the alleles match with
            # allele index.
            # So if they are different then "heterozygote" else "homozygote"
            if raw_gts[0] != raw_gts[1]:
                cmm_gts.append(CMMGT_HETEROZYGOTE)
                continue
            cmm_gts.append(CMMGT_HOMOZYGOTE)
        return cmm_gts

    @property
    def cmm_gts(self):
        if self.__cmm_gts is None:
            self.__cmm_gts = self.__cal_cmm_gts()
        return self.__cmm_gts

class _TAVcfRecord(_VcfRecord, pyCMMBase):
    """
    An encapsulated version of vcf._Record from pyVCF package to
      - determine if mutations are shared between samples
      - understand if itself is a rare mutation given frequency ratio
      - understand if itself is an intergenic mutation
      - understand if itself is an intronic mutation
    """

    def __init__(self,
                 *args,
                 **kwargs
                 ):
        super(_TAVcfRecord, self).__init__(*args, **kwargs)
        pyCMMBase.__init__(self)

        if (type(self.FILTER) is list) and (len(self.FILTER) == 0):
            self.FILTER = "PASS"
        elif type(self.FILTER) is list:
            self.FILTER = ";".join(self.FILTER)
        else:
            self.FILTER = "."

    def __get_info(self, var_name):
        """
        for internal call
        - return list of values of table_annovar column "var_name"
        - number of entry in the list = 1 + number of alternate alleles
        """
        if var_name in self.INFO:
            return self.INFO[var_name]
        if var_name in self.__added_info:
            return self.__added_info[var_name]
        return None

    def get_info(self, var_name, allele_idx=1):
        """
        parsed to be called by high-level function
        - kvot is included
        """

        info = self.__get_info(var_name)
        if info is not None:
            if (type(info) is list) and (len(info) == 1):
                info = info[0]
            elif (type(info) is list) and (len(info) > 1):
                info = info[allele_idx-1]
        if (info == "" or
            info is None or
            info == [None] or
            info == "."
            ):
           info = ""
        return info

    def __check_alt(self):
        if len(self.ALT) > 1:
            raise Exception(self.ALT + ' cannot have than one alternate allele')

    def __get_pkey_list(self):
        res = []
        res.append(self.CHROM)
        res.append(int(self.POS))
        res.append(self.REF)
        res.append(str(self.ALT[0]))
        return res

    def as_annovar_tuple(self, annovar_cols):
        self.__check_alt()
        res = self.__get_pkey_list()
        for annovar_col in annovar_cols:
            res.append(self.get_info(annovar_col))
        return tuple(res)

    def as_gtz_tuple(self, samples_id):
        self.__check_alt()
        res = self.__get_pkey_list()
        res.append(self.QUAL)
        res.append(self.FILTER)
        for sample_id in samples_id:
            res.append(self.genotype(sample_id).cmm_gts[1])
        return tuple(res)

    def as_pkey_tuple(self):
        self.__check_alt()
        return tuple(self.__get_pkey_list())

class TAVcfReader(VcfReader, pyCMMBase):
    """
    An encapsulated version of vcf.Reader from pyVCF package to
      - enable _parse_info to recognize INFO annotated by annovar
        , <see _parse_info function>
      - redirect _Call to TAVcfCall
      - redirect _Record to TAVcfRecord
    """

    def __init__(self,
                 file_name=None,
                 *args,
                 **kwargs
                 ):
        kwargs['filename'] = file_name
        super(TAVcfReader, self).__init__(*args, **kwargs)
        pyCMMBase.__init__(self)
        self.__parse_annovar_infos()

    @property
    def annovar_infos(self):
        return self.__annovar_infos

    @property
    def pkeys(self):
        pkeys = []
        pkeys.append('CHROM')
        pkeys.append('POS')
        pkeys.append('REF')
        pkeys.append('ALT')
        return pkeys

    def __parse_annovar_infos(self):
        """
        To get the list of infos field annotated by ANNOVAR
        by detecting any INFO in metaheader that end with "ANNOVAR"
        """
        self.__annovar_infos = OrderedDict()
        for info_key in self.infos:
            val = self.infos[info_key]
            if (val.desc is not None) and (val.desc.endswith("ANNOVAR")):
                self.__annovar_infos[info_key] = val
        return self.__annovar_infos

    def next(self):
        '''
        Remake version of pyVCF.parser._parse_info

        Return the next record in the file.

        ***** Remake part *****
        - use custom _Record structure, _TAVcfRecord
        '''
        line = self.reader.next()
        row = re.split(self._separator, line.rstrip())
        chrom = row[0]
        if self._prepend_chr:
            chrom = 'chr' + chrom
        pos = int(row[1])

        if row[2] != '.':
            ID = row[2]
        else:
            ID = None

        ref = row[3]
        alt = self._map(self._parse_alt, row[4].split(','))

        try:
            qual = int(row[5])
        except ValueError:
            try:
                qual = float(row[5])
            except ValueError:
                qual = None

        filt = row[6]
        if filt == '.':
            filt = None
        elif filt == 'PASS':
            filt = []
        else:
            filt = filt.split(';')
        info = self._parse_info(row[7])

        try:
            fmt = row[8]
        except IndexError:
            fmt = None
        else:
            if fmt == '.':
                fmt = None

        record = _TAVcfRecord(chrom, pos, ID, ref, alt, qual, filt,
                info, fmt,
                self._sample_indexes,
                )

        if fmt is not None:
            samples = self._parse_samples(row[9:], fmt, record)
            record.samples = samples

        return record

    def _parse_info(self, info_str):
        """
        Remake version of pyVCF.parser._parse_info

        Parse the INFO field of a VCF entry into a dictionary of Python
        types.

        ***** Remake part *****
        - convert list that unintentionally indicated by comma back to string
        - decode hex value coded for reserved character (";", "=", etc.)
        - detect if any 'info' entries are annotated by ANNOVAR more than once
          , if yes convert them into 'list'
        """
        if info_str == '.':
            return {}

        entries = info_str.split(';')
        retdict = OrderedDict()

        for entry in entries:
            entry = entry.split('=', 1)
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if (ID in self.annovar_infos) and (entry_type in ('String', 'Character')):
                val = entry[1]
                # decode hex character
                hex_idx = val.find("\\x")
                while hex_idx > 0:
                    hex_str = val[hex_idx:hex_idx+4]
                    decode_chr = val[hex_idx+2:hex_idx+4].decode('hex')
                    val = val.replace(hex_str, decode_chr)
                    hex_idx = val.find("\\x")
            elif entry_type == 'Integer':
                vals = entry[1].split(',')
                try:
                    val = self._map(int, vals)
                # Allow specified integers to be flexibly parsed as floats.
                # Handles cases with incorrectly specified header types.
                except ValueError:
                    val = self._map(float, vals)
            elif entry_type == 'Float':
                vals = entry[1].split(',')
                val = self._map(float, vals)
            elif entry_type == 'Flag':
                val = True
            elif entry_type in ('String', 'Character'):
                try:
                    vals = entry[1].split(',') # commas are reserved characters indicating multiple values
                    val = self._map(str, vals)
                except IndexError:
                    entry_type = 'Flag'
                    val = True

            try:
                if self.infos[ID].num == 1 and entry_type not in ( 'Flag', ):
                    val = val[0]
            except KeyError:
                pass

            # Check if there are infos of more than one alternate alleles annotated by ANNOVAR
            if (ID in self.annovar_infos) and (ID in retdict) and (type(retdict[ID]) is list):
                retdict[ID].append(val)
                val = retdict[ID]
            elif (ID in self.annovar_infos) and (ID in retdict):
                val = [retdict[ID], val]
            retdict[ID] = val

        return retdict

    def _parse_samples(self, samples, samp_fmt, site):
        '''
        Remake version of pyVCF.parser._parse_samples

        Parse a sample entry according to the format specified in the FORMAT
        column.

        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent

        ***** Remake part *****
        - use custom _Call structure, _TAVcfCall
        '''

        # check whether we already know how to parse this format
        if samp_fmt not in self._format_cache:
            self._format_cache[samp_fmt] = self._parse_sample_format(samp_fmt)
        samp_fmt = self._format_cache[samp_fmt]

        if cparse:
            return cparse.parse_samples(
                self.samples, samples, samp_fmt, samp_fmt._types, samp_fmt._nums, site)

        samp_data = []
        _map = self._map

        nfields = len(samp_fmt._fields)

        for name, sample in itertools.izip(self.samples, samples):

            # parse the data for this sample
            sampdat = [None] * nfields

            for i, vals in enumerate(sample.split(':')):

                # short circuit the most common
                if samp_fmt._fields[i] == 'GT':
                    sampdat[i] = vals
                    continue
                elif vals == ".":
                    sampdat[i] = None
                    continue

                entry_num = samp_fmt._nums[i]
                entry_type = samp_fmt._types[i]

                # we don't need to split single entries
                if entry_num == 1 or ',' not in vals:

                    if entry_type == 'Integer':
                        try:
                            sampdat[i] = int(vals)
                        except ValueError:
                            sampdat[i] = float(vals)
                    elif entry_type == 'Float':
                        sampdat[i] = float(vals)
                    else:
                        sampdat[i] = vals

                    if entry_num != 1:
                        sampdat[i] = (sampdat[i])

                    continue

                vals = vals.split(',')

                if entry_type == 'Integer':
                    try:
                        sampdat[i] = _map(int, vals)
                    except ValueError:
                        sampdat[i] = _map(float, vals)
                elif entry_type == 'Float' or entry_type == 'Numeric':
                    sampdat[i] = _map(float, vals)
                else:
                    sampdat[i] = vals

            # create a call object
            call = _TAVcfCall(site, name, samp_fmt(*sampdat))
            samp_data.append(call)

        return samp_data

class TAVcfInfoReader(TAVcfReader):
    """
    An encapsulated version of TAVcfReader to only read INFO field
    """

    def __init__(self,
                 *args,
                 **kwargs
                 ):
        super(TAVcfInfoReader, self).__init__(*args, **kwargs)
        self.__header_cols = self.__parse_header_cols()

    def __parse_header_cols(self):
        header_cols = []
        header_cols.append('CHROM')
        header_cols.append('POS')
        header_cols.append('REF')
        header_cols.append('ALT')
        annovar_date_idx = self.infos.keys().index("ANNOVAR_DATE")
        anno_fields = self.infos.keys()[annovar_date_idx+1:-1]
        for anno_field in anno_fields:
            header_cols.append("_"+anno_field.replace(".", "_"))
        return header_cols

    @property
    def header_cols(self):
        return self.__header_cols

    @property
    def info_cols(self):
        return self.header_cols[4:]

    @property
    def record_tuples(self):
        for item in super(TAVcfInfoReader, self).__iter__():
            yield item.as_annovar_tuple(self.annovar_infos.keys())

#    def next(self):
#        return super(TAVcfInfoReader, self).next().as_annovar_tuple(self.annovar_infos.keys())

class TAVcfGTZReader(TAVcfReader):
    """
    An encapsulated version of TAVcfReader to only read Call data
    Note: It's not necessary to inherit from TAVcfReader but to make it
    convinient for code maintainance
    """

    def __init__(self,
                 *args,
                 **kwargs
                 ):
        super(TAVcfGTZReader, self).__init__(*args, **kwargs)
        self.__header_cols = self.__parse_header_cols()

    def __parse_header_cols(self):
        header_cols = []
        header_cols.append('CHROM')
        header_cols.append('POS')
        header_cols.append('REF')
        header_cols.append('ALT')
        header_cols.append('QUAL')
        header_cols.append('FILTER')
        for sample_id in self.samples:
            header_cols.append("_"+sample_id.replace("-", "_"))
        return header_cols

    @property
    def header_cols(self):
        return self.__header_cols

    @property
    def samples_id(self):
        return self.header_cols[6:]

    @property
    def record_tuples(self):
        for item in super(TAVcfGTZReader, self).__iter__():
            yield item.as_gtz_tuple(self.samples)
#    def next(self):
#        return super(TAVcfGTZReader, self).next().as_gtz_tuple(self.samples)

    @property
    def coors(self):
        for item in super(TAVcfGTZReader, self).__iter__():
            yield item.as_pkey_tuple()
#        for item in self.next():
#            self.dbg(item)
#            self.dbg(item[0:4])
#            yield item.as_gtz_tuple()
#        return super(TAVcfGTZReader, self).next().as_pkey_tuple()
