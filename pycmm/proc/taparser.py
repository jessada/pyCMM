import itertools
import re
from vcf import Reader as VcfReader
from vcf import RESERVED_INFO
from collections import OrderedDict
from pycmm.proc.annovarlib import PredictionTranslator
from pycmm.proc.tamodel import _TAVcfRecord
from pycmm.proc.tamodel import _TAVcfCall

try:
    import cparse
except ImportError:
    cparse = None


class TAVcfReader(VcfReader):
    """
    An encapsulated version of vcf.Reader from pyVCF package to
      - enable _parse_info to recognize INFO annotated by annovar
        , <see _parse_info function>
      - redirect _Call to TAVcfCall
      - redirect _Record to TAVcfRecord
    """

    def __init__(self,
                 fsock=None,
                 filename=None,
                 compressed=None,
                 prepend_chr=False,
                 strict_whitespace=False,
                 family_infos=None,
                 ):
        super(TAVcfReader, self).__init__(fsock=fsock,
                                          filename=filename,
                                          compressed=compressed,
                                          prepend_chr=prepend_chr,
                                          strict_whitespace=strict_whitespace,
                                          )
        self.__parse_annovar_infos()
        self.__pred_tran = PredictionTranslator()
        self.__family_infos = family_infos

    @property
    def annovar_infos(self):
        return self.__annovar_infos

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
                info, fmt, self._sample_indexes, family_infos=self.__family_infos)

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

            # Translate in-silico prediction if possible
            if ID in self.__pred_tran.predictor_list:
                val = self.__pred_tran.get_prediction_info(ID, val)

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
