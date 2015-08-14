import itertools
from vcf import Reader as VcfReader
from vcf import RESERVED_INFO
from collections import OrderedDict
from pycmm.utils import mylogger
from pycmm.settings import DFLT_MAF_VAR

try:
    import cparse
except ImportError:
    cparse = None

from vcf.model import _Call as _VCFCall

CMMGT_WILDTYPE = 'wt'
CMMGT_HOMOZYGOTE = 'hom'
CMMGT_HETEROZYGOTE = 'het'
CMMGT_OTHER = 'oth'


class _CmmGt(object):
    """
    to calculate type of genotype given allele index
      - 0 -> REF
      - 1 -> first ALT
      - and so on

    ** NOTE ** the calculation here based on GT value alone
    """

    def __init__(self, gt):
        self.__gt = gt

    def __repr__(self):
        if (self.__gt == ".") or (self.__gt == "./."):
            return "."
        return str(map(lambda x: self[x], range(len(self)+1)))

    def __len__(self):
        if (self.__gt == ".") or (self.__gt == "./."):
            return 0
        return max(map(lambda x: int(x),
                       self.__gt.split("/")))

    def __getitem__(self, allele_idx):
        if (self.__gt == ".") or (self.__gt == "./."):
            return "."
        gts = self.__gt.split("/")
        # 0/0 will be translated as "wild type" no matter what allele index is
        if (gts[0] == "0") and (gts[1] == "0"):
            return CMMGT_WILDTYPE
        # The gt of which the allele index doesn't match will be considered
        # as "other", ex. allele index = 1 but the gt is 2/3
        if (gts[0] != str(allele_idx)) and (gts[1] != str(allele_idx)):
            return CMMGT_OTHER
        # The rest are in the cases one or both of the alleles match with
        # allele index.
        # So if they are different then "heterozygote" else "homozygote"
        if gts[0] != gts[1]:
            return CMMGT_HETEROZYGOTE
        return CMMGT_HOMOZYGOTE

class _ActualGt(_CmmGt):
    """
    to calculate the actual genotype based allele index and allele frequency
      - allele index
        - 0 -> REF
        - 1 -> first ALT
        - and so on

    ** NOTE1 ** the calculation here based on GT value alone
    ** NOTE2 ** the calculation may not be accurate if there are more than
    one alternate alleles, Ex Ref=A (freq=45%), Alt=T,C,G(freq=10%,10%,35%).
    But it'll be a very rare case
    """

    def __init__(self, gt, af=None):
        self.__gt = gt
        if (af == "") or (af is None) or (af == "."):
            self.__af = 0
        else:
            self.__af = af

    def __getitem__(self, allele_idx):
        def get_af(allele_idx):
            if type(self.__af) is not list:
                return self.__af
            af = self.__af[allele_idx-1]
            if af is None:
                return 0
            return af

        if (self.__gt == ".") or (self.__gt == "./."):
            return "."
        gts = self.__gt.split("/")
        # 0/0 will be translated as "wild type" if allele frequency of that
        # allele index is less than 0.5. Otherwise, it'll be considered as the
        # opposite, which is homozygote mutation
        if (gts[0] == "0") and (gts[1] == "0") and (get_af(allele_idx) < 0.5):
            return CMMGT_WILDTYPE
        if (gts[0] == "0") and (gts[1] == "0") and (get_af(allele_idx) >= 0.5):
            return CMMGT_HOMOZYGOTE
        # Regardless of frequency,
        # the gt of which the allele index doesn't match will be considered
        # as "other", ex. allele index = 1 but the gt is 2/3
        if (gts[0] != str(allele_idx)) and (gts[1] != str(allele_idx)):
            return CMMGT_OTHER
        # The rest are in the cases one or both of the alleles match with
        # allele index.
        # So if they are different then "heterozygote" else "homozygote"
        if gts[0] != gts[1]:
            return CMMGT_HETEROZYGOTE
        # Otherwise, if the frequency less than 0.5 it'll consider as
        # "homozygote"
        if get_af(allele_idx) < 0.5:
            return CMMGT_HOMOZYGOTE
        # else "wild type"
        return CMMGT_WILDTYPE

class _CmmVcfCall(_VCFCall):
    """
    An encapsulated version of vcf._Call from pyVCF package to
      - add extra genotype translation, "cmm_gt", to handle record
        with more than one alternate alleles
    """

    @property
    def cmm_gt(self):
        return _CmmGt(self.data.GT)

    @property
    def actual_gt(self):
        return _ActualGt(self.data.GT, self.site.INFO[DFLT_MAF_VAR])


class TableAnnovarVcfReader(VcfReader):
    """
    An encapsulated version of vcf.Reader from pyVCF package to fix
      - _parse_info  <see _parse_info function>
    """

    def __init__(self,
                 fsock=None,
                 filename=None,
                 compressed=None,
                 prepend_chr=False,
                 strict_whitespace=False,
                 ):
        super(TableAnnovarVcfReader, self).__init__(fsock=fsock,
                                                    filename=filename,
                                                    compressed=compressed,
                                                    prepend_chr=prepend_chr,
                                                    strict_whitespace=strict_whitespace,
                                                    )
        self.__parse_annovar_infos()

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
        - use custom _Call structure, _CmmVcfCall
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
            call = _CmmVcfCall(site, name, samp_fmt(*sampdat))
            samp_data.append(call)

        return samp_data


