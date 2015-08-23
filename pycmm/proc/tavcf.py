import itertools
import re
from vcf import Reader as VcfReader
from vcf import RESERVED_INFO
from collections import OrderedDict
from pycmm.template import pyCMMBase
from pycmm.utils import mylogger
from pycmm.settings import DFLT_MAF_VAR

try:
    import cparse
except ImportError:
    cparse = None

from vcf.model import _Record as _VcfRecord
from vcf.model import _Call as _VcfCall
from pycmm.settings import LJB_SIFT_PREDICTION_COL_NAME as SIFT_PRED_COL
from pycmm.settings import LJB_POLYPHEN2_HDIV_PREDICTION_COL_NAME as POLYPHEN2_HDIV_PRED_COL
from pycmm.settings import LJB_POLYPHEN2_HVAR_PREDICTION_COL_NAME as POLYPHEN2_HVAR_PRED_COL
from pycmm.settings import LJB_LRT_PREDICTION_COL_NAME as LRT_PRED_COL
from pycmm.settings import LJB_MUTATIONTASTER_PREDICTION_COL_NAME as MUTATIONTASTER_PRED_COL
from pycmm.settings import LJB_MUTATIONASSESSOR_PREDICTION_COL_NAME as MUTATIONASSESSOR_PRED_COL
from pycmm.settings import LJB_FATHMM_PREDICTION_COL_NAME as FATHMM_PRED_COL
from pycmm.settings import LJB_RADIALSVM_PREDICTION_COL_NAME as RADIALSVM_PRED_COL
from pycmm.settings import LJB_LR_PREDICTION_COL_NAME as LR_PRED_COL

CMMGT_WILDTYPE = 'wt'
CMMGT_HOMOZYGOTE = 'hom'
CMMGT_HETEROZYGOTE = 'het'
CMMGT_OTHER = 'oth'
CMMGT_DUMMY = 'dummy'


class PredictionTranslator(pyCMMBase):
    """
    A class to translate codes from effect predictors by using informaiton
    from http://annovar.openbioinformatics.org/en/latest/user-guide/filter/
    As 20 August 2015 the translation is as below

    SIFT (sift)                        D: Deleterious (sift<=0.05); T: tolerated (sift>0.tolerated05)
    PolyPhen 2 HDIV (pp2_hdiv)         D: Probably damaging (>=0.957),damaging P: possibly damaging (0.453<=pp2_hdiv<=0.956); B: benign (pp2_hdiv<=0.452)
    PolyPhen 2 HVar (pp2_hvar)         D: Probably damaging (>=0.909), P: possibly damaging (0.447<=pp2_hdiv<=0.909); B: benign (pp2_hdiv<=0.446)
    LRT (lrt)                          D: Deleterious; N: Neutral; U: Unknown
    MutationTaster (mt)                A" ("disease_causing_automatic"); "D" ("disease_causing"); "N" ("PolyPhenmorphism"); "P" ("polymorphism_automatic")
    MutationAssessor (ma)              H: high; M: medium; L: low; N: neutral. H/M means functional and L/N   means non-functional
    FATHMM (fathmm)                    D: Deleterious; T: Tolerated
    MetaSVM (metasvm)                  D: Deleterious; T: Tolerated
    MetaLR (metasvmalr)                D: Deleterious; T: Tolerated
    GERP++ (gerp++)                    higher scores are more deleterious
    PhyloP (phylop)                    higher scores are more deleterious
    SiPhy (siphy)                      higher scores are more deleterious
    """

    def __init__(self):
        self.__expl = self.__init_explanation()

    def __init_explanation(self):
        expl = {}
        expl[SIFT_PRED_COL] = {'D': 'Deleterious', 'T': 'Tolerated'}
        expl[POLYPHEN2_HDIV_PRED_COL] = {'D': 'Probably damaging', 'P': 'Possibly damaging', 'B': 'Benign'}
        expl[POLYPHEN2_HVAR_PRED_COL] = {'D': 'Probably damaging', 'P': 'Possibly damaging', 'B': 'Benign'}
        expl[LRT_PRED_COL] = {'D': 'Deleterious', 'N': 'Neutral', 'U': 'Unknown'}
        expl[MUTATIONTASTER_PRED_COL] = {'A': 'Disease causing automatic', 'D': 'Disease causing', 'N': 'PolyPhenmorphism', 'P': 'Polymorphism_automatic'}
        expl[MUTATIONASSESSOR_PRED_COL] = {'H': 'High', 'M': 'Medium', 'L': 'Low', 'N': 'Neutral', 'H/M': 'Functional', 'L/N': 'Non-functional'}
        expl[FATHMM_PRED_COL] = {'D': 'Deleterious', 'T': 'Tolerated'}
        expl[RADIALSVM_PRED_COL] = {'D': 'Deleterious', 'T': 'Tolerated'}
        expl[LR_PRED_COL] = {'D': 'Deleterious', 'T': 'Tolerated'}
        return expl

    def recognize(self, col_name):
        return col_name in self.__expl

    def describe(self, col_name, code):
        if code == ".":
            return "."
        if code in self.__expl[col_name]:
            return self.__expl[col_name][code]
        else:
            return "NA"


class _CmmVcfCall(_VcfCall):
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
        self.__cmm_gts = None
        self.__actual_gts = None
        self.__mutated = None

    def __cal_extra_attributes(self):
        self.__cmm_gts = self.__cal_cmm_gts()
#        self.__afs = self.__cal_afs()
        self.__actual_gts = self.__cal_actual_gts()
        self.__mutated = self.__cal_mutated()

    def __cal_cmm_gts(self):
        """
        to calculate type of genotype given allele index
          - 0 -> REF
          - 1 -> first ALT
          - and so on

        ** NOTE ** the calculation here based on GT value alone
        """
        raw_GT = self.data.GT
        cmm_gts = [CMMGT_DUMMY]
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

#    def __cal_afs(self):
#        raw_afs = self.site.INFO[DFLT_MAF_VAR]
#        if (raw_afs == "") or (raw_afs is None) or (raw_afs == "."):
#            return map(lambda x: 0, xrange(len(self.site.alleles)))
#        if type(raw_afs) is not list:
#            afs = [CMMGT_DUMMY, raw_afs]
#        else:
#            afs = [CMMGT_DUMMY]
#            for raw_af in raw_afs:
#                if raw_af is None:
#                    afs.append(0)
#                else:
#                    afs.append(raw_af)
#        return afs

    def __cal_actual_gts(self):
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
        actual_gts = []
        afs = self.site.afs
        for gt_idx in xrange(len(self.cmm_gts)):
            cmm_gt = self.cmm_gts[gt_idx]
            # Other than the problematic wild type and homozygote,
            # the actual gt should be the same as cmm gt
            if cmm_gt != CMMGT_HOMOZYGOTE and cmm_gt != CMMGT_WILDTYPE:
                actual_gts.append(cmm_gt)
                continue
            # if allele frequency less than 0.5 the actual gt remain the same
            if afs[gt_idx] < 0.5:
                actual_gts.append(cmm_gt)
                continue
            # below should be only hom and wild type with allele freq >= 0.5
            if cmm_gt == CMMGT_HOMOZYGOTE:
                actual_gts.append(CMMGT_WILDTYPE)
                continue
            actual_gts.append(CMMGT_HOMOZYGOTE)
        return actual_gts

    def __cal_mutated(self):
        """
        to identify if a genotype is actually mutated given
          - allele index
            - 0 -> REF
            - 1 -> first ALT
            - and so on
        """
        mutated = [CMMGT_DUMMY]
        for gt_idx in xrange(1, len(self.actual_gts)):
            if self.actual_gts[gt_idx] == CMMGT_HOMOZYGOTE:
                mutated.append(True)
                continue
            if self.actual_gts[gt_idx] == CMMGT_HETEROZYGOTE:
                mutated.append(True)
                continue
            mutated.append(False)
        return mutated

    @property
    def cmm_gts(self):
        if self.__cmm_gts is None:
            self.__cal_extra_attributes()
        return self.__cmm_gts

    @property
    def actual_gts(self):
        if self.__actual_gts is None:
            self.__cal_extra_attributes()
        return self.__actual_gts

    @property
    def mutated(self):
        if self.__mutated is None:
            self.__cal_extra_attributes()
        return self.__mutated

class _CmmVcfRecord(_VcfRecord):
    """
    An encapsulated version of vcf._Record from pyVCF package to
      - determine if mutations are shared between samples
      - understand if itself is a rare mutation given frequency ratio
    """

    def __init__(self, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
            sample_indexes, samples=None):
        _VcfRecord.__init__(self,
                            CHROM,
                            POS,
                            ID,
                            REF,
                            ALT,
                            QUAL,
                            FILTER,
                            INFO,
                            FORMAT,
                            sample_indexes,
                            samples=None,
                            )
        self.__afs = self.__cal_afs()

    def __cal_afs(self):
        """
        return list of allele frequencies.
        Number of entry in the list = 1 + number of alternate alleles
        """
        if DFLT_MAF_VAR in self.INFO:
            raw_afs = self.INFO[DFLT_MAF_VAR]
        else:
            raw_afs = None
        if (raw_afs == "") or (raw_afs is None) or (raw_afs == "."):
            return map(lambda x: 0, xrange(len(self.alleles)))
        if type(raw_afs) is not list:
            afs = [CMMGT_DUMMY, raw_afs]
        else:
            afs = [CMMGT_DUMMY]
            for raw_af in raw_afs:
                if raw_af is None:
                    afs.append(0)
                else:
                    afs.append(raw_af)
        return afs

    @property
    def afs(self):
        return self.__afs

    def is_shared(self, allele_idx, samples):
        """
        to identify if a genotype is shared betwen samples given
          - allele index
            - 0 -> REF
            - 1 -> first ALT
            - and so on
           - samples, list of samples to be compared if the list is
             empty, None, it will return True
        """
        if samples is None:
            return True
        if len(samples) == 0:
            return True
        for sample in samples:
            if not self.genotype(sample).mutated[allele_idx]:
                return False
        return True

    def is_rare(self, freq_ratios, allele_idx=1):
        condition, criteria = freq_ratios.items()[0]
        criteria = float(criteria)
        val = self.afs[allele_idx]
        if (val is None or
            val == "." or
            val == ""):
            return True
        mylogger.debug(condition)
        mylogger.debug(type(self.afs))
        mylogger.debug(criteria)
        freq = float(val)
        if freq < float(criteria):
            return True
        if freq > (1-float(criteria)):
            return True
        return False

class TableAnnovarVcfReader(VcfReader):
    """
    An encapsulated version of vcf.Reader from pyVCF package to
      - enable _parse_info to recognize INFO annotated by annovar
        , <see _parse_info function>
      - redirect _Call to CmmVcfCall
      - redirect _Record to CmmVcfRecord
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
        self.__pred_tran = PredictionTranslator()

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
        - use custom _Record structure, _CmmVcfRecord
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

        record = _CmmVcfRecord(chrom, pos, ID, ref, alt, qual, filt,
                info, fmt, self._sample_indexes)

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
            if self.__pred_tran.recognize(ID):
                val = self.__pred_tran.describe(ID, val)

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
