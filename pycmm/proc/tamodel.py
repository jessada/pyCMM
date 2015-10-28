from pycmm.settings import DFLT_MAF_VAR
from pycmm.settings import FUNC_REFGENE_VAR
from vcf.model import _Record as _VcfRecord
from vcf.model import _Call as _VcfCall

CMMGT_WILDTYPE = 'wt'
CMMGT_HOMOZYGOTE = 'hom'
CMMGT_HETEROZYGOTE = 'het'
CMMGT_OTHER = 'oth'
CMMGT_DUMMY = 'dummy'


class _TAVcfCall(_VcfCall):
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

class _TAVcfRecord(_VcfRecord):
    """
    An encapsulated version of vcf._Record from pyVCF package to
      - determine if mutations are shared between samples
      - understand if itself is a rare mutation given frequency ratio
      - understand if itself is an intergenic mutation
      - understand if itself is an intronic mutation
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
        self.__is_intergenic = self.__compare_infos(FUNC_REFGENE_VAR, 'intergenic')
        self.__is_intronic = self.__compare_infos(FUNC_REFGENE_VAR, 'intronic')

    @property
    def afs(self):
        return self.__afs

    @property
    def is_intergenic(self):
        return self.__is_intergenic

    @property
    def is_intronic(self):
        return self.__is_intronic

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

    def __compare_infos(self, var_name, val):
        """
        return list of booleans identified if each INFO[var_name] == val
        Number of entry in the list = 1 + number of alternate alleles
        """
        if var_name in self.INFO:
            raw_results = self.INFO[var_name]
        else:
            raw_results = None
        if (raw_results == "") or (raw_results is None) or (raw_results == "."):
            return map(lambda x: True, xrange(len(self.alleles)))
        if type(raw_results) is not list:
            raw_results = [raw_results]
        results = [CMMGT_DUMMY]
        for entry in raw_results:
            if entry == val:
                results.append(True)
            else:
                results.append(False)
        return results

    def is_shared(self, samples, allele_idx):
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
        freq = float(val)
        if freq < float(criteria):
            return True
        if freq > (1-float(criteria)):
            return True
        return False
