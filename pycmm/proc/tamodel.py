import copy
import sys
import re
from vcf.model import _Record as _VcfRecord
from vcf.model import _Call as _VcfCall
from pycmm.template import pyCMMBase
from pycmm.settings import PRIMARY_MAF_VAR
from pycmm.settings import FUNC_REFGENE_VAR
from pycmm.settings import EXONICFUNC_REFGENE_VAR
from pycmm.utils import check_equal
from pycmm.utils import check_in

CMMGT_WILDTYPE = 'wt'
CMMGT_HOMOZYGOTE = 'hom'
CMMGT_HETEROZYGOTE = 'het'
CMMGT_OTHER = 'oth'
CMM_DUMMY = 'dummy'

FUNC_INTERGENIC = 'intergenic'
FUNC_INTRONIC = 'intronic'
FUNC_UPSTREAM = 'upstream'
FUNC_DOWNSTREAM = 'downstream'
FUNC_UTR = 'UTR'
EXONICFUNC_SYNONYMOUS = 'synonymous_SNV'

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
        self.__actual_gts = None
        self.__mutated = None
        self.__shared_mutations = None

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
        afs = self.site.afss[0]
        for gt_idx in xrange(len(self.cmm_gts)):
            cmm_gt = self.cmm_gts[gt_idx]
            # Other than the problematic wild type and homozygote,
            # the actual gt should be the same as cmm gt
            if cmm_gt != CMMGT_HOMOZYGOTE and cmm_gt != CMMGT_WILDTYPE:
                actual_gts.append(cmm_gt)
                continue
            # if allele frequency less than 0.5 the actual gt remain the same
            if (afs[gt_idx] == ".") or (float(afs[gt_idx]) < 0.5):
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
        mutated = [CMM_DUMMY]
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

    @property
    def shared_mutations(self):
        return self.__shared_mutations

    @shared_mutations.setter
    def shared_mutations(self, val):
        self.__shared_mutations = val

class _TAVcfRecord(_VcfRecord, pyCMMBase):
    """
    An encapsulated version of vcf._Record from pyVCF package to
      - determine if mutations are shared between samples
      - understand if itself is a rare mutation given frequency ratio
      - understand if itself is an intergenic mutation
      - understand if itself is an intronic mutation
    """

    def __init__(self, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
            sample_indexes,
            samples=None,
            family_infos=None,
            freq_ratios=None,
            expressions=None,
            ):
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
        pyCMMBase.__init__(self)
        self.__freq_ratios = freq_ratios
        if expressions is not None:
            self.__exprs = {}
            for expr in expressions.split(";"):
                key, val = expr.split(":")
                self.__exprs[key.strip()] = val
        else:
            self.__exprs = None
        self.__afss = self.__cal_afss()
        self.__is_intergenic = None
        self.__is_intronic = None
        self.__is_synonymous = None
        self.__is_upstream = None
        self.__is_downstream = None
        self.__is_utr = None
        self.__family_infos = copy.deepcopy(family_infos)
        self.__shared_cal = False

    @property
    def freq_ratios(self):
        return self.__freq_ratios

    @property
    def afss(self):
        return self.__afss

    @property
    def is_intergenic(self):
        if self.__is_intergenic is None:
            self.__is_intergenic = self.__compare_infos(FUNC_REFGENE_VAR,
                                                        FUNC_INTERGENIC)
        return self.__is_intergenic

    @property
    def is_intronic(self):
        if self.__is_intronic is None:
            self.__is_intronic = self.__compare_infos(FUNC_REFGENE_VAR,
                                                      FUNC_INTRONIC)
        return self.__is_intronic

    @property
    def is_upstream(self):
        if self.__is_upstream is None:
            self.__is_upstream = self.__compare_infos(FUNC_REFGENE_VAR,
                                                      FUNC_UPSTREAM)
        return self.__is_upstream

    @property
    def is_downstream(self):
        if self.__is_downstream is None:
            self.__is_downstream = self.__compare_infos(FUNC_REFGENE_VAR,
                                                      FUNC_DOWNSTREAM)
        return self.__is_downstream

    @property
    def is_utr(self):
        if self.__is_utr is None:
            self.__is_utr = self.__compare_infos(FUNC_REFGENE_VAR,
                                                      FUNC_UTR)
        return self.__is_utr

    @property
    def is_synonymous(self):
        if self.__is_synonymous is None:
            self.__is_synonymous = self.__compare_infos(EXONICFUNC_REFGENE_VAR,
                                                        EXONICFUNC_SYNONYMOUS,
                                                        compare=check_equal)
        return self.__is_synonymous

    def __get_info(self, var_name):
        """
        - return list of values of table_annovar column "var_name"
        - number of entry in the list = 1 + number of alternate alleles
        """
        if var_name in self.INFO:
            return self.INFO[var_name]
        return None

    def __cal_afss(self):
        """
        - return list of allele frequencies.
        - number of entry in the list = 1 + number of alternate alleles
        """
        if self.freq_ratios is None:
            var_names = [PRIMARY_MAF_VAR]
        else:
            var_names = self.freq_ratios.keys()
        afss = []
        for var_name in var_names:
            if var_name in self.INFO:
                raw_afs = self.INFO[var_name]
            else:
                raw_afs = None
            if (raw_afs == "") or (raw_afs is None) or (raw_afs == "."):
                afs = map(lambda x: 0, xrange(len(self.alleles)))
            elif type(raw_afs) is not list:
                afs = [CMM_DUMMY, raw_afs]
            else:
                afs = [CMM_DUMMY]
                for raw_af in raw_afs:
                    if raw_af is None:
                        afs.append(0)
                    else:
                        afs.append(raw_af)
            afss.append(afs)
        return afss

    def __cal_shared(self):
        """
        - with given family informations, this function will identify shared
        mutation for each family and for each genotype
        """
        if self.__family_infos is None:
            return
        for fam_id in self.__family_infos:
            family_info = self.__family_infos[fam_id]
            shared_mutations = [CMM_DUMMY]
            for allele_idx in xrange(1, len(self.alleles)):
                shared_mutation = True
                for member in family_info.members:
                    if not self.genotype(member.sample_id).mutated[allele_idx]:
                        shared_mutation = False
                        break
                shared_mutations.append(shared_mutation)
            for member in family_info.members:
                self.genotype(member.sample_id).shared_mutations = shared_mutations
            family_info.shared_mutations = shared_mutations

    def __compare_infos(self, var_name, val, compare=check_in):
        """
        - return list of booleans identified if each INFO[var_name] == val
        - number of entry in the list = 1 + number of alternate alleles
        """
        raw_results = self.__get_info(var_name)
        if (raw_results == "") or (raw_results is None):
            return map(lambda x: True, xrange(len(self.alleles)))
        if raw_results == ".":
            return map(lambda x: False, xrange(len(self.alleles)))
        if type(raw_results) is not list:
            raw_results = [raw_results]
        results = [CMM_DUMMY]
        for entry in raw_results:
            if compare(val, entry):
                results.append(True)
            else:
                results.append(False)
        return results

    def has_mutation(self, samples, allele_idx):
        """
        to identify if a genotype mutation is found in any samples given
          - allele index
            - 0 -> REF
            - 1 -> first ALT
            - and so on
           - samples, list of samples to be compared if the list is
             empty, None, it will return True
        """
        if samples is None:
            return False
        if len(samples) == 0:
            return False
        for sample in samples:
            if self.genotype(sample).mutated[allele_idx]:
                return True
        return False

    def is_shared(self, samples, allele_idx, min_share_count=1):
        """
        to identify if a genotype mutation is shared between samples given
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
        share_count = 0
        for sample in samples:
            if not self.genotype(sample).mutated[allele_idx]:
                return False
            share_count += 1
        if share_count >= min_share_count:
            return True
        else:
            return False

    def has_shared(self, allele_idx, min_share_count=1):
        """
        to identify if a genotype mutation has been shared in a family within the given samples
          - allele index
            - 0 -> REF
            - 1 -> first ALT
            - and so on
           - samples, list of samples to be compared if the list is
             empty, None, it will return False
        """
        if not self.__shared_cal:
            self.__cal_shared()
            self.__shared_cal = True
        for fam_id in self.__family_infos:
            family_info = self.__family_infos[fam_id]
            if len(family_info.members) < min_share_count:
                continue
            if family_info.shared_mutations[allele_idx]:
                return True
        return False

    def is_rare(self, allele_idx=1):
        var_names = self.freq_ratios.keys()
        for var_idx in xrange(len(var_names)):
            criteria = self.freq_ratios[var_names[var_idx]]
            val = self.afss[var_idx][allele_idx]
            if (val is None or
                val == "." or
                val == ""):
                continue
            freq = float(val)
            if freq < float(criteria):
                continue
            if freq > (1-float(criteria)):
                continue
            return False
        return True

    def __info_repl(self, match_obj):
        repl_txt = "self.INFO["
        repl_txt += match_obj.group(0)
        repl_txt += "]"
        return repl_txt
    
    def vcf_eval(self, expr_name):
        return eval(re.sub(r'(\".+?\")',
                           self.__info_repl,
                           self.__exprs[expr_name]))
