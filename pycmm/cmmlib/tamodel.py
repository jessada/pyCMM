import copy
import sys
import re
from collections import defaultdict
from vcf.model import _Record as _VcfRecord
from vcf.model import _Call as _VcfCall
from pycmm.template import pyCMMBase
from pycmm.settings import PRIMARY_MAF_VAR
from pycmm.settings import FUNC_REFGENE_VAR
from pycmm.settings import EXONICFUNC_REFGENE_VAR
from pycmm.settings import EST_KVOT_EARLYONSET_VS_BRC_COL_NAME
from pycmm.settings import EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME
from pycmm.settings import EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME
from pycmm.settings import PATHOGENIC_COUNT_COL_NAME
from pycmm.settings import WES294_OAF_EARLYONSET_AF_COL_NAME
from pycmm.settings import WES294_OAF_BRCS_AF_COL_NAME
from pycmm.settings import EXAC_NFE_COL_NAME
from pycmm.settings import KG2014OCT_EUR_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_SYN_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_SYN_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_SYN_Z_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_MIS_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_MIS_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_MIS_Z_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_LOF_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_LOF_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_PLI_COL_NAME
from pycmm.settings import REF_MAF_COL_NAMES
from pycmm.settings import MAX_REF_MAF_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_COL_NAMES
from pycmm.settings import EXAC03_CONSTRAINT_COL_NAME
from pycmm.settings import PREDICTION_COLS
from pycmm.settings import INTERVAR_AND_EVIDENCE_COL_NAME
from pycmm.settings import INTERVAR_CLASS_COL_NAME
from pycmm.settings import INTERVAR_EVIDENCE_COL_NAME
from pycmm.utils import check_equal
from pycmm.utils import check_in
from pycmm.utils import is_number
from pycmm.cmmlib.intervarlib import parse_intervar_class
from pycmm.cmmlib.intervarlib import parse_intervar_evidence

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

EXAC_CONSTRAINT_PATTERN = re.compile(r'''(?P<var_name>.+?)=(?P<value>.+)''')

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
        for gt_idx in xrange(len(self.cmm_gts)):
            cmm_gt = self.cmm_gts[gt_idx]
            # Other than the problematic wild type and homozygote,
            # the actual gt should be the same as cmm gt
            if cmm_gt != CMMGT_HOMOZYGOTE and cmm_gt != CMMGT_WILDTYPE:
                actual_gts.append(cmm_gt)
                continue
            # if reference is not mutated the actual gt remain the same
            if not self.site.ref_is_mutated[gt_idx]:
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
        if self.__shared_mutations is None:
            self.site.cal_shared()
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
        self.__afss = self.__cal_afss()
        self.__is_intergenic = None
        self.__is_intronic = None
        self.__is_synonymous = None
        self.__is_upstream = None
        self.__is_downstream = None
        self.__is_utr = None
        self.__ref_is_mutated = None
        self.__family_infos = copy.deepcopy(family_infos)
        self.__shared_cal = False
        self.__pathogenic_counts = {}
        self.__exac_constraints = defaultdict(dict)
        self.__max_ref_maf = {}

        if (type(self.FILTER) is list) and (len(self.FILTER) == 0):
            self.FILTER = "PASS"
        elif type(self.FILTER) is list:
            self.FILTER = ";".join(self.FILTER)
        else:
            self.FILTER = "."

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

    @property
    def ref_is_mutated(self):
        if self.__ref_is_mutated is None:
            afs = self.afss[0]
            self.__ref_is_mutated = []
            for gt_idx in xrange(len(afs)):
                if ((afs[gt_idx] == ".") or
                    (afs[gt_idx] == CMM_DUMMY) or
                    (float(afs[gt_idx]) < 0.5)
                    ):
                    self.__ref_is_mutated.append(False)
                else:
                    self.__ref_is_mutated.append(True)
        return self.__ref_is_mutated

    def __get_info(self, var_name):
        """
        for internal call
        - return list of values of table_annovar column "var_name"
        - number of entry in the list = 1 + number of alternate alleles
        """
        if var_name in self.INFO:
            return self.INFO[var_name]
        return None

    def get_info(self, var_name, allele_idx=1):
        """
        parsed to be called by high-level function
        - kvot is included
        """
        def cal_est_ors(cases_freq,
                        ctrls_freq,
                        ref_is_mutated,
                        ):
            # filter out none number
            if cases_freq == "NA":
                return "NA"
            if (cases_freq is None or
                cases_freq == "" or
                not is_number(cases_freq)
                ):
                return ""
            if (ctrls_freq is None or
                ctrls_freq == "" or
                not is_number(ctrls_freq)
                ):
                return ""
            # filter "divide by 0"
            if float(ctrls_freq) == 0:
                return "INF"
            if ref_is_mutated and (float(ctrls_freq) == 1):
                return "INF"
            if ref_is_mutated:
                return "{:.4f}".format((1-float(cases_freq))/(1-float(ctrls_freq)))
            return "{:.4f}".format(float(cases_freq)/float(ctrls_freq))

        info = self.__get_info(var_name)
        if info is not None:
            if (type(info) is list) and (len(info) == 1):
                info = info[0]
            elif (type(info) is list) and (len(info) > 1):
                info = info[allele_idx-1]
        elif var_name == PATHOGENIC_COUNT_COL_NAME:
            info = self.pathogenic_count(allele_idx=allele_idx)
        elif var_name == EST_KVOT_EARLYONSET_VS_BRC_COL_NAME:
            info = cal_est_ors(cases_freq=self.get_info(WES294_OAF_EARLYONSET_AF_COL_NAME,
                                                        allele_idx),
                               ctrls_freq=self.get_info(WES294_OAF_BRCS_AF_COL_NAME,
                                                        allele_idx),
                               ref_is_mutated=self.ref_is_mutated[allele_idx],
                               )
        elif var_name == EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME:
            info = cal_est_ors(cases_freq=self.get_info(WES294_OAF_EARLYONSET_AF_COL_NAME,
                                                        allele_idx),
                               ctrls_freq=self.get_info(EXAC_NFE_COL_NAME,
                                                        allele_idx),
                               ref_is_mutated=self.ref_is_mutated[allele_idx],
                               )
        elif var_name == EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME:
            info = cal_est_ors(cases_freq=self.get_info(WES294_OAF_EARLYONSET_AF_COL_NAME,
                                                        allele_idx),
                               ctrls_freq=self.get_info(KG2014OCT_EUR_COL_NAME,
                                                        allele_idx),
                               ref_is_mutated=self.ref_is_mutated[allele_idx],
                               )
        elif var_name in EXAC03_CONSTRAINT_COL_NAMES:
            info = self.__get_exac_constraint_val(var_name, allele_idx)
        elif var_name == MAX_REF_MAF_COL_NAME:
            info = self.__get_max_ref_maf(allele_idx)
        elif var_name == INTERVAR_CLASS_COL_NAME:
            info = parse_intervar_class(self.get_info(INTERVAR_AND_EVIDENCE_COL_NAME,
                                                      allele_idx))
        elif var_name == INTERVAR_EVIDENCE_COL_NAME:
            info = parse_intervar_evidence(self.get_info(INTERVAR_AND_EVIDENCE_COL_NAME,
                                                         allele_idx))
        if (info == "" or
            info is None or
            info == [None] or
            info == "."
            ):
           info = ""
        return info

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

    def cal_shared(self):
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
            self.cal_shared()
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

    def is_pass_vqsr(self, allele_idx=1):
        return self.FILTER == "PASS"

    def pathogenic_count(self, allele_idx=1):
        if allele_idx not in self.__pathogenic_counts:
            count = 0
            for col_name in PREDICTION_COLS:
                info = self.get_info(col_name, allele_idx)
                if not info.harmful:
                    continue
                count += 1
            self.__pathogenic_counts[allele_idx] = count
        return self.__pathogenic_counts[allele_idx]

    def __get_exac_constraint_val(self, var_name, allele_idx=1):
        if var_name in self.__exac_constraints[allele_idx]:
            return self.__exac_constraints[allele_idx][var_name]
        exac_constraint_vals = self.get_info(EXAC03_CONSTRAINT_COL_NAME, allele_idx)
        if (exac_constraint_vals is None or
            exac_constraint_vals == "." or
            exac_constraint_vals == ""):
            return None
        self.__parse_exac_constraint(exac_constraint_vals, allele_idx)
        return self.__exac_constraints[allele_idx][var_name]

    def __parse_exac_constraint(self, exac_constraint_vals, allele_idx):
        exac_constaints = exac_constraint_vals.split('#')
        for exac_constaint in exac_constaints:
            match = EXAC_CONSTRAINT_PATTERN.match(exac_constaint)
            var_name = match.group('var_name')
            value = match.group('value')
            self.__exac_constraints[allele_idx][var_name] = value

    def __get_max_ref_maf(self, allele_idx=1):
        if allele_idx not in self.__max_ref_maf:
            max_ref_maf = 0
            for ref_maf_col_name in REF_MAF_COL_NAMES:
                ref_maf = self.get_info(ref_maf_col_name, allele_idx)
                if ref_maf == "":
                    continue
                ref_maf = float(ref_maf)
                if ref_maf > 0.5:
                    ref_maf = 1 - ref_maf
                if ref_maf > max_ref_maf:
                    max_ref_maf = ref_maf
            self.__max_ref_maf[allele_idx] = max_ref_maf
        return self.__max_ref_maf[allele_idx]

    def vcf_eval(self, expr, allele_idx):
        def info_repl(match_obj):
            repl_txt = "self.get_info("
            repl_txt += match_obj.group(0)
            repl_txt += ", " + str(allele_idx)
            repl_txt += ")"
            return repl_txt

        # look for annotated fields
        info_fields = {}
        for match in re.finditer(r'(\".+?\")', expr):
            info_fields[match.group(1)] = 1
        # replace the annotated fields with the actual values
        repl_expr = expr
        for info_field in info_fields:
            info_val = eval(re.sub(r'(\".+?\")',
                                   info_repl,
                                   info_field)) 
            repl_expr = re.sub(info_field,
                               "'"+str(info_val)+"'",
                               repl_expr)
        # then eval the expression
        return eval(repl_expr)
