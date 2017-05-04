# This Python file uses the following encoding: utf-8
import pyaml
from os.path import join as join_path
from collections import defaultdict
from collections import OrderedDict
from collections import namedtuple
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.settings import EST_KVOT_COLS
from pycmm.settings import DFLT_MUTREP_FREQ_RATIOS
from pycmm.settings import PREDICTION_COLS
from pycmm.settings import MUTREP_FAMILY_REPORT_BIN
from pycmm.settings import MUTREP_SUMMARY_REPORT_BIN
from pycmm.settings import EST_KVOT_EARLYONSET_VS_BRC_COL_NAME
from pycmm.settings import EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME
from pycmm.settings import EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME
#from pycmm.settings import EST_KVOT_EARLYONSET_VS_SWEGEN_COL_NAME
from pycmm.settings import WES294_OAF_EARLYONSET_AF_COL_NAME
from pycmm.settings import WES294_OAF_BRCS_AF_COL_NAME
from pycmm.settings import EXAC_NFE_COL_NAME
from pycmm.settings import KG2014OCT_EUR_COL_NAME
from pycmm.settings import GENE_REFGENE_COL_NAME
from pycmm.settings import PATHOGENIC_COUNT_COL_NAME
from pycmm.settings import INTERVAR_AND_EVIDENCE_COL_NAME
from pycmm.settings import INTERVAR_CLASS_COL_NAME
from pycmm.settings import INTERVAR_EVIDENCE_COL_NAME
from pycmm.settings import DFLT_HEADER_CORRECTIONS
from pycmm.settings import EXAC03_CONSTRAINT_COL_NAMES
from pycmm.settings import EXAC03_CONSTRAINT_COL_NAME
from pycmm.settings import FORMAT_COLS
from pycmm.settings import FORMAT_COL_FLOAT
from pycmm.settings import FORMAT_COL_INT
from pycmm.template import pyCMMBase
from pycmm.utils import DefaultOrderedDict
from pycmm.utils import is_number
from pycmm.utils.ver import VersionManager
from pycmm.cmmlib import CMMParams
from pycmm.cmmlib.samplelib import SamplesGroup
from pycmm.cmmlib.taparser import TAVcfReader as VcfReader
from pycmm.cmmlib.xlslib import CMMWorkbook as Workbook
from pycmm.cmmlib.xlslib import NO_COLOR
from pycmm.flow import get_func_arg
from pycmm.flow import init_jobs_setup_file
from pycmm.flow.cmmdb import CMMPipeline
from pycmm.cmmlib.dnalib import ALL_CHROMS
from pycmm.cmmlib.dnalib import DNARegion
from pycmm.cmmlib.tamodel import CMMGT_HOMOZYGOTE
from pycmm.cmmlib.tamodel import CMMGT_HETEROZYGOTE

ACTION_DELETE_ROW = "DELETE_ROW"
ACTION_COLOR_ROW = "COLOR_ROW"
ACTION_COLOR_COL = "COLOR_COLUMN"

DFLT_COLOR_HET_SHARED = 'SILVER'
DFLT_COLOR_HOM_SHARED = 'GRAY25'
CELL_TYPE_HET_SHARED = 'HET_SHARED'
CELL_TYPE_HOM_SHARED = 'HOM_SHARED'
DFLT_COLOR_HET_ZYGO = 'XLS_LIGHT_GREEN'
DFLT_COLOR_HOM_ZYGO = 'XLS_DARK_GREEN'
CELL_TYPE_HET_ZYGO = 'HET_ZYGO'
CELL_TYPE_HOM_ZYGO = 'HOM_ZYGO'
DFLT_COLOR_HET_RECESSIVE = 'XLS_GREEN'
DFLT_COLOR_HOM_RECESSIVE = 'XLS_ORANGE'
CELL_TYPE_HET_RECESSIVE = 'HET_RECESSIVE'
CELL_TYPE_HOM_RECESSIVE = 'HOM_RECESSIVE'
DFLT_COLOR_SEPARATOR = 'XLS_GRAY1_2'
CELL_TYPE_SEPARATOR = 'SEPARATOR'

RECORDS_LOG_INTERVAL = 1000

# *************** report layout section ***************
JOBS_SETUP_RPT_LAYOUT_SECTION = "REPORT_LAYOUT"
JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX = "ANNOTATED_VCF_TABIX"
JOBS_SETUP_RPT_ANNO_COLS_KEY = "COLUMNS"
JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY = "ANNOTATION_EXCLUSION_TAGS"
JOBS_SETUP_RPT_HEADER_CORRECTIONS_KEY = "HEADER_CORRECTIONS"
JOBS_SETUP_RPT_OLD_HEADER_KEY = "OLD_HEADER"
JOBS_SETUP_RPT_NEW_HEADER_KEY = "NEW_HEADER"
JOBS_SETUP_RPT_REGIONS_KEY = "REGIONS"
JOBS_SETUP_RPT_FREQ_RATIOS_KEY = "FREQUENCY_RATIOS"
JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY = "COLUMN"
JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY = "FREQUENCY"
JOBS_SETUP_RPT_EXPRESSIONS_KEY = "EXPRESSIONS"
JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY = "NAME"
JOBS_SETUP_RPT_EXPRESSIONS_PATTERN_KEY = "PATTERN"
JOBS_SETUP_RPT_EXPRESSIONS_USAGES_KEY = "USAGES"
JOBS_SETUP_RPT_EXPRESSIONS_ACTION_KEY = "ACTION"
JOBS_SETUP_RPT_EXPRESSIONS_INFO_KEY = "INFO"
JOBS_SETUP_RPT_SPLIT_CHROM_KEY = "SPLIT_CHROM"
JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY = "SUMMARY_FAMILIES"
JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY = "EXTRA_ANNOTATION_COLUMNS"
JOBS_SETUP_RPT_EXTRA_ANNO_CALL_DETAIL = "Calling_detail"
JOBS_SETUP_RPT_EXTRA_ANNO_CALL_GQ = "Calling_GQ"
JOBS_SETUP_RPT_SHOW_SHARED_MUTATIONS = "Show_shared_mutations"
JOBS_SETUP_RPT_MT_KEY = "Mitochondria"
JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY = "ROWS_FILTER_ACTIONS_CRITERIA"
JOBS_SETUP_RPT_FILTER_RARE = "Rare"
JOBS_SETUP_RPT_FILTER_PASS_VQSR = "PASS-VQSR"
JOBS_SETUP_RPT_FILTER_NON_INTERGENIC = "Non-Intergenic"
JOBS_SETUP_RPT_FILTER_NON_INTRONIC = "Non-Intronic"
JOBS_SETUP_RPT_FILTER_NON_UPSTREAM = "Non-Upstream"
JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM = "Non-Downstream"
JOBS_SETUP_RPT_FILTER_NON_UTR = "Non-UTR"
JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS = "Non-Synonymous"
JOBS_SETUP_RPT_FILTER_HAS_MUTATION = "Has-Mutation"
JOBS_SETUP_RPT_FILTER_HAS_SHARED = "Has-Shared"
JOBS_SETUP_RPT_FILTER_NON_RECESSIVE_GENE = "Non-Recessive-Gene"
#JOBS_SETUP_RPT_ONLY_SUMMARY_KEY = "ONLY_SUMMARY"
#JOBS_SETUP_RPT_ONLY_FAMILIES_KEY = "ONLY_FAMILIES"
JOBS_SETUP_RPT_FILTER_GENES_KEY = "FILTER_GENES"
JOBS_SETUP_RPT_COLORING_SAMPLES_KEY = "COLORING_SAMPLES"
JOBS_SETUP_RPT_COLORING_SHARED = "Shared"
JOBS_SETUP_RPT_COLORING_ZYGOSITY = "Zygosity"

RPT_LAYOUT_CAPTION_ANNOATED_COLS = "annotated columns"
RPT_LAYOUT_CAPTION_FREQ_RATIOS = "frequency ratio(s)"
RPT_LAYOUT_CAPTION_RPT_REGIONS = "report region(s)"
RPT_LAYOUT_CAPTION_FILTER_ACTIONS = "filter actions"

FILTER_RARE = JOBS_SETUP_RPT_FILTER_RARE
FILTER_PASS_VQSR = JOBS_SETUP_RPT_FILTER_PASS_VQSR
FILTER_NON_INTERGENIC = JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
FILTER_NON_INTRONIC = JOBS_SETUP_RPT_FILTER_NON_INTRONIC
FILTER_NON_UPSTREAM = JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
FILTER_NON_DOWNSTREAM = JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
FILTER_NON_UTR = JOBS_SETUP_RPT_FILTER_NON_UTR
FILTER_NON_SYNONYMOUS = JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
FILTER_HAS_MUTATION = JOBS_SETUP_RPT_FILTER_HAS_MUTATION
FILTER_HAS_SHARED = JOBS_SETUP_RPT_FILTER_HAS_SHARED
FILTER_NON_RECESSIVE_GENE = JOBS_SETUP_RPT_FILTER_NON_RECESSIVE_GENE

XLS_CHROM_COL_IDX = 0
XLS_POS_COL_IDX = 1
XLS_REF_COL_IDX = 2
XLS_ALT_COL_IDX = 3
XLS_FILTER_COL_IDX = 4
LAYOUT_VCF_COLS = XLS_FILTER_COL_IDX + 1

RECESSIVE_CPD_HET_CASES_COUNT_COL_IDX = 0
RECESSIVE_CPD_HET_CTRLS_COUNT_COL_IDX = 1
RECESSIVE_CPD_HET_FREQ_RATIO_COL_IDX = 2
RECESSIVE_HOM_CASES_COUNT_COL_IDX = 3
RECESSIVE_HOM_CTRLS_COUNT_COL_IDX = 4
RECESSIVE_HOM_FREQ_RATIO_COL_IDX = 5
RECESSIVE_LAST_COL_IDX = RECESSIVE_HOM_FREQ_RATIO_COL_IDX

CPD_HET_COUNT = 'cpd_het_count'
HOM_COUNT = 'hom_count'
RECESSIVE_SAMPLES_ID = 'recessive_samples_id'

TXT_COMPOUND_HETEROZYGOTE_CASES_COUNT = "compound heterozygote cases count"
TXT_COMPOUND_HETEROZYGOTE_FREQ_RATIO = "compound heterozygote frequency ratio"
TXT_HOMOZYGOTE_CASES_COUNT = "homozygote cases count"
TXT_HOMOZYGOTE_FREQ_RATIO = "homozygote frequency ratio"

VcfRecordBuffer = namedtuple('VcfRecordBuffer',
                             'vcf_record allele_idx')
RecessiveAnalysisResult = namedtuple('RecessiveAnalysisResult',
                                     'affected unaffected filtered_buffer_idxs')

class ActionDelRow(pyCMMBase):
    """ A structure to parse action to delete a row from a mutation report sheet """

    def __init__(self,
                 pattern,
                 *args,
                 **kwargs
                 ):
        super(ActionDelRow, self).__init__(*args, **kwargs)
        self.__pattern = pattern

    def get_raw_obj_str(self):
        return {"pattern": self.pattern}

    @property
    def pattern(self):
        return self.__pattern

class ActionColorRow(pyCMMBase):
    """ A structure to parse action to color a row in a mutation report sheet """

    def __init__(self,
                 pattern,
                 info,
                 *args,
                 **kwargs
                 ):
        super(ActionColorRow, self).__init__(*args, **kwargs)
        self.__pattern = pattern
        self.__info = info

    def get_raw_obj_str(self):
        return {"pattern": self.pattern,
                "color": self.color}

    @property
    def pattern(self):
        return self.__pattern

    @property
    def color(self):
        return self.__info

class ActionColorCol(pyCMMBase):
    """ A structure to parse action to color a column in a row in a mutation report sheet """

    def __init__(self,
                 pattern,
                 info,
                 *args,
                 **kwargs
                 ):
        super(ActionColorCol, self).__init__(*args, **kwargs)
        self.__pattern = pattern
        self.__info = info

    def get_raw_obj_str(self):
        return {"pattern": self.pattern,
                "column name": self.col_name,
                "color": self.color}

    @property
    def pattern(self):
        return self.__pattern

    @property
    def col_name(self):
        return self.__info.split(":")[0]

    @property
    def color(self):
        return self.__info.split(":")[1]

class VcfExpressions(pyCMMBase):
    """
    An user-friendly structure encapsulating
    VCF expression configuration
    """

    def __init__(self,
                 expressions,
                 *args,
                 **kwargs
                 ):
        super(VcfExpressions, self).__init__(*args, **kwargs)
        self.__expressions = expressions
        self.__patterns = None
        self.__actions = None

    @property
    def name(self):
        return self.__expression[JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY]

    @property
    def patterns(self):
        if self.__patterns is None:
            self.__patterns = {}
            for expr in self.__expressions:
                name = expr[JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY]
                pattern = expr[JOBS_SETUP_RPT_EXPRESSIONS_PATTERN_KEY]
                self.__patterns[name] = pattern
        return self.__patterns

    @property
    def actions(self):
        if self.__actions is None:
            # all kind of possible actions are parsed here
            # possible actions are
            #   - delete row
            #   - color a column at the row
            #   - color a row
            # each action is structured as a list
            # each item in the list consist of at least a pattern
            # -> if the evaluation of the pattern is True the do the action
            # and the item can also have an info field
            self.__actions = defaultdict(list)
            actions_color_col = defaultdict(list)
            for expr in self.__expressions:
                if JOBS_SETUP_RPT_EXPRESSIONS_USAGES_KEY in expr:
                    pattern = expr[JOBS_SETUP_RPT_EXPRESSIONS_PATTERN_KEY]
                    for usage in expr[JOBS_SETUP_RPT_EXPRESSIONS_USAGES_KEY]:
                        action = usage[JOBS_SETUP_RPT_EXPRESSIONS_ACTION_KEY]
                        if JOBS_SETUP_RPT_EXPRESSIONS_INFO_KEY in usage:
                            info = usage[JOBS_SETUP_RPT_EXPRESSIONS_INFO_KEY]
                        if action == ACTION_DELETE_ROW:
                            self.__actions[ACTION_DELETE_ROW].append(ActionDelRow(pattern))
                        if action == ACTION_COLOR_ROW:
                            self.__actions[ACTION_COLOR_ROW].append(ActionColorRow(pattern, info=info))
                        if action == ACTION_COLOR_COL:
                            acc = ActionColorCol(pattern, info=info)
                            actions_color_col[acc.col_name].append(acc)
            self.__actions[ACTION_COLOR_COL] = actions_color_col
        return self.__actions

class ReportLayout(CMMParams):
    """ A structure to parse and keep mutation report layout """

    def __init__(self, *args, **kwargs):
        super(ReportLayout, self).__init__(*args, **kwargs)
        self.__init_properties()
        self.__init_cell_colors()

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(ReportLayout, self).get_raw_obj_str(*args, **kwargs)
        raw_repr["annotated tabix file"] = self.annotated_vcf_tabix
        raw_repr["annotation columns"] = self.anno_cols
        if self.anno_excl_tags is not None and len(self.anno_excl_tags) > 0:
            raw_repr["annotation exclusion tags"] = self.anno_excl_tags
        raw_repr["genotyping calling detail"] = self.call_detail
        raw_repr["genotyping calling quality"] = self.call_gq
        raw_repr[RPT_LAYOUT_CAPTION_FREQ_RATIOS] = self.freq_ratios
        if self.exprs is not None:
            raw_repr["expression actions"] = self.exprs.actions
        raw_repr["genes list"] = self.filter_genes
        if self.report_regions is None:
            raw_repr[RPT_LAYOUT_CAPTION_RPT_REGIONS] = "ALL"
        else:
            raw_repr[RPT_LAYOUT_CAPTION_RPT_REGIONS] = self.report_regions
        raw_repr["split chromosome"] = self.split_chrom
        filter_actions = OrderedDict()
        filter_actions[FILTER_RARE] = self.filter_rare
        filter_actions[FILTER_PASS_VQSR] = self.filter_pass_vqsr
        filter_actions[FILTER_NON_INTERGENIC] = self.filter_non_intergenic
        filter_actions[FILTER_NON_INTRONIC] = self.filter_non_intronic
        filter_actions[FILTER_NON_UPSTREAM] = self.filter_non_upstream
        filter_actions[FILTER_NON_DOWNSTREAM] = self.filter_non_downtream
        filter_actions[FILTER_NON_UTR] = self.filter_non_utr
        filter_actions[FILTER_NON_SYNONYMOUS] = self.filter_non_synonymous
        filter_actions[FILTER_HAS_MUTATION] = self.filter_has_mutation
        filter_actions[FILTER_HAS_SHARED] = self.filter_has_shared
        filter_actions[FILTER_NON_RECESSIVE_GENE] = self.filter_non_recessive_gene
        raw_repr[RPT_LAYOUT_CAPTION_FILTER_ACTIONS] = filter_actions
#        report_exclusion = {}
#        report_exclusion['only summary'] = self.only_summary
#        report_exclusion['only families'] = self.only_families
#        raw_repr['report exclusion'] = report_exclusion
        return raw_repr

    def __init_properties(self):
        self.__anno_cols = None
        self.__freq_ratios = None
        self.__header_corrections = self.__get_header_corrections()
        self.__exprs = self._get_job_config(JOBS_SETUP_RPT_EXPRESSIONS_KEY)
        if self.__exprs is not None:
            self.__exprs = VcfExpressions(self.__exprs)

    def __get_header_corrections(self):
        raw_header_corrections = self._get_job_config(JOBS_SETUP_RPT_HEADER_CORRECTIONS_KEY)
        header_corrections = {}
        if raw_header_corrections is not None:
            for header_correction in raw_header_corrections:
                old_header = header_correction[JOBS_SETUP_RPT_OLD_HEADER_KEY]
                new_header = header_correction[JOBS_SETUP_RPT_NEW_HEADER_KEY]
                header_corrections[old_header] = new_header
        return header_corrections

    def __cal_anno_cols(self):
        # open the annotated vcf tabix file for checking
        # if the expected output columns are exist
        vcf_reader = VcfReader(filename=self.annotated_vcf_tabix)
        vcf_record = vcf_reader.next()
        # generate columns list
        anno_cols = []
        for col_name in self._get_job_config(JOBS_SETUP_RPT_ANNO_COLS_KEY,
                                             required=True):
            # excluce columns based on configuration
            excluded = False
            for excl_tag in self.anno_excl_tags:
                if excl_tag in ALL_MUTREP_ANNO_COLS[col_name]:
                    excluded = True
                    break
            if excluded:
                continue
            # exclude columns that are not actually in the annotated tabix file
            # the underlying assumption is that all the variants must have the 
            # same number of fields annotated by ANNOVAR
            if (col_name not in vcf_record.INFO.keys() and
                col_name not in EXAC03_CONSTRAINT_COL_NAMES and
                col_name not in EST_KVOT_COLS and
                col_name != PATHOGENIC_COUNT_COL_NAME and
                col_name != INTERVAR_CLASS_COL_NAME and
                col_name != INTERVAR_EVIDENCE_COL_NAME
                ):
                self.warning("Columns " + col_name + " is missing")
                continue
            if ((col_name == EST_KVOT_EARLYONSET_VS_BRC_COL_NAME) and
                ((WES294_OAF_EARLYONSET_AF_COL_NAME not in vcf_record.INFO.keys()) or
                 (WES294_OAF_BRCS_AF_COL_NAME not in vcf_record.INFO.keys())
                 )
                ):
                self.warning(col_name + " cannot be calculated")
                continue
            if ((col_name == EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME) and
                ((WES294_OAF_EARLYONSET_AF_COL_NAME not in vcf_record.INFO.keys()) or
                 (EXAC_NFE_COL_NAME not in vcf_record.INFO.keys())
                 )
                ):
                self.warning(col_name + " cannot be calculated")
                continue
            if ((col_name == EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME) and
                ((WES294_OAF_EARLYONSET_AF_COL_NAME not in vcf_record.INFO.keys()) or
                 (KG2014OCT_EUR_COL_NAME not in vcf_record.INFO.keys())
                 )
                ):
                self.warning(col_name + " cannot be calculated")
                continue
            if (col_name in EXAC03_CONSTRAINT_COL_NAMES and
                EXAC03_CONSTRAINT_COL_NAME not in vcf_record.INFO.keys()
                ):
                self.warning(col_name + " cannot be calculated")
                continue
            if ((col_name == INTERVAR_CLASS_COL_NAME or
                col_name == INTERVAR_EVIDENCE_COL_NAME) and
                INTERVAR_AND_EVIDENCE_COL_NAME not in vcf_record.INFO.keys()
                ):
                self.warning(col_name + " cannot be calculated")
                continue
            # here are the columns that can be shown without errors
            anno_cols.append(col_name)
        return anno_cols

    @property
    def annotated_vcf_tabix(self):
        return self._get_job_config(JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX,
                                    required=True)

    @property
    def anno_excl_tags(self):
        return self._get_job_config(JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY,
                                    default_val=[])

    @property
    def anno_cols(self):
        if self.__anno_cols is None:
            self.__anno_cols = self.__cal_anno_cols()
        return self.__anno_cols

    @property
    def report_regions(self):
        regions = self._get_job_config(JOBS_SETUP_RPT_REGIONS_KEY)
        if regions is not None:
            regions = map(lambda x: DNARegion(str(x)),
                          regions)
        return regions

    @property
    def exprs(self):
        return self.__exprs

    @property
    def split_chrom(self):
        return self._get_job_config(JOBS_SETUP_RPT_SPLIT_CHROM_KEY,
                                    default_val=False)

    @property
    def call_detail(self):
        return JOBS_SETUP_RPT_EXTRA_ANNO_CALL_DETAIL in self._get_job_config(JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY,
                                                                      default_val=[])

    @property
    def call_gq(self):
        return JOBS_SETUP_RPT_EXTRA_ANNO_CALL_GQ in self._get_job_config(JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY,
                                                                  default_val=[])

    @property
    def filter_rare(self):
        return JOBS_SETUP_RPT_FILTER_RARE in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                  default_val=[])

    @property
    def filter_pass_vqsr(self):
        return JOBS_SETUP_RPT_FILTER_PASS_VQSR in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                       default_val=[])

    @property
    def filter_non_intergenic(self):
        return JOBS_SETUP_RPT_FILTER_NON_INTERGENIC in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                            default_val=[])

    @property
    def filter_non_intronic(self):
        return JOBS_SETUP_RPT_FILTER_NON_INTRONIC in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                          default_val=[])

    @property
    def filter_non_upstream(self):
        return JOBS_SETUP_RPT_FILTER_NON_UPSTREAM in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                          default_val=[])

    @property
    def filter_non_downtream(self):
        return JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                            default_val=[])

    @property
    def filter_non_utr(self):
        return JOBS_SETUP_RPT_FILTER_NON_UTR in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                     default_val=[])

    @property
    def filter_non_synonymous(self):
        return JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                            default_val=[])

    @property
    def filter_has_mutation(self):
        return JOBS_SETUP_RPT_FILTER_HAS_MUTATION in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                          default_val=[])

    @property
    def filter_has_shared(self):
        return JOBS_SETUP_RPT_FILTER_HAS_SHARED in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                        default_val=[])

    @property
    def filter_non_recessive_gene(self):
        return JOBS_SETUP_RPT_FILTER_NON_RECESSIVE_GENE in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
                                                                                default_val=[])

    @property
    def coloring_shared(self):
        return JOBS_SETUP_RPT_COLORING_SHARED in self._get_job_config(JOBS_SETUP_RPT_COLORING_SAMPLES_KEY,
                                                                      default_val=[])

    @property
    def coloring_zygosity(self):
        return JOBS_SETUP_RPT_COLORING_ZYGOSITY in self._get_job_config(JOBS_SETUP_RPT_COLORING_SAMPLES_KEY,
                                                                        default_val=[])

    @property
    def header_corrections(self):
        return self.__header_corrections

    @property
    def freq_ratios(self):
        freq_ratios = self._get_job_config(JOBS_SETUP_RPT_FREQ_RATIOS_KEY)
        if ((self.__freq_ratios is None) and
            (freq_ratios is not None)
            ):
            tmp_ratios = OrderedDict()
            for freq_ratio in freq_ratios:
                col = freq_ratio[JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY]
                freq = freq_ratio[JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY]
                tmp_ratios[col] = freq
            self.__freq_ratios = tmp_ratios
        elif self.__freq_ratios is None:
            self.__freq_ratios = []
        return self.__freq_ratios

    @property
    def filter_genes(self):
        return self._get_job_config(JOBS_SETUP_RPT_FILTER_GENES_KEY)

    def __init_cell_colors(self):
        self.__cell_colors = {}
        self.__cell_colors[CELL_TYPE_HET_SHARED] = DFLT_COLOR_HET_SHARED
        self.__cell_colors[CELL_TYPE_HOM_SHARED] = DFLT_COLOR_HOM_SHARED
        self.__cell_colors[CELL_TYPE_HET_ZYGO] = DFLT_COLOR_HET_ZYGO
        self.__cell_colors[CELL_TYPE_HOM_ZYGO] = DFLT_COLOR_HOM_ZYGO
        self.__cell_colors[CELL_TYPE_HET_RECESSIVE] = DFLT_COLOR_HET_RECESSIVE
        self.__cell_colors[CELL_TYPE_HOM_RECESSIVE] = DFLT_COLOR_HOM_RECESSIVE
        self.__cell_colors[CELL_TYPE_SEPARATOR] = DFLT_COLOR_SEPARATOR

    @property
    def cell_color_het_shared(self):
        return self.__cell_colors[CELL_TYPE_HET_SHARED]

    @property
    def cell_color_hom_shared(self):
        return self.__cell_colors[CELL_TYPE_HOM_SHARED]

    @property
    def cell_color_het_zygo(self):
        return self.__cell_colors[CELL_TYPE_HET_ZYGO]

    @property
    def cell_color_hom_zygo(self):
        return self.__cell_colors[CELL_TYPE_HOM_ZYGO]

    @property
    def cell_color_het_recessive(self):
        return self.__cell_colors[CELL_TYPE_HET_RECESSIVE]

    @property
    def cell_color_hom_recessive(self):
        return self.__cell_colors[CELL_TYPE_HOM_RECESSIVE]

    @property
    def cell_color_separator(self):
        return self.__cell_colors[CELL_TYPE_SEPARATOR]

    @property
    def show_shared_mutations(self):
        return self._get_job_config(JOBS_SETUP_RPT_SHOW_SHARED_MUTATIONS,
                                    default_val=False)

#    @property
#    def only_summary(self):
#        return self._get_job_config(JOBS_SETUP_RPT_ONLY_SUMMARY_KEY,
#                                    default_val=False)
#
#    @property
#    def only_families(self):
#        return self._get_job_config(JOBS_SETUP_RPT_ONLY_FAMILIES_KEY,
#                                    default_val=False)

class MutRepPipeline(CMMPipeline):
    """ A class to control CMMDB best practice pipeline """

    def __init__(self, *args, **kwargs):
        super(MutRepPipeline, self).__init__(*args, **kwargs)
        self.__init_properties()

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(MutRepPipeline, self).get_raw_obj_str(*args, **kwargs)
        return raw_repr

    def get_third_party_software_version(self):
        vm = VersionManager()
        versions = OrderedDict()
        versions['pyvcf'] = vm.pyvcf_version
        versions['pyaml'] = vm.pyaml_version
        versions['xlsxwriter'] = vm.xlsxwriter_version
        return versions

    def __init_properties(self):
        self.__report_layout = None

    @property
    def report_layout(self):
        if self.__report_layout is None:
            self.__report_layout = ReportLayout(entries=self._get_job_config(JOBS_SETUP_RPT_LAYOUT_SECTION,
                                                                             required=True))
        return self.__report_layout

    @property
    def summary_rpt_file(self):
        return join_path(self.rpts_out_dir,
                         self.dataset_name+"_summary.xlsx")

    @property
    def plain_fmts(self):
        return self.__plain_fmts

    def __set_report_format(self):
        self.__plain_fmts = self.__wb.add_colors_format({})

    def __add_sheet(self, sheet_name):
        ws = self.__wb.add_worksheet(sheet_name)
        ws.set_default_row(12)
        return ws

    def correct_header(self, old_header):
        if old_header in self.report_layout.header_corrections:
            return self.report_layout.header_corrections[old_header]
        return old_header

    def __write_header(self, ws, samples_id, samples_header):
        def write_samples_header(start_ncol_idx,
                                 samples_header,
                                 cell_fmt,
                                 coloring_samples_idx=[],
                                 ):
            last_col_idx = start_ncol_idx
            for sample_idx in xrange(len(samples_header)):
                sample_id = samples_header[sample_idx]
                if sample_idx in coloring_samples_idx:
                    unaffected_sample_color = self.report_layout.cell_color_separator
                    sample_cell_fmt = self.plain_fmts[unaffected_sample_color]
                else:
                    sample_cell_fmt = cell_fmt
                ws.write(0, last_col_idx, sample_id, sample_cell_fmt)
                last_col_idx += 1
                if self.report_layout.call_gq:
                    ws.write(0, last_col_idx, sample_id+"_GQ", cell_fmt)
                    last_col_idx += 1
                if self.report_layout.call_detail:
                    ws.write(0, last_col_idx, sample_id+"(detail)", cell_fmt)
                    last_col_idx += 1
            return last_col_idx
        anno_cols = self.report_layout.anno_cols
        cell_fmt = self.plain_fmts[NO_COLOR]
        ws.write(0, XLS_CHROM_COL_IDX, "CHROM", cell_fmt)
        ws.write(0, XLS_POS_COL_IDX, "POS", cell_fmt)
        ws.write(0, XLS_REF_COL_IDX, "REF", cell_fmt)
        ws.write(0, XLS_ALT_COL_IDX, "ALT", cell_fmt)
        ws.write(0, XLS_FILTER_COL_IDX, "FILTER", cell_fmt)
        len_anno_cols = len(anno_cols)
        for anno_idx in xrange(len_anno_cols):
            col_name = self.correct_header(anno_cols[anno_idx])
            ws.write(0, anno_idx+LAYOUT_VCF_COLS, col_name, cell_fmt)
        sample_start_idx = LAYOUT_VCF_COLS + len_anno_cols
        if (self.samples_groups is not None and
            (len(self.samples_groups) > 1 or 
             0 not in self.samples_groups)):
            last_col_idx = sample_start_idx
            for group_no in self.samples_groups.keys():
                if group_no == 0:
                    continue
                group_samples = self.samples_groups[group_no]
                samples = SamplesGroup(filter(lambda x: x.sample_id in samples_id,
                                              group_samples))
                last_col_idx += 1
                last_col_idx = write_samples_header(last_col_idx,
                                                    samples.ids_w_fam_pref,
                                                    cell_fmt,
                                                    )
        else:
            last_col_idx = write_samples_header(sample_start_idx,
                                                samples_header,
                                                cell_fmt,
                                                )
        return sample_start_idx, last_col_idx

    def __format_call_detail(self, call, call_fmt):
        call_detail = []
        for fmt_idx in xrange(len(call_fmt)):
            fmt = call_fmt[fmt_idx]
            data = call.data[fmt_idx]
            if type(data) is list:
                data = map(lambda x: str(x), data)
                call_detail.append(fmt + "=" + ",".join(data))
            else:
                if data is None:
                    data = "."
                call_detail.append(fmt + "=" + str(data))
        return ":".join(call_detail)

    def __write_zygosities(self,
                           ws,
                           row,
                           start_col,
                           vcf_record,
                           allele_idx,
                           samples_id,
                           dflt_cell_fmt,
                           recessive_samples_id=None,
                           ):
        last_col_idx = start_col
        for sample_idx in xrange(len(samples_id)):
            sample_id = samples_id[sample_idx]
            call = vcf_record.genotype(sample_id)
            zygo = call.cmm_gts[allele_idx]
            # determine cell(s) format
            if self.report_layout.coloring_shared:
                if type(call.shared_mutations) is not list:
                    zygo_fmt = dflt_cell_fmt
                elif not call.shared_mutations[allele_idx]:
                    zygo_fmt = dflt_cell_fmt
                # then assuming that they are shared
                elif call.actual_gts[allele_idx] == CMMGT_HOMOZYGOTE:
                    hom_shared_color = self.report_layout.cell_color_hom_shared
                    zygo_fmt = self.plain_fmts[hom_shared_color]
                else:
                    het_shared_color = self.report_layout.cell_color_het_shared
                    zygo_fmt = self.plain_fmts[het_shared_color]
            elif self.report_layout.coloring_zygosity:
                if call.actual_gts[allele_idx] == CMMGT_HOMOZYGOTE:
                    hom_zygo_color = self.report_layout.cell_color_hom_zygo
                    zygo_fmt = self.plain_fmts[hom_zygo_color]
                elif call.actual_gts[allele_idx] == CMMGT_HETEROZYGOTE:
                    het_zygo_color = self.report_layout.cell_color_het_zygo
                    zygo_fmt = self.plain_fmts[het_zygo_color]
                else:
                    zygo_fmt = dflt_cell_fmt
            else:
                zygo_fmt = dflt_cell_fmt

            # write content to cell(s)
            ws.write(row, last_col_idx, zygo, zygo_fmt)
            last_col_idx += 1
            if self.report_layout.call_gq:
                zygo_gq = call.data.GQ
                if zygo_gq is None:
                    zygo_gq = "."
                ws.write(row, last_col_idx, zygo_gq, zygo_fmt)
                last_col_idx += 1
            if self.report_layout.call_detail:
                formatted_call = self.__format_call_detail(call, vcf_record.FORMAT.split(":"))
                ws.write(row, last_col_idx, formatted_call, zygo_fmt)
                last_col_idx += 1
        return last_col_idx

    def __write_content(self,
                        ws,
                        row,
                        vcf_record,
                        allele_idx,
                        samples_id,
                        ):
        anno_cols = self.report_layout.anno_cols
        # use input expression to determine if row will have background color
        row_color = None
        if self.report_layout.exprs is not None:
            color_row_actions = self.report_layout.exprs.actions[ACTION_COLOR_ROW]
            for cra in color_row_actions:
                if vcf_record.vcf_eval(cra.pattern, allele_idx):
                    row_color = cra.color
                    break
        if row_color is not None:
            dflt_cell_fmt = self.plain_fmts[row_color]
        else:
            dflt_cell_fmt = self.plain_fmts[NO_COLOR]
        # start writing content
        alt_allele = vcf_record.alleles[allele_idx]
        if is_number(vcf_record.CHROM):
            ws.write(row, XLS_CHROM_COL_IDX, int(vcf_record.CHROM), dflt_cell_fmt)
        else:
            ws.write(row, XLS_CHROM_COL_IDX, vcf_record.CHROM, dflt_cell_fmt)
        ws.write(row, XLS_POS_COL_IDX, vcf_record.POS, dflt_cell_fmt)
        ws.write(row, XLS_REF_COL_IDX, vcf_record.REF, dflt_cell_fmt)
        ws.write(row, XLS_ALT_COL_IDX, str(alt_allele), dflt_cell_fmt)
        ws.write(row, XLS_FILTER_COL_IDX, vcf_record.FILTER, dflt_cell_fmt)
        len_anno_cols = len(anno_cols)
        # annotate INFO columns
        if self.report_layout.exprs is not None:
            color_col_actions = self.report_layout.exprs.actions[ACTION_COLOR_COL]
        else:
            color_col_actions = {}
        # for each INFO column
        for anno_idx in xrange(len_anno_cols):
            anno_col_name = anno_cols[anno_idx]
            info = vcf_record.get_info(anno_col_name, allele_idx)
            if (info == "" or
                info is None
                ):
                info = ""
            # determine if the column needed to be colored even though it's empty
            info_cell_fmt = dflt_cell_fmt
            color_col = False
            if anno_col_name in color_col_actions:
                anno_color_col_actions = color_col_actions[anno_col_name]
                for cca in anno_color_col_actions:
                    if vcf_record.vcf_eval(cca.pattern, allele_idx):
                        info_cell_fmt = self.plain_fmts[cca.color]
                        color_col = True
            if info == "" and not color_col:
                continue
            if anno_col_name in PREDICTION_COLS:
                info = info.description
            # determine cell format
            if (anno_col_name in FORMAT_COLS and
                FORMAT_COLS[anno_col_name] == FORMAT_COL_FLOAT and
                info != "INF" and
                is_number(info)
                ):
                ws.write(row, anno_idx+LAYOUT_VCF_COLS, float(info), info_cell_fmt)
            elif (anno_col_name in FORMAT_COLS and
                FORMAT_COLS[anno_col_name] == FORMAT_COL_INT and
                is_number(info)
                ):
                ws.write(row, anno_idx+LAYOUT_VCF_COLS, int(info), info_cell_fmt)
            elif is_number(info):
                ws.write(row, anno_idx+LAYOUT_VCF_COLS, info, info_cell_fmt)
            else:
                ws.write(row, anno_idx+LAYOUT_VCF_COLS, str(info).decode('utf-8'), info_cell_fmt)
        last_col_idx = LAYOUT_VCF_COLS + len_anno_cols
        # annotate samples information
        if (self.samples_groups is not None and
            (len(self.samples_groups) > 1 or 
             0 not in self.samples_groups)):
            sep_color = self.report_layout.cell_color_separator
            sep_fmt = self.plain_fmts[sep_color]
            for group_no in self.samples_groups.keys():
                if group_no == 0:
                    continue
                group_samples = self.samples_groups[group_no]
                samples = SamplesGroup(filter(lambda x: x.sample_id in samples_id,
                                              group_samples))
                ws.write(row, last_col_idx, "", sep_fmt)
                last_col_idx += 1
                last_col_idx = self.__write_zygosities(ws,
                                                       row,
                                                       last_col_idx,
                                                       vcf_record,
                                                       allele_idx,
                                                       samples.ids,
                                                       dflt_cell_fmt,
                                                       )
        else:
            last_col_idx = self.__write_zygosities(ws,
                                                   row,
                                                   last_col_idx,
                                                   vcf_record,
                                                   allele_idx,
                                                   samples_id,
                                                   dflt_cell_fmt,
                                                   )
        # logging
        if row % RECORDS_LOG_INTERVAL == 0:
            log_msg = str(row)
            log_msg += " records were written to the sheet"
            self.info(log_msg)
        return row + 1

    def __write_contents(self,
                         ws,
                         row,
                         vcf_records,
                         samples_id,
                         check_shared,
                         ):
        for vcf_record in vcf_records:
            for allele_idx in xrange(1, len(vcf_record.alleles)):
                if (check_shared and
                    not vcf_record.is_shared(samples_id, allele_idx)):
                    continue
                if (self.report_layout.filter_rare and
                    not vcf_record.is_rare(allele_idx=allele_idx)):
                    continue
                if (self.report_layout.filter_pass_vqsr and
                    not vcf_record.is_pass_vqsr(allele_idx=allele_idx)):
                    continue
                if (self.report_layout.filter_non_intergenic and
                    vcf_record.is_intergenic[allele_idx]):
                    continue
                if (self.report_layout.filter_non_intronic and
                    vcf_record.is_intronic[allele_idx]):
                    continue
                if (self.report_layout.filter_non_upstream and
                    vcf_record.is_upstream[allele_idx]):
                    continue
                if (self.report_layout.filter_non_downtream and
                    vcf_record.is_downstream[allele_idx]):
                    continue
                if (self.report_layout.filter_non_utr and
                    vcf_record.is_utr[allele_idx]):
                    continue
                if (self.report_layout.filter_non_synonymous and
                    vcf_record.is_synonymous[allele_idx]):
                    continue
                if (self.report_layout.filter_has_mutation and
                    not vcf_record.has_mutation(samples_id, allele_idx)):
                    continue
                if (self.report_layout.filter_has_shared and
                    not vcf_record.has_shared(allele_idx)):
                    continue
                if vcf_record.alleles[allele_idx] == "*":
                    continue
                delete_row = False
                if self.report_layout.exprs is not None:
                    del_row_actions = self.report_layout.exprs.actions[ACTION_DELETE_ROW]
                    for dra in del_row_actions:
                        if vcf_record.vcf_eval(dra.pattern, allele_idx):
                            delete_row = True
                            break
                if delete_row:
                    continue
                # If filter_genes is defined, only the filter_genes can be
                # in the report. Otherwise, no filtering applied
                if self.report_layout.filter_genes is not None: 
                    gene_found = False
                    gene_refgenes = vcf_record.get_info(GENE_REFGENE_COL_NAME,
                                                        allele_idx).split(',')
                    for filter_gene in self.report_layout.filter_genes:
                        if filter_gene in gene_refgenes:
                            gene_found = True
                            break
                    if not gene_found:
                        continue
                row = self.__write_content(ws,
                                           row,
                                           vcf_record,
                                           allele_idx,
                                           samples_id)
        return row

    def __set_layout(self, ws, sample_start_idx, record_size):
#        ws.autofilter(0, 0, 0, sample_start_idx-1)
        ws.autofilter(0, 0, 0, record_size-1)
        ws.set_column(XLS_CHROM_COL_IDX, XLS_CHROM_COL_IDX, 2.5)
        ws.set_column(XLS_POS_COL_IDX, XLS_POS_COL_IDX, 9)
        ws.set_column(XLS_REF_COL_IDX, XLS_REF_COL_IDX, 3.5)
        ws.set_column(XLS_ALT_COL_IDX, XLS_ALT_COL_IDX, 3.5)
        ws.set_column(XLS_FILTER_COL_IDX, XLS_FILTER_COL_IDX, 3.5)
        ws.set_column(sample_start_idx, record_size-1, 3)
        row_freeze_idx = 1
        anno_cols = self.report_layout.anno_cols
        if GENE_REFGENE_COL_NAME in anno_cols:
            col_freeze_idx = anno_cols.index(GENE_REFGENE_COL_NAME) + LAYOUT_VCF_COLS + 1
        else:
            col_freeze_idx = 0
        ws.freeze_panes(row_freeze_idx, col_freeze_idx)

    def __add_muts_sheet(self,
                         sheet_name,
                         report_regions,
                         samples_id=None,
                         samples_header=None,
                         check_shared=False,
                         ):
        wb = self.__wb
        annotated_vcf_tabix = self.report_layout.annotated_vcf_tabix
        if annotated_vcf_tabix.endswith('.vcf.gz'):
            vcf_reader = VcfReader(filename=annotated_vcf_tabix,
                                   family_infos=self.families_info,
                                   freq_ratios=self.report_layout.freq_ratios,
                                   )
        else:
            self.thrown(annotated_vcf_tabix + ' does not endswith .vcf.gz')
        row = 1
        if samples_id is None:
            samples_id = vcf_reader.samples
        if samples_header is None:
            samples_header = samples_id
        ws = self.__add_sheet(sheet_name)
        sample_start_idx, ncol = self.__write_header(ws,
                                                     samples_id,
                                                     samples_header,
                                                     )
        if report_regions is None:
            row = self.__write_contents(ws,
                                        row,
                                        vcf_reader,
                                        samples_id,
                                        check_shared)
        else:
            for report_region in report_regions:
                if report_region.start_pos is None:
                    vcf_records = vcf_reader.fetch(report_region.chrom,
                                                   0,
                                                   1000000000)
                else:
                    vcf_records = vcf_reader.fetch(report_region.chrom,
                                                   int(report_region.start_pos),
                                                   int(report_region.end_pos),
                                                   )
                row = self.__write_contents(ws,
                                            row,
                                            vcf_records,
                                            samples_id,
                                            check_shared)
        log_msg = "Finish .. "
        log_msg += " total of " + str(row-1)
        log_msg += " records were written to the sheet"
        self.info(log_msg)
        self.__set_layout(ws, sample_start_idx, ncol)

    def __add_criteria_sheet(self,
                             ):
        def write_filter(row_idx,
                         filter_text,
                         ):
            ws.write(row_idx, 0, filter_text)
            return row_idx + 1
        
        row_idx = 0
        ws = self.__add_sheet('criteria')
        rpt_layout = self.report_layout.get_raw_obj_str()
        row_idx = write_filter(row_idx, "Criteria applied")
        caption_indent = " " * 5
        criteria_indent = " " * 10
        for rpt_caption in rpt_layout:
            if rpt_caption == RPT_LAYOUT_CAPTION_RPT_REGIONS:
                row_idx = write_filter(row_idx, caption_indent+"DNA region(s)")
                rpt_regions = rpt_layout[rpt_caption]
                if type(rpt_regions) is list:
                    for rpt_region in rpt_regions:
                        row_idx = write_filter(row_idx, criteria_indent+"- "+str(rpt_region))
                else:
                    row_idx = write_filter(row_idx, criteria_indent+"- "+rpt_regions)
            if rpt_caption == RPT_LAYOUT_CAPTION_FILTER_ACTIONS:
                filter_actions = rpt_layout[rpt_caption]
                # check if there is any filter action applied
                has_action = False
                for filter_action in filter_actions:
                    if filter_actions[filter_action] is True:
                        has_action = True
                        break
                # if there is any action, then log the info
                if has_action:
                    row_idx = write_filter(row_idx, caption_indent+"filter")
                    for filter_action in filter_actions:
                        if filter_actions[filter_action] is False:
                            continue
                        if filter_action == FILTER_RARE:
                            row_idx = write_filter(row_idx, criteria_indent+"- Only rare variants (MAF 1000G EUR < 20%")
                        if filter_action == FILTER_NON_INTERGENIC:
                            row_idx = write_filter(row_idx, criteria_indent+"- No intergenic variants")
                        if filter_action == FILTER_NON_INTRONIC:
                            row_idx = write_filter(row_idx, criteria_indent+"- No intronic variants")
                        if filter_action == FILTER_NON_UPSTREAM:
                            row_idx = write_filter(row_idx, criteria_indent+"- No upstream variants")
                        if filter_action == FILTER_NON_DOWNSTREAM:
                            row_idx = write_filter(row_idx, criteria_indent+"- No downstream variants")
                        if filter_action == FILTER_NON_UTR:
                            row_idx = write_filter(row_idx, criteria_indent+"- No UTR3 and UTR5 variants")
                        if filter_action == FILTER_NON_SYNONYMOUS:
                            row_idx = write_filter(row_idx, criteria_indent+"- No synonymous variants")
                        if filter_action == FILTER_HAS_MUTATION:
                            row_idx = write_filter(row_idx, criteria_indent+"- Variants must be present in at least one of the samples in the sheet")
                        if filter_action == FILTER_HAS_SHARED:
                            row_idx = write_filter(row_idx, criteria_indent+"- Variants must be shared by all of the samples in the sheet")
                    
    def gen_summary_report(self,
                           report_regions,
                           out_file=None,
                           ):
        if out_file is None:
            out_file = self.summary_rpt_file
        self.info("")
        self.info(" >>>> generating summary report")
        self.info(" >>>> report file: " + out_file)
        self.__wb = Workbook(filename=out_file)
        self.__set_report_format()
        self.info("")
        self.info(" >> add 'summary_all' sheet")
        self.__add_muts_sheet("summary_all",
                              report_regions,
                              samples_id=self.samples_id,
                              samples_header=self.samples_id_w_fam_pref,
                              )
        self.__add_criteria_sheet()
        self.__wb.close()

    def __submit_report_jobs(self,
                             job_script,
                             job_params_prefix,
                             job_name_prefix,
                             slurm_log_prefix,
                             out_file_prefix,
                             report_regions,
                             ):
        if self.report_layout.split_chrom:
            if report_regions is None:
                report_regions = map(lambda x: DNARegion(x), ALL_CHROMS)
            for report_region in report_regions:
                chr_slurm_log_prefix = slurm_log_prefix
                chr_slurm_log_prefix += "_" + report_region.region_key
                job_name = job_name_prefix
                job_name += "_" + report_region.region_key
                out_file = out_file_prefix
                out_file += "_" + report_region.region_key
                out_file += ".xlsx"
                job_params = job_params_prefix
                job_params += " -r " + report_region.raw_region
                job_params += " -o " + out_file
                self.submit_job(job_name,
                                self.project_code,
                                "core",
                                "4",
                                self.job_alloc_time,
                                chr_slurm_log_prefix,
                                job_script,
                                job_params,
                                )
        else:
            job_name = job_name_prefix
            out_file = out_file_prefix + ".xlsx"
            job_params = job_params_prefix
            if report_regions is not None:
                job_params += " -r " + ",".join(map(lambda x: x.raw_region,
                                                    report_regions))
            job_params += " -o " + out_file
            self.submit_job(job_name,
                            self.project_code,
                            "core",
                            "4",
                            self.job_alloc_time,
                            slurm_log_prefix,
                            job_script,
                            job_params,
                            )

    def gen_summary_reports(self):
        report_regions = self.report_layout.report_regions
        if self.project_code is None:
            self.gen_summary_report(report_regions)
        else:
            slurm_log_prefix = join_path(self.slurm_log_dir,
                                         self.dataset_name)
            slurm_log_prefix += '_rpts_smy'
            job_name_prefix = self.dataset_name + '_rpts_smy'
            job_script = MUTREP_SUMMARY_REPORT_BIN
            job_params_prefix = " -j " + self.jobs_setup_file
            out_file_prefix = join_path(self.rpts_out_dir,
                                        self.dataset_name+"_summary")
            self.__submit_report_jobs(job_script,
                                      job_params_prefix,
                                      job_name_prefix,
                                      slurm_log_prefix,
                                      out_file_prefix,
                                      report_regions,
                                      )

    def gen_family_report(self,
                          fam_id,
                          report_regions,
                          out_file=None,
                          ):
        fam_info = self.families_info[fam_id]
        if out_file is None:
            out_file = join_path(self.rpts_out_dir,
                                 self.dataset_name+"_fam"+fam_info.fam_id+".xlsx")
        self.info("")
        self.info(" >>>> generating family report for family " + fam_info.fam_id)
        self.info(" >>>> report file: " + out_file)
        samples_id = map(lambda x: x.sample_id, fam_info.members)
        self.__wb = Workbook(filename=out_file)
        self.__set_report_format()
        if len(samples_id) == 1:
            self.info("")
            self.info(" >> add '" + samples_id[0] + "' sheet")
            self.__add_muts_sheet(samples_id[0],
                                  report_regions,
                                  samples_id=samples_id,
                                  check_shared=True,
                                  )
        else:
            self.info("")
            self.info(" >> add 'shared' sheet")
            self.__add_muts_sheet("shared",
                                  report_regions,
                                  samples_id=samples_id,
                                  check_shared=True)
            for sample_name in samples_id:
                self.info("")
                self.info(" >> add '" + sample_name + "' sheet")
                self.__add_muts_sheet(sample_name,
                                      report_regions,
                                      samples_id=[sample_name],
                                      check_shared=True)
        self.__wb.close()

    def __gen_family_reports(self, fam_id):
        report_regions = self.report_layout.report_regions
        if self.project_code is None:
            self.gen_family_report(fam_id, report_regions)
        else:
            slurm_log_prefix = join_path(self.slurm_log_dir,
                                         self.dataset_name)
            slurm_log_prefix += '_rpts_fam'
            slurm_log_prefix += fam_id
            job_name_prefix = self.dataset_name + '_rpts_fam' + fam_id
            job_script = MUTREP_FAMILY_REPORT_BIN
            job_params_prefix = " -j " + self.jobs_setup_file
            job_params_prefix += " -f " + fam_id
            out_file_prefix = join_path(self.rpts_out_dir,
                                        self.dataset_name+"_fam"+fam_id)
            self.__submit_report_jobs(job_script,
                                      job_params_prefix,
                                      job_name_prefix,
                                      slurm_log_prefix,
                                      out_file_prefix,
                                      report_regions,
                                      )

    def gen_families_reports(self):
        if self.families_info is None:
            return
        for fam_id in self.families_info:
            self.__gen_family_reports(fam_id)

    def __gen_reports(self):
        self.gen_summary_reports()
#        if self.report_layout.only_summary:
#            self.gen_summary_reports()
#        elif self.report_layout.only_families:
#            self.gen_families_reports()
#        else:
#            self.gen_families_reports()
#            self.gen_summary_reports()

    def monitor_init(self, *args, **kwargs):
        self.__gen_reports()
        super(MutRepPipeline, self).monitor_init(*args, **kwargs)

    def run_offline_pipeline(self):
        self.__gen_reports()

def create_jobs_setup_file(*args, **kwargs):
    job_setup_document, stream = init_jobs_setup_file(*args, **kwargs)
    rpt_cfg = {}
    anno_cols = get_func_arg('anno_cols', kwargs)
    rpt_cfg[JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX] = get_func_arg('annotated_vcf_tabix',
                                                               kwargs)
    if anno_cols is None:
        rpt_cfg[JOBS_SETUP_RPT_ANNO_COLS_KEY] = ALL_MUTREP_ANNO_COLS.keys()
    else:
        rpt_cfg[JOBS_SETUP_RPT_ANNO_COLS_KEY] = anno_cols.split(",")
    anno_excl_tags = get_func_arg('anno_excl_tags', kwargs)
    if anno_excl_tags is not None:
        rpt_cfg[JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY] = anno_excl_tags.split(",")
    header_corrections = get_func_arg('header_corrections', kwargs, default_val=DFLT_HEADER_CORRECTIONS)
    if header_corrections is not None:
        hcs = []
        for raw_header_correction in header_corrections.split(","):
            header_correction = {} 
            old_header, new_header = raw_header_correction.split(":")
            header_correction[JOBS_SETUP_RPT_OLD_HEADER_KEY] = old_header
            header_correction[JOBS_SETUP_RPT_NEW_HEADER_KEY] = new_header
            hcs.append(header_correction)
        rpt_cfg[JOBS_SETUP_RPT_HEADER_CORRECTIONS_KEY] = hcs
    report_regions = get_func_arg('report_regions', kwargs)
    if report_regions is not None:
        rpt_cfg[JOBS_SETUP_RPT_REGIONS_KEY] = report_regions.split(",")
    expression_patterns = get_func_arg('expression_patterns', kwargs)
    if expression_patterns is not None:
        exprs = []
        for raw_pattern in  expression_patterns.split(","):
            name, pattern = raw_pattern.split(":")
            expr = DefaultOrderedDict(list)
            expr[JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY] = name.strip()
            expr[JOBS_SETUP_RPT_EXPRESSIONS_PATTERN_KEY] = "'" + pattern.strip().replace("'", "''") + "'"
            exprs.append(expr)
        expression_usages = get_func_arg('expression_usages', kwargs)
        if expression_usages is not None:
            for raw_usage in expression_usages.split(","):
                usage = {}
                usage_expr_name = raw_usage.split(":")[0]
                usage[JOBS_SETUP_RPT_EXPRESSIONS_ACTION_KEY] = raw_usage.split(":")[1]
                if len(raw_usage.split(":")) > 2:
                    usage[JOBS_SETUP_RPT_EXPRESSIONS_INFO_KEY] = ":".join(raw_usage.split(":")[2:])
                for expr in exprs:
                    if expr[JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY] == usage_expr_name:
                        expr[JOBS_SETUP_RPT_EXPRESSIONS_USAGES_KEY].append(usage)
        rpt_cfg[JOBS_SETUP_RPT_EXPRESSIONS_KEY] = exprs
    rpt_cfg[JOBS_SETUP_RPT_SPLIT_CHROM_KEY] = get_func_arg('split_chrom',
                                                           kwargs,
                                                           default_val=False,
                                                           )
#    rpt_cfg[JOBS_SETUP_RPT_ONLY_SUMMARY_KEY] = get_func_arg('only_summary',
#                                                            kwargs,
#                                                            default_val=False,
#                                                            )
#    rpt_cfg[JOBS_SETUP_RPT_ONLY_FAMILIES_KEY] = get_func_arg('only_families',
#                                                             kwargs,
#                                                             default_val=False,
#                                                             )
    call_detail = get_func_arg('call_detail', kwargs, default_val=False)
    call_gq = get_func_arg('call_gq', kwargs, default_val=False)
    if call_detail or call_gq:
        extra_anno_cols = []
        if call_detail:
            extra_anno_cols.append(JOBS_SETUP_RPT_EXTRA_ANNO_CALL_DETAIL)
        if call_gq:
            extra_anno_cols.append(JOBS_SETUP_RPT_EXTRA_ANNO_CALL_GQ)
        rpt_cfg[JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY] = extra_anno_cols
    rpt_cfg[JOBS_SETUP_RPT_SHOW_SHARED_MUTATIONS] = get_func_arg('show_shared_mutations',
                                                                 kwargs,
                                                                 default_val=False,
                                                                 )
    coloring_shared = get_func_arg('coloring_shared', kwargs, default_val=False)
    coloring_zygosity = get_func_arg('coloring_zygosity', kwargs, default_val=False)
    if coloring_shared or coloring_zygosity:
        coloring_samples = []
        if coloring_shared:
            coloring_samples.append(JOBS_SETUP_RPT_COLORING_SHARED)
        if coloring_zygosity:
            coloring_samples.append(JOBS_SETUP_RPT_COLORING_ZYGOSITY)
        rpt_cfg[JOBS_SETUP_RPT_COLORING_SAMPLES_KEY] = coloring_samples
    rows_filter_actions = get_func_arg('rows_filter_actions', kwargs)
    if rows_filter_actions is not None:
        filter_criterias = []
        for filter_criteria in rows_filter_actions.split(","):
            if filter_criteria == JOBS_SETUP_RPT_FILTER_RARE:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_PASS_VQSR:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_INTERGENIC:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_INTRONIC:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_UPSTREAM:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_UTR:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_HAS_MUTATION:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_HAS_SHARED:
                filter_criterias.append(filter_criteria)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_RECESSIVE_GENE:
                filter_criterias.append(filter_criteria)
        rpt_cfg[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY] = filter_criterias
    frequency_ratios = get_func_arg('frequency_ratios',
                                    kwargs,
                                    default_val=DFLT_MUTREP_FREQ_RATIOS,
                                    )
    job_freq_ratios = []
    for frequency_ratio in frequency_ratios.split(","):
        (col, freq) = frequency_ratio.split(":")
        job_freq_ratios.append({JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY: col,
                                JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY: freq,
                                })
    rpt_cfg[JOBS_SETUP_RPT_FREQ_RATIOS_KEY] = job_freq_ratios
    filter_genes = get_func_arg('filter_genes', kwargs)
    if filter_genes is not None:
        rpt_cfg[JOBS_SETUP_RPT_FILTER_GENES_KEY] = filter_genes.split(",")
    job_setup_document[JOBS_SETUP_RPT_LAYOUT_SECTION] = rpt_cfg
    pyaml.dump(job_setup_document, stream)
