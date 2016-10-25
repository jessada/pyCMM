# This Python file uses the following encoding: utf-8
import re
import pyaml
from os.path import join as join_path
from collections import defaultdict
from collections import OrderedDict
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.settings import DFLT_MUTREP_FREQ_RATIOS
from pycmm.settings import PREDICTION_COLS
from pycmm.template import pyCMMBase
from pycmm.utils.ver import VersionManager
from pycmm.cmmlib import CMMParams
from pycmm.cmmlib.taparser import TAVcfReader as VcfReader
from pycmm.cmmlib.xlslib import CMMWorkbook as Workbook
from pycmm.cmmlib.xlslib import NO_COLOR
from pycmm.flow import get_func_arg
from pycmm.flow.cmmdb import CMMPipeline
from pycmm.flow import init_jobs_setup_file
from pycmm.settings import MUTREP_FAMILY_REPORT_BIN
from pycmm.settings import MUTREP_SUMMARY_REPORT_BIN
from pycmm.cmmlib.dnalib import ALL_CHROMS
from pycmm.cmmlib.tamodel import CMMGT_HOMOZYGOTE
from pycmm.cmmlib.tamodel import CMMGT_HETEROZYGOTE

ACTION_DELETE_ROW = "del_row"
ACTION_COLOR_ROW = "color_row"
ACTION_COLOR_COL = "color_col"

DFLT_COLOR_HET_SHARED = 'SILVER'
DFLT_COLOR_HOM_SHARED = 'GRAY25'
CELL_TYPE_HET_SHARED = 'HET_SHARED'
CELL_TYPE_HOM_SHARED = 'HOM_SHARED'

CHROM_POS_PATTERN = re.compile(r'''(?P<chrom>.+?):(?P<start_pos>.+?)-(?P<end_pos>.+)''')
LAYOUT_VCF_COLS = 4
RECORDS_LOG_INTERVAL = 1000

# *************** report layout section ***************
JOBS_SETUP_RPT_LAYOUT_SECTION = "REPORT_LAYOUT"
JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX = "ANNOTATED_VCF_TABIX"
JOBS_SETUP_RPT_ANNO_COLS_KEY = "COLUMNS"
JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY = "ANNOTATION_EXCLUSION_TAGS"
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
JOBS_SETUP_RPT_CALL_DETAIL_KEY = "Calling_detail"
JOBS_SETUP_RPT_MT_KEY = "Mitochondria"
JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY = "ROWS_FILTER_ACTIONS_CRITERIA"
JOBS_SETUP_RPT_FILTER_RARE = "Rare"
JOBS_SETUP_RPT_FILTER_NON_INTERGENIC = "Non-Intergenic"
JOBS_SETUP_RPT_FILTER_NON_INTRONIC = "Non-Intronic"
JOBS_SETUP_RPT_FILTER_NON_UPSTREAM = "Non-Upstream"
JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM = "Non-Downstream"
JOBS_SETUP_RPT_FILTER_NON_UTR = "Non-UTR"
JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS = "Non-Synonymous"
JOBS_SETUP_RPT_FILTER_HAS_MUTATION = "Has-Mutation"
JOBS_SETUP_RPT_FILTER_HAS_SHARED = "Has-Shared"
JOBS_SETUP_RPT_ONLY_SUMMARY_KEY = "ONLY_SUMMARY"
JOBS_SETUP_RPT_ONLY_FAMILIES_KEY = "ONLY_FAMILIES"

class ReportRegion(pyCMMBase):
    """ A structure to parse and keep mutation report region """

    def __init__(self,
                 raw_region,
                 **kwargs
                 ):
        super(ReportRegion, self).__init__(**kwargs)
        self.__parse_region(raw_region)
        self.__raw_region = raw_region

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["chrom"] = self.chrom
        if self.start_pos is not None:
            raw_repr["start position"] = self.start_pos
            raw_repr["end position"] = self.end_pos
        return raw_repr

    @property
    def chrom(self):
        return self.__chrom

    @property
    def start_pos(self):
        return self.__start_pos

    @property
    def end_pos(self):
        return self.__end_pos

    @property
    def raw_region(self):
        return self.__raw_region

    def __parse_region(self, raw_region):
        match = CHROM_POS_PATTERN.match(raw_region)
        if match is not None:
            self.__chrom = match.group('chrom')
            self.__start_pos = match.group('start_pos')
            self.__end_pos = match.group('end_pos')
        else:
            self.__chrom = raw_region
            self.__start_pos = None
            self.__end_pos = None

class ActionDelRow(pyCMMBase):
    """ A structure to parse action to delete a row from a mutation report sheet """

    def __init__(self,
                 pattern,
                 **kwargs
                 ):
        super(ActionDelRow, self).__init__(**kwargs)
        self.__pattern = pattern

    def get_raw_repr(self):
        return {"pattern": self.pattern}

    @property
    def pattern(self):
        return self.__pattern

class ActionColorRow(pyCMMBase):
    """ A structure to parse action to color a row in a mutation report sheet """

    def __init__(self,
                 pattern,
                 info,
                 **kwargs
                 ):
        super(ActionColorRow, self).__init__(**kwargs)
        self.__pattern = pattern
        self.__info = info

    def get_raw_repr(self):
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
                 **kwargs
                 ):
        super(ActionColorCol, self).__init__(**kwargs)
        self.__pattern = pattern
        self.__info = info

    def get_raw_repr(self):
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
                 **kwargs
                 ):
        super(VcfExpressions, self).__init__(**kwargs)
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
                            self.__actions[ACTION_COLOR_COL].append(ActionColorCol(pattern, info=info))
        return self.__actions

class ReportLayout(CMMParams):
    """ A structure to parse and keep mutation report layout """

    def __init__(self, **kwargs):
        super(ReportLayout, self).__init__(**kwargs)
        self.__init_properties()
        self.__init_cell_colors()

    def get_raw_repr(self, **kwargs):
        raw_repr = super(ReportLayout, self).get_raw_repr(**kwargs)
        raw_repr["annotated tabix file"] = self.annotated_vcf_tabix
        raw_repr["annotated columns"] = self.anno_cols
        if self.anno_excl_tags is not None and len(self.anno_excl_tags) > 0:
            raw_repr["annotation exclusion tags"] = self.anno_excl_tags
        raw_repr["genotyping calling detail"] = self.call_detail
        raw_repr["frequency ratio(s)"] = self.freq_ratios
        if self.report_regions is None:
            raw_repr["report region(s)"] = "ALL"
        else:
            raw_repr["report region(s)"] = self.report_regions
        raw_repr["split chromosome"] = self.split_chrom
        filter_actions = OrderedDict()
        filter_actions['Rare'] = self.filter_rare
        filter_actions['Non-Intergenic'] = self.filter_non_intergenic
        filter_actions['Non-Intronic'] = self.filter_non_intronic
        filter_actions['Non-Upstream'] = self.filter_non_upstream
        filter_actions['Non-Downstream'] = self.filter_non_downtream
        filter_actions['Non-UTR'] = self.filter_non_utr
        filter_actions['Non-Synonymous'] = self.filter_non_synonymous
        filter_actions['Has-mutation'] = self.filter_has_mutation
        filter_actions['Has-shared'] = self.filter_has_shared
        raw_repr['filter actions'] = filter_actions
        report_exclusion = {}
        report_exclusion['only summary'] = self.only_summary
        report_exclusion['only families'] = self.only_families
        raw_repr['report exclusion'] = report_exclusion
        raw_repr["summary_families sheet"] = self.summary_families_sheet
        return raw_repr

    def __init_properties(self):
        self.__anno_cols = None
        self.__freq_ratios = None
        self.__exprs = self._get_job_config(JOBS_SETUP_RPT_EXPRESSIONS_KEY)
        if self.__exprs is not None:
            self.__exprs = VcfExpressions(self.__exprs)

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
            if col_name not in vcf_record.INFO.keys():
                self.warning("Columns " + col_name + " is missing")
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
            regions = map(lambda x: ReportRegion(str(x)),
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
    def summary_families_sheet(self):
        return self._get_job_config(JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY,
                                    default_val=False)

    @property
    def call_detail(self):
        return JOBS_SETUP_RPT_CALL_DETAIL_KEY in self._get_job_config(JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY,
                                                                      default_val=[])

    @property
    def filter_rare(self):
        return JOBS_SETUP_RPT_FILTER_RARE in self._get_job_config(JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY,
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


        return self._get_job_config(JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY,
                                    default_val=False)
    @property
    def only_summary(self):
        return self._get_job_config(JOBS_SETUP_RPT_ONLY_SUMMARY_KEY,
                                    default_val=False)

    @property
    def only_families(self):
        return self._get_job_config(JOBS_SETUP_RPT_ONLY_FAMILIES_KEY,
                                    default_val=False)

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

    def __init_cell_colors(self):
        self.__cell_colors = {}
        self.__cell_colors[CELL_TYPE_HET_SHARED] = DFLT_COLOR_HET_SHARED
        self.__cell_colors[CELL_TYPE_HOM_SHARED] = DFLT_COLOR_HOM_SHARED

    @property
    def cell_color_het_shared(self):
        return self.__cell_colors[CELL_TYPE_HET_SHARED]

    @property
    def cell_color_hom_shared(self):
        return self.__cell_colors[CELL_TYPE_HOM_SHARED]

class MutRepPipeline(CMMPipeline):
    """ A class to control CMMDB best practice pipeline """

    def __init__(self, **kwargs):
        super(MutRepPipeline, self).__init__(**kwargs)
        self.__init_properties()

    def get_raw_repr(self, **kwargs):
        raw_repr = super(MutRepPipeline, self).get_raw_repr(**kwargs)
        return raw_repr

    def get_third_party_software_version(self):
        vm = VersionManager()
        versions = OrderedDict()
        versions['pyvcf'] = vm.pyvcf_version
        versions['pyaml'] = vm.pyaml_version
        versions['xlsxwriter'] = vm.xlsxwriter_version
        return versions

    def __init_properties(self):
        pass

    @property
    def report_layout(self):
        return ReportLayout(entries=self._get_job_config(JOBS_SETUP_RPT_LAYOUT_SECTION,
                                                         required=True))

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

    def __write_header(self, ws, samples_id):
        anno_cols = self.report_layout.anno_cols
        cell_fmt = self.plain_fmts[NO_COLOR]
        ws.write(0, 0, "CHROM", cell_fmt)
        ws.write(0, 1, "POS", cell_fmt)
        ws.write(0, 2, "REF", cell_fmt)
        ws.write(0, 3, "ALT", cell_fmt)
        len_anno_cols = len(anno_cols)
        for anno_idx in xrange(len_anno_cols):
            ws.write(0, anno_idx+LAYOUT_VCF_COLS, anno_cols[anno_idx], cell_fmt)
        sample_start_idx = LAYOUT_VCF_COLS + len_anno_cols
        ncol = sample_start_idx
        for sample_idx in xrange(len(samples_id)):
            sample_id = samples_id[sample_idx]
            if self.report_layout.call_detail:
                col_idx = (2*sample_idx) + (sample_start_idx)
                ws.write(0, col_idx, sample_id, cell_fmt)
                ws.write(0, col_idx+1, sample_id+"(detail)", cell_fmt)
                ncol += 2
            else:
                ws.write(0, sample_idx+sample_start_idx, sample_id, cell_fmt)
                ncol += 1
        return ncol

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
                if vcf_record.vcf_eval(cra.pattern):
                    row_color = cra.color
                    break
        if row_color is not None:
            dflt_cell_fmt = self.plain_fmts[row_color]
        else:
            dflt_cell_fmt = self.plain_fmts[NO_COLOR]
        # start writing content
        alt_allele = vcf_record.alleles[allele_idx]
        ws.write(row, 0, vcf_record.CHROM, dflt_cell_fmt)
        ws.write(row, 1, str(vcf_record.POS), dflt_cell_fmt)
        ws.write(row, 2, vcf_record.REF, dflt_cell_fmt)
        ws.write(row, 3, str(alt_allele), dflt_cell_fmt)
        len_anno_cols = len(anno_cols)
        # annotate INFO columns
        if self.report_layout.exprs is not None:
            color_col_names = map(lambda x: x.col_name,
                                  self.report_layout.exprs.actions[ACTION_COLOR_COL])
        else:
            color_col_names = []
        # for each INFO column
        for anno_idx in xrange(len_anno_cols):
            anno_col_name = anno_cols[anno_idx]
            info = vcf_record.INFO[anno_col_name]
            if (type(info) is list) and (len(info) == 1):
                info = info[0]
            elif (type(info) is list) and (len(info) > 1):
                info = info[allele_idx-1]
            if anno_col_name in PREDICTION_COLS:
                info = info.description
            if info is None:
                continue
#                info = ""
            if info == [None]:
                continue
            if info == ".":
                continue
#                info = ""
            # determine cell format
            info_cell_fmt = dflt_cell_fmt
            if anno_col_name in color_col_names:
                color_col_actions = self.report_layout.exprs.actions[ACTION_COLOR_COL]
                for cca in color_col_actions:
                    if ((anno_col_name == cca.col_name) and
                        vcf_record.vcf_eval(cca.pattern)
                        ):
                        info_cell_fmt = self.plain_fmts[cca.color]
                        break
            ws.write(row, anno_idx+LAYOUT_VCF_COLS, str(info).decode('utf-8'), info_cell_fmt)
        # annotate samples information
        sample_start_idx = LAYOUT_VCF_COLS + len_anno_cols
        for sample_idx in xrange(len(samples_id)):
            call = vcf_record.genotype(samples_id[sample_idx])
            zygo = call.cmm_gts[allele_idx]
            # determine cell(s) format
            if type(call.shared_mutations) is not list:
                zygo_fmt = dflt_cell_fmt
            elif not call.shared_mutations[allele_idx]:
                zygo_fmt = dflt_cell_fmt
            elif call.actual_gts[allele_idx] == CMMGT_HOMOZYGOTE:
                hom_shared_color = self.report_layout.cell_color_hom_shared
                zygo_fmt = self.plain_fmts[hom_shared_color]
            else:
                het_shared_color = self.report_layout.cell_color_het_shared
                zygo_fmt = self.plain_fmts[het_shared_color]
            # write content to cell(s)
            if self.report_layout.call_detail:
                col_idx = (2*sample_idx) + (sample_start_idx)
                ws.write(row, col_idx, zygo, zygo_fmt)
                formatted_call = self.__format_call_detail(call, vcf_record.FORMAT.split(":"))
                ws.write(row, col_idx+1, formatted_call, zygo_fmt)
            else:
                ws.write(row, sample_idx+sample_start_idx, zygo, zygo_fmt)

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
                        if vcf_record.vcf_eval(dra.pattern):
                            delete_row = True
                            break
                if delete_row:
                    continue
                self.__write_content(ws,
                                     row,
                                     vcf_record,
                                     allele_idx,
                                     samples_id)
                if row % RECORDS_LOG_INTERVAL == 0:
                    log_msg = str(row)
                    log_msg += " records were written to the sheet"
                    self.info(log_msg)
                row += 1
        return row

    def __set_layout(self, ws, record_size):
        ws.autofilter(0, 0, 0, record_size-1)

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
            samples = vcf_reader.samples
        else:
            samples = samples_id
        ws = self.__add_sheet(sheet_name)
        if samples_header is not None:
            ncol = self.__write_header(ws, samples_header)
        else:
            ncol = self.__write_header(ws, samples)
        if report_regions is None:
            row = self.__write_contents(ws,
                                        row,
                                        vcf_reader,
                                        samples,
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
                                            samples,
                                            check_shared)
        # freeze panes
        log_msg = "Finish .. "
        log_msg += " total of " + str(row-1)
        log_msg += " records were written to the sheet"
        self.info(log_msg)
        self.__set_layout(ws, ncol)
        ws.freeze_panes(1, 0)

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
#        if (self.report_layout.summary_families_sheet and
#            self.families_info is not None):
#            self.info("")
#            self.info(" >> add 'summary_families' sheet")
#            self.__add_muts_sheet("summary_families",
#                                  report_regions,
#                                  samples_id=self.samples_id,
#                                  samples_header=self.samples_id_w_fam_pref,
#                                  )
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
                region_params = ALL_CHROMS
            else:
                region_params = map(lambda x: x.raw_region,
                                    report_regions)
            for region_param in region_params:
                chr_slurm_log_prefix = slurm_log_prefix
                chr_slurm_log_prefix += "_chr" + region_param.split(":")[0]
                job_name = job_name_prefix
                job_name += "_chr" + region_param.split(":")[0]
                out_file = out_file_prefix
                out_file += "_chr" + region_param.split(":")[0]
                out_file += ".xlsx"
                job_params = job_params_prefix
                job_params += " -r " + region_param
                job_params += " -o " + out_file
                self.dbg(job_script + job_params)
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

    def monitor_init(self, **kwargs):
        if self.report_layout.only_summary:
            self.gen_summary_reports()
        elif self.report_layout.only_families:
            self.gen_families_reports()
        else:
            self.gen_families_reports()
            self.gen_summary_reports()
        super(MutRepPipeline, self).monitor_init(**kwargs)

#class MutRepPipeline(CMMDBPipeline):
#    """ A class to control mutation report pipeline """
#
#    def __init__(self,
#                 jobs_setup_file,
#                 ):
#        CMMDBPipeline.__init__(self,
#                               jobs_setup_file=jobs_setup_file
#                               )
#        self.__parse_report_layout()
#
#    def get_raw_repr(self):
#        return {"dataset name": self.dataset_name,
#                "project code": self.project_code,
#                "jobs report file": self.jobs_report_file,
#                }
#
#    def __parse_report_layout(self):
#        self.__report_layout = ReportLayout(self._jobs_info[JOBS_SETUP_RPT_LAYOUT_SECTION])
#
#    @property
#    def job_alloc_time(self):
#        return self._jobs_info[JOBS_SETUP_RPT_ALLOC_TIME_KEY]
#
#    @property
#    def annotated_vcf_tabix(self):
#        return self.report_layout.annotated_vcf_tabix
#
#    @property
#    def report_layout(self):
#        return self.__report_layout
#
#    @property
#    def cell_fmt_mgr(self):
#        return self.__cell_fmt_mgr
#
#    def __cal_gt(self, gt, allele_idx):
#        if gt == ".":
#            return "."
#        if gt == "./.":
#            return "."
#        gts = gt.split("/")
#        if (gts[0] == "0") and (gts[1] == "0"):
#            return "wt"
#        if (gts[0] == str(allele_idx)) or (gts[1] == str(allele_idx)):
#            if gts[0] == gts[1]:
#                return "hom"
#            else:
#                return "het"
#        return "."
#
#    def __init_cells_format(self, wb):
#        self.__cell_fmt_mgr = CellFormatManager(wb, COLOR_RGB)
#
#
#    def gen_reports(self):
#        pass
#
#    def __garbage_collecting(self):
#        pass
#
#    def monitor_action(self):
#        CMMDBPipeline.monitor_action(self)
#        self.__garbage_collecting()

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
    report_regions = get_func_arg('report_regions', kwargs)
    if report_regions is not None:
        rpt_cfg[JOBS_SETUP_RPT_REGIONS_KEY] = report_regions.split(",")
    expression_patterns = get_func_arg('expression_patterns', kwargs)
    if expression_patterns is not None:
        exprs = []
        for raw_pattern in  expression_patterns.split(";"):
            name, pattern = raw_pattern.split(":")
            expr = defaultdict(list)
            expr[JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY] = name.strip()
            expr[JOBS_SETUP_RPT_EXPRESSIONS_PATTERN_KEY] = pattern.strip()
            exprs.append(expr)
        expression_usages = get_func_arg('expression_usages', kwargs)
        if expression_usages is not None:
            for raw_usage in expression_usages.split(";"):
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
    rpt_cfg[JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY] = get_func_arg('summary_families_sheet',
                                                                kwargs,
                                                                default_val=False,
                                                                )
    rpt_cfg[JOBS_SETUP_RPT_ONLY_SUMMARY_KEY] = get_func_arg('only_summary',
                                                            kwargs,
                                                            default_val=False,
                                                            )
    rpt_cfg[JOBS_SETUP_RPT_ONLY_FAMILIES_KEY] = get_func_arg('only_families',
                                                             kwargs,
                                                             default_val=False,
                                                             )
    call_detail = get_func_arg('call_detail', kwargs, default_val=False)
    if call_detail:
        extra_anno_cols = []
        if call_detail:
            extra_anno_cols.append(JOBS_SETUP_RPT_CALL_DETAIL_KEY)
        rpt_cfg[JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY] = extra_anno_cols
    rows_filter_actions = get_func_arg('rows_filter_actions', kwargs)
    if rows_filter_actions is not None:
        filter_criterias = []
        for filter_criteria in rows_filter_actions.split(","):
            if filter_criteria == JOBS_SETUP_RPT_FILTER_RARE:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_RARE)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_INTERGENIC:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_INTERGENIC)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_INTRONIC:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_INTRONIC)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_UPSTREAM:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_UPSTREAM)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_UTR:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_UTR)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_HAS_MUTATION:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_HAS_MUTATION)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_HAS_SHARED:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_HAS_SHARED)
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
    job_setup_document[JOBS_SETUP_RPT_LAYOUT_SECTION] = rpt_cfg
    pyaml.dump(job_setup_document, stream)
