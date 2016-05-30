# This Python file uses the following encoding: utf-8
import sys
import re
import datetime
import pyaml
import yaml
import vcf
import xlsxwriter
from os import listdir
from os.path import join as join_path
from os.path import isdir
from os.path import isfile
from collections import OrderedDict
from collections import defaultdict
from pycmm.settings import PREDICTION_COLS
from pycmm.settings import DFLT_ANNOVAR_DB_FOLDER
from pycmm.settings import DFLT_ANNOVAR_DB_NAMES
from pycmm.settings import DFLT_ANNOVAR_DB_OPS
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.settings import MUTREP_FAMILY_REPORT_BIN
from pycmm.settings import MUTREP_SUMMARY_REPORT_BIN
from pycmm.template import pyCMMBase
from pycmm.utils import exec_sh
from pycmm.cmmlib.dnalib import ALL_CHROMS
from pycmm.cmmlib.taparser import TAVcfReader as VcfReader
from pycmm.cmmlib.tamodel import CMMGT_HOMOZYGOTE
from pycmm.cmmlib.tamodel import CMMGT_HETEROZYGOTE
from pycmm.flow.cmmdb import CMMDBPipeline
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ALLOC_TIME_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_LAYOUT_SECTION
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ANNO_COLS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_REGIONS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FREQ_RATIOS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXPRESSIONS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXPRESSIONS_PATTERN_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXPRESSIONS_USAGES_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXPRESSIONS_ACTION_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXPRESSIONS_INFO_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_SPLIT_CHROM_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_CALL_DETAIL_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_MT_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_RARE
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_NON_INTERGENIC
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_NON_INTRONIC
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_NON_UPSTREAM
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_NON_UTR
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_HAS_MUTATION
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FILTER_HAS_SHARED
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ONLY_SUMMARY_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ONLY_FAMILIES_KEY

ACTION_DELETE_ROW = "del_row"
ACTION_COLOR_ROW = "color_row"
ACTION_COLOR_COL = "color_col"

# *** color definition sections ***
COLOR_RGB = OrderedDict()
COLOR_RGB['GREEN_ANNIKA'] = '#CCFFCC'
COLOR_RGB['PINK_ANNIKA'] = '#E6B9B8'
COLOR_RGB['GRAY25'] = '#DCDCDC'
COLOR_RGB['PLUM'] = '#8E4585'
COLOR_RGB['GREEN'] = '#008000'
COLOR_RGB['LIGHT_GREEN'] = '#90EE90'
COLOR_RGB['SKY_BLUE'] = '#87CEEB'
COLOR_RGB['GRAY40'] = '#808080'
COLOR_RGB['LIGHT_YELLOW'] = '#FFFFE0'
COLOR_RGB['OLIVE'] = '#808000'
COLOR_RGB['ORANGE'] = '#FF6600'
COLOR_RGB['DARK_SLATE_GRAY'] = '#2F4F4F'
COLOR_RGB['PURPLE'] = '#800080'
COLOR_RGB['RED'] = '#FF0000'
COLOR_RGB['ROSY_BROWN'] = '#BC8F8F'
COLOR_RGB['SILVER'] = '#C0C0C0'
COLOR_RGB['SKY_BLUE'] = '#87CEEB'
COLOR_RGB['TAN'] = '#D2B48C'
COLOR_RGB['TEAL'] = '#008080'
COLOR_RGB['TURQUOISE'] = '#40E0D0'
COLOR_RGB['YELLOW'] = '#FFFF00'
COLOR_RGB['MEDIUM_AQUA_MARINE'] = '#66CDAA'
COLOR_RGB['BLUE'] = '#0000FF'
COLOR_RGB['SLATE_GRAY'] = '#708090'
COLOR_RGB['LIME_GREEN'] = '#32CD32'
COLOR_RGB['BROWN'] = '#800000'
COLOR_RGB['CORAL'] = '#FF7F50'
COLOR_RGB['DARK_BLUE'] = '#00008B'
COLOR_RGB['YELLOW_GREEN'] = '#9ACD32'
COLOR_RGB['DODGER_BLUE'] = '#1E90FF'
COLOR_RGB['GOLDEN_ROD'] = '#DAA520'
COLOR_RGB['CYAN'] = '#00FFFF'
COLOR_RGB['ROYAL_BLUE'] = '#4169E1'
COLOR_RGB['GOLD'] = '#FFD700'
COLOR_RGB['LIME'] = '#00FF00'
COLOR_RGB['MAGENTA'] = '#FF00FF'
COLOR_RGB['ICEBLUE'] = '#A5F2F3'
COLOR_RGB['LIGHT_BLUE'] = '#ADD8E6'

DFLT_COLOR_HET_SHARED = 'SILVER'
DFLT_COLOR_HOM_SHARED = 'GRAY25'
CELL_TYPE_HET_SHARED = 'HET_SHARED'
CELL_TYPE_HOM_SHARED = 'HOM_SHARED'

DFLT_FMT = 'default_format'
# *********************************

CHROM_POS_PATTERN = re.compile(r'''(?P<chrom>.+?):(?P<start_pos>.+?)-(?P<end_pos>.+)''')
LAYOUT_VCF_COLS = 4
RECORDS_LOG_INTERVAL = 1000

class CellFormatManager(pyCMMBase):
    """ A class to define cell format property """

    def __init__(self, work_book, color_dict):
        self.__wb = work_book
        self.__color_dict = color_dict
        self.__dflt_hash_fmt = {'font_name': 'Arial', 'font_size': 10}
        self.__dflt_fmt = self.__add_fmt(self.__dflt_hash_fmt)
        self.__init_colors_formats()

    def get_raw_repr(self):
        return {"color dict": self.__color_dict,
                "number of color": self.n_colors,
                }

    def __add_fmt(self, fmt_dict):
        return self.__wb.add_format(fmt_dict)

    @property
    def default_format(self):
        return self.__dflt_fmt

    def __init_colors_format(self, default_hash_format):
        fmts = OrderedDict()
        fmts[DFLT_FMT] = self.__add_fmt(default_hash_format)
        colors = self.__color_dict.keys()
        for color_idx in xrange(len(colors)):
            color = colors[color_idx]
            color_hash_fmt = default_hash_format.copy()
            color_hash_fmt['bg_color'] = self.__color_dict[color]
            fmts[color_idx] = self.__add_fmt(color_hash_fmt)
            fmts[color] = fmts[color_idx]
        return fmts

    def __init_colors_formats(self):
        dflt_cell_hash_fmt = self.__dflt_hash_fmt.copy()
        self.__cell_fmts = self.__init_colors_format(dflt_cell_hash_fmt)

    @property
    def n_colors(self):
        return len(self.__color_dict)

    @property
    def cell_fmts(self):
        return self.__cell_fmts

class ReportRegion(pyCMMBase):
    """ A structure to parse and keep mutation report region """

    def __init__(self,
                 raw_region,
                 ):
        self.__parse_region(raw_region)
        self.__raw_region = raw_region

    def get_raw_repr(self):
        if self.start_pos is not None:
            return {"chrom": self.chrom,
                    "start position": self.start_pos,
                    "end_pos": self.end_pos,
                    }
        else:
            return {"chrom": self.chrom,
                    }

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
                 ):
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
                 ):
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
                 ):
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
                 ):
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

class ReportLayout(pyCMMBase):
    """ A structure to parse and keep mutation report layout """

    def __init__(self,
                 layout_params,
                 ):
        self.__layout_params = layout_params
        self.__freq_ratios = None
        self.__anno_cols = None
        self.__exprs = self.__parse_exprs()
        self.__init_cell_colors()

    def __cal_anno_cols(self):
        # open the annotated vcf tabix file for checking
        # if the expected output columns are exist
        vcf_reader = VcfReader(filename=self.annotated_vcf_tabix)
        vcf_record = vcf_reader.next()
        # generate columns list
        anno_cols = []
        for col_name in self.__layout_params[JOBS_SETUP_RPT_ANNO_COLS_KEY]:
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

    @property
    def anno_cols(self):
        if self.__anno_cols is None:
            self.__anno_cols = self.__cal_anno_cols()
        return self.__anno_cols

    @property
    def anno_excl_tags(self):
        if JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY in self.__layout_params:
            return self.__layout_params[JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY]
        else:
            return []

    @property
    def annotated_vcf_tabix(self):
        if JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX in self.__layout_params:
            return self.__layout_params[JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX]
        else:
            return None

    @property
    def report_regions(self):
        if JOBS_SETUP_RPT_REGIONS_KEY in self.__layout_params:
            return map(lambda x: ReportRegion(str(x)),
                       self.__layout_params[JOBS_SETUP_RPT_REGIONS_KEY])
        else:
            return None

    @property
    def freq_ratios(self):
        if ((self.__freq_ratios is None) and
            (JOBS_SETUP_RPT_FREQ_RATIOS_KEY in self.__layout_params)
            ):
            freq_ratios = OrderedDict()
            for freq_ratio in self.__layout_params[JOBS_SETUP_RPT_FREQ_RATIOS_KEY]:
                col = freq_ratio[JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY]
                freq = freq_ratio[JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY]
                freq_ratios[col] = freq
            self.__freq_ratios = freq_ratios
        elif self.__freq_ratios is None:
            self.__freq_ratios = []
        return self.__freq_ratios

    @property
    def exprs(self):
        return self.__exprs

    def __parse_exprs(self):
        if JOBS_SETUP_RPT_EXPRESSIONS_KEY in self.__layout_params:
            return VcfExpressions(self.__layout_params[JOBS_SETUP_RPT_EXPRESSIONS_KEY])
        else:
            return None

    @property
    def split_chrom(self):
        if JOBS_SETUP_RPT_SPLIT_CHROM_KEY in self.__layout_params:
            return self.__layout_params[JOBS_SETUP_RPT_SPLIT_CHROM_KEY]
        else:
            return False

    @property
    def summary_families_sheet(self):
        if JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY in self.__layout_params:
            return self.__layout_params[JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY]
        else:
            return False

    @property
    def call_detail(self):
        if JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_CALL_DETAIL_KEY in self.__layout_params[JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY]

    @property
    def filter_rare(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_RARE in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def filter_non_intergenic(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_NON_INTERGENIC in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def filter_non_intronic(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_NON_INTRONIC in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def filter_non_upstream(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_NON_UPSTREAM in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def filter_non_downtream(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def filter_non_utr(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_NON_UTR in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def filter_non_synonymous(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def filter_has_mutation(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_HAS_MUTATION in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def filter_has_shared(self):
        if JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_FILTER_HAS_SHARED in self.__layout_params[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY]

    @property
    def only_summary(self):
        if JOBS_SETUP_RPT_ONLY_SUMMARY_KEY in self.__layout_params:
            return self.__layout_params[JOBS_SETUP_RPT_ONLY_SUMMARY_KEY]
        else:
            return False

    @property
    def only_families(self):
        if JOBS_SETUP_RPT_ONLY_FAMILIES_KEY in self.__layout_params:
            return self.__layout_params[JOBS_SETUP_RPT_ONLY_FAMILIES_KEY]
        else:
            return False

class MutRepPipeline(CMMDBPipeline):
    """ A class to control mutation report pipeline """

    def __init__(self,
                 jobs_setup_file,
                 ):
        CMMDBPipeline.__init__(self,
                               jobs_setup_file=jobs_setup_file
                               )
        self.__parse_report_layout()

    def get_raw_repr(self):
        return {"dataset name": self.dataset_name,
                "project code": self.project_code,
                "jobs report file": self.jobs_report_file,
                }

    def __parse_report_layout(self):
        self.__report_layout = ReportLayout(self._jobs_info[JOBS_SETUP_RPT_LAYOUT_SECTION])

    @property
    def rpt_alloc_time(self):
        return self._jobs_info[JOBS_SETUP_RPT_ALLOC_TIME_KEY]

    @property
    def annotated_vcf_tabix(self):
        return self.report_layout.annotated_vcf_tabix

    @property
    def report_layout(self):
        return self.__report_layout

    @property
    def summary_rpt_file(self):
        return join_path(self.rpts_out_dir,
                         self.dataset_name+"_summary.xlsx")

    @property
    def cell_fmt_mgr(self):
        return self.__cell_fmt_mgr

    def __cal_gt(self, gt, allele_idx):
        if gt == ".":
            return "."
        if gt == "./.":
            return "."
        gts = gt.split("/")
        if (gts[0] == "0") and (gts[1] == "0"):
            return "wt"
        if (gts[0] == str(allele_idx)) or (gts[1] == str(allele_idx)):
            if gts[0] == gts[1]:
                return "hom"
            else:
                return "het"
        return "."

    def __set_layout(self, ws, record_size):
        ws.autofilter(0, 0, 0, record_size-1)

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

    def __init_cells_format(self, wb):
        self.__cell_fmt_mgr = CellFormatManager(wb, COLOR_RGB)

    def __add_sheet(self, wb, sheet_name):
        ws = wb.add_worksheet(sheet_name)
        ws.set_default_row(12)
        return ws

    def __write_header(self, ws, samples_list):
        anno_cols = self.report_layout.anno_cols
        cell_fmt = self.cell_fmt_mgr.cell_fmts[DFLT_FMT]
        ws.write(0, 0, "CHROM", cell_fmt)
        ws.write(0, 1, "POS", cell_fmt)
        ws.write(0, 2, "REF", cell_fmt)
        ws.write(0, 3, "ALT", cell_fmt)
        len_anno_cols = len(anno_cols)
        for anno_idx in xrange(len_anno_cols):
            ws.write(0, anno_idx+LAYOUT_VCF_COLS, anno_cols[anno_idx], cell_fmt)
        sample_start_idx = LAYOUT_VCF_COLS + len_anno_cols
        ncol = sample_start_idx
        for sample_idx in xrange(len(samples_list)):
            sample_id = samples_list[sample_idx]
            if self.report_layout.call_detail:
                col_idx = (2*sample_idx) + (sample_start_idx)
                ws.write(0, col_idx, sample_id, cell_fmt)
                ws.write(0, col_idx+1, sample_id+"(detail)", cell_fmt)
                ncol += 2
            else:
                ws.write(0, sample_idx+sample_start_idx, sample_id, cell_fmt)
                ncol += 1
        return ncol

    def __write_content(self,
                        ws,
                        row,
                        vcf_record,
                        allele_idx,
                        samples_list,
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
            dflt_cell_fmt = self.cell_fmt_mgr.cell_fmts[row_color]
        else:
            dflt_cell_fmt = self.cell_fmt_mgr.cell_fmts[DFLT_FMT]
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
                info = ""
            if info == [None]:
                info = ""
            # determine cell format
            info_cell_fmt = dflt_cell_fmt
            if anno_col_name in color_col_names:
                color_col_actions = self.report_layout.exprs.actions[ACTION_COLOR_COL]
                for cca in color_col_actions:
                    if ((anno_col_name == cca.col_name) and
                        vcf_record.vcf_eval(cca.pattern)
                        ):
                        info_cell_fmt = self.cell_fmt_mgr.cell_fmts[cca.color]
                        break
            ws.write(row, anno_idx+LAYOUT_VCF_COLS, str(info).decode('utf-8'), info_cell_fmt)
        # annotate samples information
        sample_start_idx = LAYOUT_VCF_COLS + len_anno_cols
        for sample_idx in xrange(len(samples_list)):
            call = vcf_record.genotype(samples_list[sample_idx])
            zygo = call.cmm_gts[allele_idx]
            # determine cell(s) format
            if type(call.shared_mutations) is not list:
                zygo_fmt = dflt_cell_fmt
            elif not call.shared_mutations[allele_idx]:
                zygo_fmt = dflt_cell_fmt
            elif call.actual_gts[allele_idx] == CMMGT_HOMOZYGOTE:
                hom_shared_color = self.report_layout.cell_color_hom_shared
                zygo_fmt = self.cell_fmt_mgr.cell_fmts[hom_shared_color]
            else:
                het_shared_color = self.report_layout.cell_color_het_shared
                zygo_fmt = self.cell_fmt_mgr.cell_fmts[het_shared_color]
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
                         samples_list,
                         check_shared,
                         ):
        for vcf_record in vcf_records:
            for allele_idx in xrange(1, len(vcf_record.alleles)):
                if (check_shared and
                    not vcf_record.is_shared(samples_list, allele_idx)):
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
                    not vcf_record.has_mutation(samples_list, allele_idx)):
                    continue
                if (self.report_layout.filter_has_shared and
                    not vcf_record.has_shared(allele_idx)):
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
                                     samples_list)
                if row % RECORDS_LOG_INTERVAL == 0:
                    log_msg = str(row)
                    log_msg += " records were written to the sheet"
                    self.info(log_msg)
                row += 1
        return row

    def __add_muts_sheet(self,
                         wb,
                         sheet_name,
                         report_regions,
                         samples_list=None,
                         samples_header=None,
                         check_shared=False,
                         ):
        if self.annotated_vcf_tabix.endswith('.vcf.gz'):
            vcf_reader = VcfReader(filename=self.annotated_vcf_tabix,
                                   family_infos=self.family_infos,
                                   freq_ratios=self.report_layout.freq_ratios,
                                   )
        else:
            self.thrown(self.annotated_vcf_tabix + ' does not endswith .vcf.gz')
        row = 1
        if samples_list is None:
            samples = vcf_reader.samples
        else:
            samples = samples_list
        ws = self.__add_sheet(wb, sheet_name)
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
                slurm_log_file = slurm_log_prefix
                slurm_log_file += "_chr" + region_param.split(":")[0]
                slurm_log_file += "_" + self.time_stamp.strftime("%Y%m%d%H%M%S")
                slurm_log_file += ".log"
                job_name = job_name_prefix
                job_name += "_chr" + region_param.split(":")[0]
                out_file = out_file_prefix
                out_file += "_chr" + region_param.split(":")[0]
                out_file += ".xlsx"
                job_params = job_params_prefix
                job_params += " -r " + region_param
                job_params += " -o " + out_file
                self.dbg(job_script + job_params)
# *********************************************************************************************** Need refactoring ***********************************************************************************************
                self.submit_job(job_name,
                                self.project_code,
                                "core",
                                "1",
                                self.rpt_alloc_time,
                                slurm_log_file,
                                job_script,
                                job_params,
                                )
# *********************************************************************************************** Need refactoring ***********************************************************************************************
        else:
            slurm_log_file = slurm_log_prefix
            slurm_log_file += "_" + self.time_stamp.strftime("%Y%m%d%H%M%S")
            slurm_log_file += ".log"
            job_name = job_name_prefix
            out_file = out_file_prefix + ".xlsx"
            job_params = job_params_prefix
            if report_regions is not None:
                job_params += " -r " + ",".join(map(lambda x: x.raw_region,
                                                    report_regions))
            job_params += " -o " + out_file
# *********************************************************************************************** Need refactoring ***********************************************************************************************
            self.submit_job(job_name,
                            self.project_code,
                            "core",
                            "1",
                            self.rpt_alloc_time,
                            slurm_log_file,
                            job_script,
                            job_params,
                            )
# *********************************************************************************************** Need refactoring ***********************************************************************************************

    def gen_summary_report(self,
                           report_regions,
                           out_file=None,
                           ):
        if out_file is None:
            out_file = self.summary_rpt_file
        self.info("")
        self.info(" >>>> generating summary report")
        self.info(" >>>> report file: " + out_file)
        wb = xlsxwriter.Workbook(out_file)
        self.__init_cells_format(wb)
        self.info("")
        self.info(" >> add 'summary_all' sheet")
        self.__add_muts_sheet(wb,
                              "summary_all",
                              report_regions,
                              )
        if (self.report_layout.summary_families_sheet and
            self.family_infos is not None):
            self.info("")
            self.info(" >> add 'summary_families' sheet")
            self.__add_muts_sheet(wb,
                                  "summary_families",
                                  report_regions,
                                  samples_list=self.samples_list,
                                  samples_header=self.samples_list_w_fam_pref,
                                  )
        wb.close()

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
        fam_info = self.family_infos[fam_id]
        if out_file is None:
            out_file = join_path(self.rpts_out_dir,
                                 self.dataset_name+"_fam"+fam_info.fam_id+".xlsx")
        self.info("")
        self.info(" >>>> generating family report for family " + fam_info.fam_id)
        self.info(" >>>> report file: " + out_file)
        samples_list = map(lambda x: x.sample_id, fam_info.members)
        wb = xlsxwriter.Workbook(out_file)
        self.__init_cells_format(wb)
        if len(samples_list) == 1:
            self.info("")
            self.info(" >> add '" + samples_list[0] + "' sheet")
            self.__add_muts_sheet(wb,
                                  samples_list[0],
                                  report_regions,
                                  samples_list=samples_list,
                                  check_shared=True,
                                  )
        else:
            self.info("")
            self.info(" >> add 'shared' sheet")
            self.__add_muts_sheet(wb,
                                  "shared",
                                  report_regions,
                                  samples_list=samples_list,
                                  check_shared=True)
            for sample_name in samples_list:
                self.info("")
                self.info(" >> add '" + sample_name + "' sheet")
                self.__add_muts_sheet(wb,
                                      sample_name,
                                      report_regions,
                                      samples_list=[sample_name],
                                      check_shared=True)
        wb.close()

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
        if self.family_infos is None:
            return
        for fam_id in self.family_infos:
            self.__gen_family_reports(fam_id)

    def gen_reports(self):
        pass

    def __garbage_collecting(self):
        pass

    def monitor_action(self):
        CMMDBPipeline.monitor_action(self)
        self.__garbage_collecting()
