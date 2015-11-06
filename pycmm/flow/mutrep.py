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
from pycmm.settings import PREDICTION_COLS
from pycmm.settings import DFLT_ANNOVAR_DB_FOLDER
from pycmm.settings import DFLT_ANNOVAR_DB_NAMES
from pycmm.settings import DFLT_ANNOVAR_DB_OPS
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.settings import MUTREP_FAMILY_REPORT_BIN
from pycmm.settings import MUTREP_SUMMARY_REPORT_BIN
#from pycmm.settings import MT_ANNO_COLS
from pycmm.template import pyCMMBase
from pycmm.utils import exec_sh
from pycmm.utils import mylogger
from pycmm.proc.taparser import TAVcfReader as VcfReader
from pycmm.flow.cmmdb import CMMDBPipeline
from pycmm.flow.cmmdb import ALL_CHROMS
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ALLOC_TIME_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_LAYOUT_SECTION
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ANNO_COLS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_REGIONS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FREQ_RATIOS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_SPLIT_CHROM_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_CALL_DETAIL_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_MT_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ROW_EXCLUSION_CRITERIA_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXCLUDE_COMMON
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXCLUDE_INTERGENIC
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_EXCLUDE_INTRONIC
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ONLY_SUMMARY_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_RPT_ONLY_FAMILIES_KEY

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

class ReportLayout(pyCMMBase):
    """ A structure to parse and keep mutation report layout """

    def __init__(self,
                 layout_params,
                 ):
        self.__layout_params = layout_params
        self.__anno_cols = self.__cal_anno_cols()

    def __cal_anno_cols(self):
        anno_cols = []
        for col_name in self.__layout_params[JOBS_SETUP_RPT_ANNO_COLS_KEY]:
            excluded = False
            for excl_tag in self.anno_excl_tags:
                if excl_tag in ALL_MUTREP_ANNO_COLS[col_name]:
                    excluded = True
                    break
            if not excluded:
                anno_cols.append(col_name)
        return anno_cols

    @property
    def anno_cols(self):
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
        if JOBS_SETUP_RPT_FREQ_RATIOS_KEY in self.__layout_params:
            freq_ratios = {}
            for freq_ratio in self.__layout_params[JOBS_SETUP_RPT_FREQ_RATIOS_KEY]:
                col = freq_ratio[JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY]
                freq = freq_ratio[JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY]
                freq_ratios[col] = freq
            return freq_ratios
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
    def exclude_common(self):
        if JOBS_SETUP_RPT_ROW_EXCLUSION_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_EXCLUDE_COMMON in self.__layout_params[JOBS_SETUP_RPT_ROW_EXCLUSION_CRITERIA_KEY]

    @property
    def exclude_intergenic(self):
        if JOBS_SETUP_RPT_ROW_EXCLUSION_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_EXCLUDE_INTERGENIC in self.__layout_params[JOBS_SETUP_RPT_ROW_EXCLUSION_CRITERIA_KEY]

    @property
    def exclude_intronic(self):
        if JOBS_SETUP_RPT_ROW_EXCLUSION_CRITERIA_KEY not in self.__layout_params:
            return False
        else:
            return JOBS_SETUP_RPT_EXCLUDE_INTRONIC in self.__layout_params[JOBS_SETUP_RPT_ROW_EXCLUSION_CRITERIA_KEY]

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
        mylogger.getLogger(__name__)
        CMMDBPipeline.__init__(self,
                               jobs_setup_file=jobs_setup_file
                               )
        self.__parse_report_layout()
        mylogger.getLogger(__name__)

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
        if self.report_layout.annotated_vcf_tabix is not None:
            return self.report_layout.annotated_vcf_tabix
        else:
            return self.annovar_config.annotated_vcf + ".gz"

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
        cell_fmt = self.cell_fmt_mgr.cell_fmts[DFLT_FMT]
        alt_allele = vcf_record.alleles[allele_idx]
        ws.write(row, 0, vcf_record.CHROM, cell_fmt)
        ws.write(row, 1, str(vcf_record.POS), cell_fmt)
        ws.write(row, 2, vcf_record.REF, cell_fmt)
        ws.write(row, 3, str(alt_allele), cell_fmt)
        len_anno_cols = len(anno_cols)
        # annotate INFO columns
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
            ws.write(row, anno_idx+LAYOUT_VCF_COLS, str(info).decode('utf-8'), cell_fmt)
        # annotate samples information
        sample_start_idx = LAYOUT_VCF_COLS + len_anno_cols
        for sample_idx in xrange(len(samples_list)):
            call = vcf_record.genotype(samples_list[sample_idx])
            zygo = call.cmm_gts[allele_idx]
            if self.report_layout.call_detail:
                col_idx = (2*sample_idx) + (sample_start_idx)
                ws.write(row, col_idx, zygo, cell_fmt)
                formatted_call = self.__format_call_detail(call, vcf_record.FORMAT.split(":"))
                ws.write(row, col_idx+1, formatted_call, cell_fmt)
            else:
                ws.write(row, sample_idx+sample_start_idx, zygo, cell_fmt)

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
                if (self.report_layout.exclude_common and
                    not vcf_record.is_rare(self.report_layout.freq_ratios,
                                           allele_idx=allele_idx)):
                    continue
                if (self.report_layout.exclude_intergenic and
                    vcf_record.is_intergenic[allele_idx]):
                    continue
                if (self.report_layout.exclude_intronic and
                    vcf_record.is_intronic[allele_idx]):
                    continue
                self.__write_content(ws,
                                     row,
                                     vcf_record,
                                     allele_idx,
                                     samples_list)
                if row % RECORDS_LOG_INTERVAL == 0:
                    log_msg = str(row)
                    log_msg += " records were written to the sheet"
                    mylogger.info(log_msg)
                row += 1
        return row

    def __add_muts_sheet(self,
                         wb,
                         sheet_name,
                         report_regions,
                         samples_list=None,
                         check_shared=False,
                         ):
        if self.annotated_vcf_tabix.endswith('.vcf.gz'):
            vcf_reader = VcfReader(filename=self.annotated_vcf_tabix)
        else:
            self.thrown(self.annotated_vcf_tabix + ' does not endswith .vcf.gz')
        row = 1
        if samples_list is None:
            samples = vcf_reader.samples
        else:
            samples = samples_list
        ws = self.__add_sheet(wb, sheet_name)
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
        mylogger.info(log_msg)
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
                mylogger.debug(job_script + job_params)
                self.submit_job(job_name,
                                self.project_code,
                                "core",
                                "1",
                                self.rpt_alloc_time,
                                slurm_log_file,
                                job_script,
                                job_params,
                                )
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
            self.submit_job(job_name,
                            self.project_code,
                            "core",
                            "1",
                            self.rpt_alloc_time,
                            slurm_log_file,
                            job_script,
                            job_params,
                            )

    def gen_summary_report(self,
                           report_regions,
                           out_file=None,
                           ):
        if out_file is None:
            out_file = self.summary_rpt_file
        mylogger.info("")
        mylogger.info(" >>>> generating summary report")
        mylogger.info(" >>>> report file: " + out_file)
        wb = xlsxwriter.Workbook(out_file)
        self.__init_cells_format(wb)
        mylogger.info("")
        mylogger.info(" >> add 'summary_all' sheet")
        self.__add_muts_sheet(wb,
                              "summary_all",
                              report_regions,
                              )
        if (self.report_layout.summary_families_sheet and
            self.family_infos is not None):
            mylogger.info("")
            mylogger.info(" >> add 'summary_families' sheet")
            self.__add_muts_sheet(wb,
                                  "summary_families",
                                  report_regions,
                                  samples_list=self.samples_list,
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
        mylogger.info("")
        mylogger.info(" >>>> generating family report for family " + fam_info.fam_id)
        mylogger.info(" >>>> report file: " + out_file)
        samples_list = map(lambda x: x.sample_id, fam_info.members)
        wb = xlsxwriter.Workbook(out_file)
        self.__init_cells_format(wb)
        if len(samples_list) == 1:
            mylogger.info("")
            mylogger.info(" >> add '" + samples_list[0] + "' sheet")
            self.__add_muts_sheet(wb,
                                  samples_list[0],
                                  report_regions,
                                  samples_list=samples_list,
                                  check_shared=True,
                                  )
        else:
            mylogger.info("")
            mylogger.info(" >> add 'shared' sheet")
            self.__add_muts_sheet(wb,
                                  "shared",
                                  report_regions,
                                  samples_list=samples_list,
                                  check_shared=True)
            for sample_name in samples_list:
                mylogger.info("")
                mylogger.info(" >> add '" + sample_name + "' sheet")
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
