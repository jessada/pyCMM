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
from pycmm.proc.tavcf import TableAnnovarVCFReader as VCFReader
from pycmm.template import pyCMMBase
from pycmm.settings import DFLT_ANNOVAR_DB_FOLDER
from pycmm.settings import DFLT_ANNOVAR_DB_NAMES
from pycmm.settings import DFLT_ANNOVAR_DB_OPS
from pycmm.utils import exec_sh
from pycmm.utils import mylogger
from pycmm.flow.cmmdb import CMMDBPipeline
from pycmm.settings import MUTREP_ALLOC_TIME
from pycmm.flow.cmmdb import JOBS_SETUP_REPORT_LAYOUT_SECTION
from pycmm.flow.cmmdb import JOBS_SETUP_REPORT_ANNO_COLS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_REPORT_REGIONS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_REPORT_CALL_INFO_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_REPORT_FREQ_RATIOS_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_REPORT_FREQ_RATIOS_COL_KEY
from pycmm.flow.cmmdb import JOBS_SETUP_REPORT_FREQ_RATIOS_FREQ_KEY


CHROM_POS_PATTERN = re.compile(r'''(?P<chrom>.+?):(?P<start_pos>.+?)-(?P<end_pos>.+)''')
LAYOUT_VCF_COLS = 4

class ReportRegion(pyCMMBase):
    """ A structure to parse and keep mutation report layout """

    def __init__(self,
                 raw_region,
                 ):
        self.__parse_region(raw_region)

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

    @property
    def anno_cols(self):
        return self.__layout_params[JOBS_SETUP_REPORT_ANNO_COLS_KEY]

    @property
    def report_regions(self):
        if JOBS_SETUP_REPORT_REGIONS_KEY in self.__layout_params:
            return map(lambda x: ReportRegion(str(x)),
                       self.__layout_params[JOBS_SETUP_REPORT_REGIONS_KEY])
        else:
            return None

    @property
    def call_info(self):
        return self.__layout_params[JOBS_SETUP_REPORT_CALL_INFO_KEY]

    @property
    def freq_ratios(self):
        if JOBS_SETUP_REPORT_FREQ_RATIOS_KEY in self.__layout_params:
            freq_ratios = {}
            for freq_ratio in self.__layout_params[JOBS_SETUP_REPORT_FREQ_RATIOS_KEY]:
                col = freq_ratio[JOBS_SETUP_REPORT_FREQ_RATIOS_COL_KEY]
                freq = freq_ratio[JOBS_SETUP_REPORT_FREQ_RATIOS_FREQ_KEY]
                freq_ratios[col] = freq
            return freq_ratios
        else:
            return None

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
        self.__report_layout = ReportLayout(self._jobs_info[JOBS_SETUP_REPORT_LAYOUT_SECTION])

    @property
    def report_layout(self):
        return self.__report_layout

    @property
    def summary_rpt_file(self):
        return join_path(self.rpts_out_dir,
                         self.dataset_name+"_summary.xlsx")

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

    def __format_call_info(self, call, call_fmt):
        call_info = []
        for fmt_idx in xrange(len(call_fmt)):
            fmt = call_fmt[fmt_idx]
            data = call.data[fmt_idx]
            if type(data) is list:
                data = map(lambda x: str(x), data)
                call_info.append(fmt + "=" + ",".join(data))
            else:
                if data is None:
                    data = "."
                call_info.append(fmt + "=" + str(data))
        return ":".join(call_info)

    def __add_sheet(self, sheet_name):
        ws = self.__wb.add_worksheet(sheet_name)
        ws.set_default_row(12)
        return ws

    def __write_header(self, ws, samples_list):
        anno_cols = self.report_layout.anno_cols
        ws.write(0, 0, "CHROM")
        ws.write(0, 1, "POS")
        ws.write(0, 2, "REF")
        ws.write(0, 3, "ALT")
        len_anno_cols = len(anno_cols)
        for anno_idx in xrange(len_anno_cols):
            ws.write(0, anno_idx+LAYOUT_VCF_COLS, anno_cols[anno_idx])
        sample_start_idx = LAYOUT_VCF_COLS + len_anno_cols
        for sample_idx in xrange(len(samples_list)):
            sample_id = samples_list[sample_idx]
            if self.report_layout.call_info:
                col_idx = (2*sample_idx) + (sample_start_idx)
                ws.write(0, col_idx, sample_id)
                ws.write(0, col_idx+1, sample_id+"(call info)")
            else:
                ws.write(0, sample_idx+sample_start_idx, sample_id)

    def __write_content(self,
                        ws,
                        row,
                        vcf_record,
                        allele_idx,
                        samples_list,
                        ):
        anno_cols = self.report_layout.anno_cols
        alt_allele = vcf_record.alleles[allele_idx]
        ws.write(row, 0, vcf_record.CHROM)
        ws.write(row, 1, vcf_record.POS)
        ws.write(row, 2, vcf_record.REF)
        ws.write(row, 3, str(alt_allele))
        len_anno_cols = len(anno_cols)
        for anno_idx in xrange(len_anno_cols):
            info = vcf_record.INFO[anno_cols[anno_idx]]
            if (type(info) is list) and (len(info) == 1):
                info = info[0]
            elif (type(info) is list) and (len(info) > 1):
                info = info[allele_idx-1]
            if info is None:
                info = ""
            if info == [None]:
                info = ""
            ws.write(row, anno_idx+LAYOUT_VCF_COLS, str(info))
        sample_start_idx = LAYOUT_VCF_COLS + len_anno_cols
        for sample_idx in xrange(len(samples_list)):
            call = vcf_record.genotype(samples_list[sample_idx])
            zygo = self.__cal_gt(call["GT"], allele_idx)
            if self.report_layout.call_info:
                col_idx = (2*sample_idx) + (sample_start_idx)
                ws.write(row, col_idx, zygo)
                formatted_call = self.__format_call_info(call, vcf_record.FORMAT.split(":"))
                ws.write(row, col_idx+1, formatted_call)
            else:
                ws.write(row, sample_idx+sample_start_idx, zygo)

    def __add_muts_sheet(self, sheet_name, samples_list=None):
        if self.input_vcf_tabix.endswith('.vcf.gz'):
            vcf_reader = VCFReader(filename=self.input_vcf_tabix)
        else:
            self.thrown(self.input_vcf_tabix + ' does not endswith .vcf.gz')
        row = 1
        if samples_list is None:
            samples_list = vcf_reader.samples
        ws = self.__add_sheet(sheet_name)
        self.__write_header(ws, samples_list)
        if self.report_layout.report_regions is None:
            for vcf_record in vcf_reader:
                for allele_idx in xrange(1, len(vcf_record.alleles)):
                    self.__write_content(ws, row, vcf_record, allele_idx, samples_list)
                    row += 1
        else:
            for report_region in self.report_layout.report_regions:
                if report_region.start_pos is None:
                    vcf_records = vcf_reader.fetch(report_region.chrom)
                else:
                    vcf_records = vcf_reader.fetch(report_region.chrom,
                                                   int(report_region.start_pos),
                                                   int(report_region.end_pos),
                                                   )
                for vcf_record in vcf_records:
                    for allele_idx in xrange(1, len(vcf_record.alleles)):
                        self.__write_content(ws, row, vcf_record, allele_idx, samples_list)
                        row += 1
        # freeze panes
        ws.freeze_panes(1, 0)

    def gen_summary_report(self, out_file=None):
        if out_file is None:
            out_file = self.summary_rpt_file
        self.__wb = xlsxwriter.Workbook(out_file)
        self.__add_muts_sheet("summary")
        self.__wb.close()

    def gen_family_reports(self,
                           family_info,
                           out_file=None
                           ):
        self.__wb = xlsxwriter.Workbook(out_file)
        if self.families_list is None:
            self.__add_muts_sheet("summary")
        else:
            pass
        self.__wb.close()

    def __garbage_collecting(self):
        pass

    def monitor_action(self):
        CMMDBPipeline.monitor_action(self)
        self.__garbage_collecting()
