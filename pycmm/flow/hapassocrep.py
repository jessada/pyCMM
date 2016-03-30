# This Python file uses the following encoding: utf-8
#import xlsxwriter
from os.path import join as join_path
from pycmm.template import pyCMMBase
from pycmm.cmmlib.xlslib import CMMWorkbook as Workbook
from pycmm.cmmlib.xlslib import NO_COLOR
from pycmm.flow import CMMPipeline
from pycmm.flow.plink import CMMPipeline
from pycmm.flow.plink import JOBS_SETUP_REPORT_PARAMS_SECTION
from pycmm.flow.plink import JOBS_SETUP_CUTOFF_PVALUE_KEY
from pycmm.cmmlib.plinklib import SnpInfoReader
from pycmm.cmmlib.plinklib import HapAssocReader

HAPLO_INFO_SIZE = 4
SNP_INFO_SIZE = 5

COL_SNP_IDX = 0
COL_F_MISS_A_IDX = 1
COL_F_MISS_U_IDX = 2
COL_CHROM_IDX = 3
COL_POS_IDX = 4


class RptParams(pyCMMBase):
    """  To handle and parse PLINK parameters  """

    def __init__(self, params, **kwargs):
        self.__params = params
        super(RptParams, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["cutoff p-value"] = self.cutoff_pvalue
        return raw_repr

    @property
    def cutoff_pvalue(self):
        if JOBS_SETUP_CUTOFF_PVALUE_KEY in self.__params:
            return self.__params[JOBS_SETUP_CUTOFF_PVALUE_KEY]
        return DFLT_CUTOFF_PVALUE

class HapAssocRepPipeline(CMMPipeline):
    """ To control production of haplotype association study reports """

    def __init__(self,
                 hap_assoc_file,
                 snp_info_file,
                 **kwargs
                 ):
        self.__hap_assoc_file = hap_assoc_file
        self.__load_snp_info(snp_info_file)
        self.__n_haplos = None
        self.__rpt_params = None
        super(HapAssocRepPipeline, self).__init__(**kwargs)

    @property
    def rpt_params(self):
        if self.__rpt_params is None:
            self.__rpt_params = RptParams(self._jobs_info[JOBS_SETUP_REPORT_PARAMS_SECTION])
        return self.__rpt_params

    @property
    def snp_info_dict(self):
        return self.__snp_info_dict

    @property
    def plain_fmts(self):
        return self.__plain_fmts

    @property
    def snp_fmts(self):
        return self.__plain_fmts

    @property
    def bp_fmts(self):
        return self.__bp_fmts

    @property
    def haplo_fmts(self):
        return self.__haplo_fmts

    @property
    def n_snps(self):
        return self.__n_snps

    @property
    def n_haplos(self):
        if self.__n_haplos is None:
            hap_assoc_reader = HapAssocReader(file_name=self.__hap_assoc_file)
            self.__n_haplos = len(list(hap_assoc_reader)) - 1
        return self.__n_haplos

    def __load_snp_info(self, snp_info_file):
        # list all snps that are in the assoc.hap file
        hap_assoc_reader = HapAssocReader(file_name=self.__hap_assoc_file)
        hap_assoc_snps = {}
        hap_assoc_snps["SNP"] = 1
        for hap_assoc_rec in hap_assoc_reader:
            snp_list = hap_assoc_rec.snps.split("|")
            for snp in snp_list:
                hap_assoc_snps[snp] = 1
        # then only pick snp in the info file that are in the list above
        self.__snp_info_dict = {}
        snp_info_reader = SnpInfoReader(file_name=snp_info_file)
        snp_info_count = 0
        for rec in snp_info_reader:
            if rec.snp not in hap_assoc_snps:
                continue
            self.__snp_info_dict[rec.snp] = rec
            rec.snp_idx = snp_info_count
            snp_info_count += 1
        self.__n_snps = snp_info_count

    def __add_sheet(self, wb, sheet_name):
        ws = wb.add_worksheet(sheet_name)
        ws.set_default_row(12)
        ws.set_row(1, 35, None, {})
        ws.set_row(2, 35, None, {})
        ws.set_row(3, 35, None, {})
        ws.set_row(4, 45, None, {})
        return ws

    def __add_snp_to_ws(self, ws, row_idx, snp_info):
        cell_fmt = self.snp_fmts[NO_COLOR]
        ws.write(row_idx, COL_SNP_IDX, snp_info.snp, cell_format=cell_fmt)
        ws.write(row_idx, COL_F_MISS_A_IDX, snp_info.f_miss_a, cell_format=cell_fmt)
        ws.write(row_idx, COL_F_MISS_U_IDX, snp_info.f_miss_u, cell_format=cell_fmt)
        ws.write(row_idx, COL_CHROM_IDX, snp_info.chrom, cell_format=cell_fmt)
        ws.write(row_idx, COL_POS_IDX, snp_info.pos, cell_format=cell_fmt)

    def __add_snps_to_ws(self,
                         ws,
                         ):
        for snp in self.snp_info_dict:
            snp_info = self.snp_info_dict[snp]
            row_idx = HAPLO_INFO_SIZE + snp_info.snp_idx
            self.__add_snp_to_ws(ws, row_idx, snp_info)
        ws.set_column(COL_SNP_IDX, COL_SNP_IDX, 7)
        ws.set_column(COL_F_MISS_A_IDX, COL_F_MISS_U_IDX, 6)
        ws.set_column(COL_CHROM_IDX, COL_CHROM_IDX, 2)
        ws.set_column(COL_POS_IDX, COL_POS_IDX, 9.5)

    def __add_haplotype_stat_to_ws(self,
                                   ws,
                                   col,
                                   haplo_info,
                                   color,
                                   ):
        cell_fmt = self.haplo_fmts[color]
        ws.write(0, col, col+1)
        ws.write(1, col, haplo_info.f_a, cell_format=cell_fmt)
        ws.write(2, col, haplo_info.f_u, cell_format=cell_fmt)
        ws.write(3, col, haplo_info.ors, cell_format=cell_fmt)
        ws.write(4, col, haplo_info.pvalue, cell_format=cell_fmt)

    def __color_column(self, ws, col, color):
        cell_fmt = self.bp_fmts[color]
        for row_idx in xrange(HAPLO_INFO_SIZE, self.n_snps+HAPLO_INFO_SIZE):
            ws.write(row_idx, col, "", cell_format=cell_fmt)

    def __add_haplotypes_to_ws(self, ws):
        hap_assoc_reader = HapAssocReader(file_name=self.__hap_assoc_file)
        haplo_idx = 0
        for hap_assoc_rec in hap_assoc_reader:
            haplo_col = haplo_idx + SNP_INFO_SIZE - 1
            if hap_assoc_rec.pvalue < self.rpt_params.cutoff_pvalue:
                color = "YELLOW"
                self.__color_column(ws, haplo_col, color)
            else:
                color = NO_COLOR
            self.__add_haplotype_stat_to_ws(ws, haplo_col, hap_assoc_rec, color)
            # Map haplo to the corresponding markers
            bps_list = hap_assoc_rec.haplotype
            snps_list = hap_assoc_rec.snps.split("|")
            if bps_list == "HAPLOTYPE":
                haplo_idx += 1
                continue
            bp_fmt = self.bp_fmts[color]
            for bp_idx in xrange(len(bps_list)):
                bp = bps_list[bp_idx]
                snp = snps_list[bp_idx]
                if snp in self.snp_info_dict:
                    row_idx = self.snp_info_dict[snp].snp_idx + HAPLO_INFO_SIZE
                    ws.write(row_idx, haplo_col, bp, cell_format=bp_fmt)
            haplo_idx += 1
        ws.set_column(SNP_INFO_SIZE, self.n_haplos+SNP_INFO_SIZE-1, 1.6)

    def __set_report_format(self, wb):
        self.__plain_fmts = wb.add_colors_format({})
        self.__bp_fmts = wb.add_colors_format({'align': 'center'})
        self.__haplo_fmts = wb.add_colors_format({'rotation': 90})

    def __add_overview_haplos_sheet(self,
                                    wb,
                                    sheet_name,
                                    ):
        ws = self.__add_sheet(wb, sheet_name)
        self.__add_snps_to_ws(ws)
        self.__add_haplotypes_to_ws(ws)

    def gen_report(self, out_file):
        wb = Workbook(filename=out_file)
        self.__set_report_format(wb)
        self.__add_overview_haplos_sheet(wb, "overview")
        wb.close()
        self.info(" >>>> report file: " + out_file)
        return out_file
