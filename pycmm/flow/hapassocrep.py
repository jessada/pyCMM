# This Python file uses the following encoding: utf-8
import xlsxwriter
from os.path import join as join_path
from pycmm.flow import CMMPipeline
from pycmm.cmmlib.plinklib import SnpInfoReader
from pycmm.cmmlib.plinklib import HapAssocReader

HAPLO_INFO_SIZE = 5
SNP_INFO_SIZE = 5

COL_SNP_IDX = 0
COL_F_MISS_A_IDX = 1
COL_F_MISS_U_IDX = 2
COL_CHROM_IDX = 3
COL_POS_IDX = 4


class HapAssocRepPipeline(CMMPipeline):
    """ To control production of haplotype association study reports """

    def __init__(self,
                 hap_assoc_file,
                 snp_info_file,
                 **kwargs
                 ):
        self.__hap_assoc_file = hap_assoc_file
        self.__load_snp_info(snp_info_file)
        super(HapAssocRepPipeline, self).__init__(**kwargs)

    @property
    def hap_assoc_report_file(self):
        return join_path(self.rpts_out_dir,
                         self.project_name+".xlsx")

    @property
    def snp_info_dict(self):
        return self.__snp_info_dict

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
#        rec = snp_info_reader.next()
#        self.__snp_info_dict[rec.snp] = rec
#        rec.snp_idx = snp_info_count
#        snp_info_count += 1
        for rec in snp_info_reader:
            if rec.snp not in hap_assoc_snps:
                continue
            self.__snp_info_dict[rec.snp] = rec
            rec.snp_idx = snp_info_count
            snp_info_count += 1

    def __add_sheet(self, wb, sheet_name):
        ws = wb.add_worksheet(sheet_name)
        ws.set_default_row(12)
        return ws

    def __add_snp_to_ws(self, ws, row_idx, snp_info):
        ws.write(row_idx, COL_SNP_IDX, snp_info.snp)
        ws.write(row_idx, COL_F_MISS_A_IDX, snp_info.f_miss_a)
        ws.write(row_idx, COL_F_MISS_U_IDX, snp_info.f_miss_u)
        ws.write(row_idx, COL_CHROM_IDX, snp_info.chrom)
        ws.write(row_idx, COL_POS_IDX, snp_info.pos)

    def __add_snps_to_ws(self,
                         ws,
                         ):
        for snp in self.snp_info_dict:
            snp_info = self.snp_info_dict[snp]
            row_idx = HAPLO_INFO_SIZE + snp_info.snp_idx
            self.__add_snp_to_ws(ws, row_idx, snp_info)
        ws.set_column(COL_SNP_IDX, COL_SNP_IDX, 7)
        ws.set_column(COL_F_MISS_A_IDX, COL_F_MISS_U_IDX, 6)
        ws.set_column(COL_CHROM_IDX, COL_CHROM_IDX, 0.8)
        ws.set_column(COL_POS_IDX, COL_POS_IDX, 10)

    def __add_haplotype_stat_to_ws(self,
                                   ws,
                                   col,
                                   haplo_info,
                                   ):
        ws.write(0, col, col+1)
        ws.write(1, col, haplo_info.f_a)
        ws.write(2, col, haplo_info.f_u)
        ws.write(3, col, haplo_info.chisq)
        ws.write(4, col, haplo_info.ors)
        ws.write(5, col, haplo_info.pvalue)

    def __add_haplotypes_to_ws(self,
                                ws,
                                ):
        hap_assoc_reader = HapAssocReader(file_name=self.__hap_assoc_file)
        haplo_idx = 0
        for hap_assoc_rec in hap_assoc_reader:
            haplo_col = haplo_idx + SNP_INFO_SIZE - 1
            self.__add_haplotype_stat_to_ws(ws, haplo_col, hap_assoc_rec)
            # Map haplo to the corresponding markers
            bps_list = hap_assoc_rec.haplotype
            snps_list = hap_assoc_rec.snps.split("|")
            if bps_list == "HAPLOTYPE":
                haplo_idx += 1
                continue
            for bp_idx in xrange(len(bps_list)):
                bp = bps_list[bp_idx]
                snp = snps_list[bp_idx]
                if snp in self.snp_info_dict:
                    row_idx = self.snp_info_dict[snp].snp_idx + HAPLO_INFO_SIZE
                    ws.write(row_idx, haplo_col, bp)
            haplo_idx += 1

    def __add_overview_haplos_sheet(self,
                                    wb,
                                    sheet_name,
                                    ):
        ws = self.__add_sheet(wb, sheet_name)
        self.__add_snps_to_ws(ws)
        self.__add_haplotypes_to_ws(ws)

    def gen_report(self):
        wb = xlsxwriter.Workbook(self.hap_assoc_report_file)
        self.__add_overview_haplos_sheet(wb, "overview")
        wb.close()
        self.info(" >>>> report file: " + self.hap_assoc_report_file)
        return self.hap_assoc_report_file
