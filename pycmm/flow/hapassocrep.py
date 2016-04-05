# This Python file uses the following encoding: utf-8
import operator
from os.path import join as join_path
from pycmm.template import pyCMMBase
from pycmm.cmmlib.xlslib import CMMWorkbook as Workbook
from pycmm.cmmlib.xlslib import NO_COLOR
from pycmm.flow import CMMPipeline
from pycmm.flow.plink import CMMPipeline
from pycmm.flow.plink import JOBS_SETUP_REPORT_PARAMS_SECTION
from pycmm.flow.plink import JOBS_SETUP_CUTOFF_PVALUE_KEY
from pycmm.flow.plink import JOBS_SETUP_CUTOFF_ORS_KEY
from pycmm.flow.plink import JOBS_SETUP_FAMLIY_HAPLOTYPE_PREFIX_KEY
from pycmm.flow.plink import DFLT_CUTOFF_PVALUE
from pycmm.flow.plink import DFLT_CUTOFF_ORS
from pycmm.cmmlib.plinklib import SnpInfoReader
from pycmm.cmmlib.plinklib import HapAssocReader
from pycmm.cmmlib.plinklib import FamReader
from pycmm.cmmlib.plinklib import TPedReader

HAPLO_INFO_SIZE = 4
SNP_INFO_SIZE = 5

COL_SNP_IDX = 0
COL_F_MISS_A_IDX = 1
COL_F_MISS_U_IDX = 2
COL_CHROM_IDX = 3
COL_POS_IDX = 4

GT0_COLOR = "XLS_GREEN"
GT1_COLOR = "XLS_PINK"
SIG_ORS_COLOR = "ICEBLUE"
BEST_ORS_COLOR = "ROYAL_BLUE"
SIG_PVALUE_COLOR = "YELLOW"
BEST_PVALUE_COLOR = "GOLD"

HAPLO_COUNT_LOG_INTERVAL = 1000


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

    @property
    def cutoff_ors(self):
        if JOBS_SETUP_CUTOFF_ORS_KEY in self.__params:
            return self.__params[JOBS_SETUP_CUTOFF_ORS_KEY]
        return DFLT_CUTOFF_ORS

    @property
    def fam_hap_prefix(self):
        if JOBS_SETUP_FAMLIY_HAPLOTYPE_PREFIX_KEY in self.__params:
            return self.__params[JOBS_SETUP_FAMLIY_HAPLOTYPE_PREFIX_KEY]
        return None

class HapAssocRepPipeline(CMMPipeline):
    """ To control production of haplotype association study reports """

    def __init__(self,
                 hap_assoc_file,
                 snp_info_file,
                 **kwargs
                 ):
        self.__init_properties()
        self.__hap_assoc_file = hap_assoc_file
        self.__snp_info_file = snp_info_file
        super(HapAssocRepPipeline, self).__init__(**kwargs)
        self.__load_families_haplotypes()

    def __init_properties(self):
        self.__rpt_params = None
        self.__fams_hap_idx = None
        self.__hap_assoc_snps = None
        self.__snps_info = None
        self.__n_snps = None
        self.__sorted_ors_haplotypes = None
        self.__sorted_pvalue_haplotypes = None

    @property
    def rpt_params(self):
        if self.__rpt_params is None:
            self.__rpt_params = RptParams(self._jobs_info[JOBS_SETUP_REPORT_PARAMS_SECTION])
        return self.__rpt_params

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
    def hap_stat_fmts(self):
        return self.__hap_stat_fmts

    @property
    def n_snps(self):
        if self.__n_snps is None:
            self.snps_info
        return self.__n_snps

    @property
    def hap_assoc_snps(self):
        if self.__hap_assoc_snps is None:
            # list all snps that are in the assoc.hap file
            hap_assoc_reader = HapAssocReader(file_name=self.__hap_assoc_file)
            self.__hap_assoc_snps = {}
            self.__hap_assoc_snps["SNP"] = 1
            for hap_assoc_rec in hap_assoc_reader:
                snp_list = hap_assoc_rec.snps.split("|")
                for snp in snp_list:
                    self.__hap_assoc_snps[snp] = 1
        return self.__hap_assoc_snps

    @property
    def snps_info(self):
        if self.__snps_info is None:
            # load information of snps from snp.info file but only pick snps
            # that are in the haplotype association study result
            self.__snps_info = {}
            snp_info_reader = SnpInfoReader(file_name=self.__snp_info_file)
            snp_info_count = 0
            for snp_info_rec in snp_info_reader:
                if snp_info_rec.snp not in self.hap_assoc_snps:
                    continue
                self.__snps_info[snp_info_rec.snp] = snp_info_rec
                snp_info_rec.snp_idx = snp_info_count
                snp_info_count += 1
            self.__n_snps = snp_info_count - 1
        return self.__snps_info

    @property
    def sorted_ors_haplotypes(self):
        if self.__sorted_ors_haplotypes is None:
            self.__sorted_ors_haplotypes = self.__get_sorted_haplotypes(attr="ors",
                                                                        reverse=True,
                                                                        )
        return self.__sorted_ors_haplotypes

    @property
    def sorted_pvalue_haplotypes(self):
        if self.__sorted_pvalue_haplotypes is None:
            self.__sorted_pvalue_haplotypes = self.__get_sorted_haplotypes(attr="pvalue")
        return self.__sorted_pvalue_haplotypes

    def __compare_haplotypes(self, hap_assoc_rec, fam_gts):
        hap_assoc_bps_list = hap_assoc_rec.haplotype
        hap_assoc_snps_list = hap_assoc_rec.snps.split("|")
        matched_allele = -1
        for allele_idx in xrange(2):
            matched = True
            compared = False
            for bp_idx in xrange(len(hap_assoc_bps_list)):
                hap_assoc_bp = hap_assoc_bps_list[bp_idx]
                snp = hap_assoc_snps_list[bp_idx]
                # allow matching with if no genotype data
                if snp not in fam_gts:
                    continue
                # allow matching with empty (uncertain) family genotype
                if fam_gts[snp][allele_idx] == "0":
                    continue
                # compared
                if hap_assoc_bp != fam_gts[snp][allele_idx]:
                    matched = False
                    break
                compared = True
            if matched and compared:
                matched_allele = allele_idx
                break
        return matched_allele

    def __get_best_hap_idx(self, fam_gts, hap_assoc_recs):
        for hap_assoc_rec in hap_assoc_recs:
            if len(hap_assoc_rec.haplotype) == 1:
                continue
            matched_allele = self.__compare_haplotypes(hap_assoc_rec, fam_gts)
            if matched_allele != -1:
                return hap_assoc_rec.col_idx
        return None

    def __get_best_ors_hap_idx(self, fam_gts):
        return self.__get_best_hap_idx(fam_gts, self.sorted_ors_haplotypes)

    def __get_best_pvalue_hap_idx(self, fam_gts):
        return self.__get_best_hap_idx(fam_gts, self.sorted_pvalue_haplotypes)

    def __get_sorted_haplotypes(self, attr, reverse=False):
        hap_assoc_reader = HapAssocReader(file_name=self.__hap_assoc_file)
        hap_assoc_reader.next()
        hap_assoc_recs = list(hap_assoc_reader)
        for idx in xrange(len(hap_assoc_recs)):
            hap_assoc_recs[idx].col_idx = idx
        hap_assoc_recs.sort(key=operator.attrgetter(attr),
                            reverse=reverse)
        return hap_assoc_recs

    def __load_families_haplotypes(self):
        # load genotyping information from tped and tfam files and add them 
        # to FamilyInfo object
        if ((self.__fams_hap_idx is None) and
            (self.rpt_params.fam_hap_prefix is not None) and
            (self.families_info is not None)
            ):
            self.__fam_hap_idxs = {}
            tfam_file_name = self.rpt_params.fam_hap_prefix + ".tfam"
            tped_file_name = self.rpt_params.fam_hap_prefix + ".tped"
            # read the file using TFamReader because they have
            # the same structures
            tfam_reader = FamReader(file_name=tfam_file_name)
            fam_hap_idx = 0
            # get row idx(0) from tfam file
            for fam_rec in tfam_reader:
                fam_hap_key = fam_rec.fam_id + "_" + fam_rec.indv_id
                self.__fam_hap_idxs[fam_hap_key] = fam_hap_idx
                fam_hap_idx += 1
            # load complete genotyping information to each sample(family)
            # specified in 'sample_info' parameters
            for fam_id in self.families_info:
                fam_info = self.families_info[fam_id]
                for member_info in fam_info.members_info:
                    fam_hap_key = fam_id + "_" + member_info.sample_id
                    if fam_hap_key not in self.__fam_hap_idxs:
                        msg = "Member '" + member_info.sample_id + "' of"
                        msg += " family '" + fam_id + "' is not found in"
                        msg += " family haplotype information"
                        msg += " (" + tfam_file_name + ")"
                        self.warning(msg)
                        continue
                    fam_hap_idx = self.__fam_hap_idxs[fam_hap_key]
                    tped_reader = TPedReader(file_name=tped_file_name)
                    gts = {}
                    for tped_rec in tped_reader:
                        # Only load the genoytpes that are in the haplotype
                        # association study result
                        if tped_rec.snp not in self.hap_assoc_snps:
                            continue
                        gts[tped_rec.snp] = tped_rec.gts[fam_hap_idx].split()
                    member_info.gts = gts

    def __add_sheet(self, sheet_name):
        ws = self.__wb.add_worksheet(sheet_name)
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
        for snp in self.snps_info:
            snp_info = self.snps_info[snp]
            row_idx = HAPLO_INFO_SIZE + snp_info.snp_idx
            self.__add_snp_to_ws(ws, row_idx, snp_info)
        ws.set_column(COL_SNP_IDX, COL_SNP_IDX, 7)
        ws.set_column(COL_F_MISS_A_IDX, COL_F_MISS_U_IDX, 6)
        ws.set_column(COL_CHROM_IDX, COL_CHROM_IDX, 2)
        ws.set_column(COL_POS_IDX, COL_POS_IDX, 9.5)

    def __add_haplotype_stat_to_ws(self,
                                   ws,
                                   col,
                                   hap_info,
                                   pvalue_color,
                                   ors_color,
                                   ):
        if pvalue_color == NO_COLOR:
            pvalue_color = ors_color
        pvalue_fmt = self.hap_stat_fmts[pvalue_color]
        ors_fmt = self.hap_stat_fmts[ors_color]
        ws.write(0, col, col+1)
        ws.write(1, col, hap_info.f_a, cell_format=ors_fmt)
        ws.write(2, col, hap_info.f_u, cell_format=ors_fmt)
        ws.write(3, col, hap_info.ors, cell_format=ors_fmt)
        ws.write(4, col, hap_info.pvalue, cell_format=pvalue_fmt)

    def __color_column(self, ws, col, color):
        cell_fmt = self.bp_fmts[color]
        for row_idx in xrange(HAPLO_INFO_SIZE+1, self.n_snps+HAPLO_INFO_SIZE+1):
            ws.write(row_idx, col, "", cell_format=cell_fmt)

    def __add_haplotypes_to_ws(self,
                               ws,
                               hap_start_col=SNP_INFO_SIZE,
                               fam_gts=None,
                               best_ors_hap_idx=None,
                               best_pvalue_hap_idx=None,
                               ):
        hap_assoc_reader = HapAssocReader(file_name=self.__hap_assoc_file)
        # writer header
        header_rec = hap_assoc_reader.next()
        self.__add_haplotype_stat_to_ws(ws,
                                        col=SNP_INFO_SIZE-1,
                                        hap_info=header_rec,
                                        pvalue_color=NO_COLOR,
                                        ors_color=NO_COLOR,
                                        )
        # write content
        hap_col = hap_start_col
        hap_count = 0
        for hap_assoc_rec in hap_assoc_reader:
            # define color of the column using cut-off p-value and 
            # cut-off odds ratio
            if ((best_pvalue_hap_idx is not None) and
                (hap_count == best_pvalue_hap_idx)
                ):
                pvalue_color = BEST_PVALUE_COLOR
            elif hap_assoc_rec.pvalue < self.rpt_params.cutoff_pvalue:
                pvalue_color = SIG_PVALUE_COLOR
            else:
                pvalue_color = NO_COLOR
            if ((best_ors_hap_idx is not None) and
                (hap_count == best_ors_hap_idx)
                ):
                ors_color = BEST_ORS_COLOR
            elif hap_assoc_rec.ors > self.rpt_params.cutoff_ors:
                ors_color = SIG_ORS_COLOR
            else:
                ors_color = pvalue_color
            self.__add_haplotype_stat_to_ws(ws,
                                            col=hap_col,
                                            hap_info=hap_assoc_rec,
                                            pvalue_color=pvalue_color,
                                            ors_color=ors_color,
                                            )
            if fam_gts is not None:
                matched_allele = self.__compare_haplotypes(hap_assoc_rec,
                                                           fam_gts)
                if matched_allele == -1:
                    bp_color = ors_color
                elif matched_allele == 0:
                    bp_color = GT0_COLOR
                elif matched_allele == 1:
                    bp_color = GT1_COLOR
            else:
                bp_color = ors_color
            if bp_color != NO_COLOR:
                self.__color_column(ws, hap_col, bp_color)
            # Map haplo to the corresponding markers
            hap_assoc_bps_list = hap_assoc_rec.haplotype
            hap_assoc_snps_list = hap_assoc_rec.snps.split("|")
            bp_fmt = self.bp_fmts[bp_color]
            for bp_idx in xrange(len(hap_assoc_bps_list)):
                bp = hap_assoc_bps_list[bp_idx]
                snp = hap_assoc_snps_list[bp_idx]
                if snp in self.snps_info:
                    row_idx = self.snps_info[snp].snp_idx + HAPLO_INFO_SIZE
                    ws.write(row_idx, hap_col, bp, cell_format=bp_fmt)
            hap_col += 1
            hap_count += 1
            if hap_count % HAPLO_COUNT_LOG_INTERVAL == 0:
                log_msg = str(hap_count)
                log_msg += " haplotypes were written to the sheet"
                self.info(log_msg)
        log_msg = "Finish .. "
        log_msg += " total of " + str(hap_count)
        log_msg += " haplotypes were written to the sheet"
        self.info(log_msg)
        ws.set_column(hap_start_col, hap_start_col+hap_count-1, 1.6)

    def __set_report_format(self):
        self.__plain_fmts = self.__wb.add_colors_format({})
        self.__bp_fmts = self.__wb.add_colors_format({'align': 'center'})
        self.__hap_stat_fmts = self.__wb.add_colors_format({'rotation': 90})

    def __add_fam_hap_to_ws(self,
                            ws,
                            col_name,
                            gts,
                            allele_idx,
                            col_idx,
                            color=NO_COLOR):
        gt_header_row = HAPLO_INFO_SIZE
        gt_header_fmt = self.hap_stat_fmts[color]
        ws.write(gt_header_row,
                 col_idx,
                 col_name,
                 cell_format=gt_header_fmt)
        bp_fmt = self.plain_fmts[color]
        for snp in self.snps_info:
            if snp == "SNP":
                continue
            row_idx = self.snps_info[snp].snp_idx + HAPLO_INFO_SIZE
            if snp in gts:
                bp = gts[snp][allele_idx]
                if bp == "0":
                    ws.write(row_idx, col_idx, "", cell_format=bp_fmt)
                else:
                    ws.write(row_idx, col_idx, bp, cell_format=bp_fmt)
            else:
                ws.write(row_idx, col_idx, "X", cell_format=bp_fmt)

    def __add_fam_haps_to_ws(self,
                             ws,
                             fam_id,
                             gts,
                             gt_start_col=SNP_INFO_SIZE):
        self.__add_fam_hap_to_ws(ws,
                                 col_name=fam_id+" (shared)",
                                 gts=gts,
                                 allele_idx=0,
                                 col_idx=gt_start_col,
                                 color=GT0_COLOR,
                                 )
        self.__add_fam_hap_to_ws(ws,
                                 col_name=fam_id+" (unshared)",
                                 gts=gts,
                                 allele_idx=1,
                                 col_idx=gt_start_col+1,
                                 color=GT1_COLOR,
                                 )
        ws.set_column(gt_start_col, gt_start_col+1, 1.6)

    def __add_fam_haps_sheet(self,
                             fam_id,
                             gts,
                             best_ors_hap_idx=None,
                             best_pvalue_hap_idx=None,
                             ):
        sheet_name = "Family " + fam_id
        ws = self.__add_sheet(sheet_name)
        self.info(" >> add '" + sheet_name + "' sheet")
        self.__add_snps_to_ws(ws)
        self.__add_fam_haps_to_ws(ws,
                                  fam_id=fam_id,
                                  gts=gts,
                                  )
        self.__add_haplotypes_to_ws(ws,
                                    hap_start_col=SNP_INFO_SIZE+2,
                                    fam_gts=gts,
                                    best_ors_hap_idx=best_ors_hap_idx,
                                    best_pvalue_hap_idx=best_pvalue_hap_idx,
                                    )
        ws.freeze_panes(HAPLO_INFO_SIZE+1, SNP_INFO_SIZE+2)

    def __add_fams_haps_sheet(self):
        for fam_id in self.families_info:
            fam_info = self.families_info[fam_id]
            for member_info in fam_info.members_info:
                if member_info.gts is None:
                    continue
                best_ors_hap_idx = self.__get_best_ors_hap_idx(member_info.gts)
                best_pvalue_hap_idx = self.__get_best_pvalue_hap_idx(member_info.gts)
                self.__add_fam_haps_sheet(fam_id=fam_id,
                                          gts=member_info.gts,
                                          best_ors_hap_idx=best_ors_hap_idx,
                                          best_pvalue_hap_idx=best_pvalue_hap_idx,
                                          )

    def __add_legend_sheet(self):
        ws = self.__wb.add_worksheet("legend")
        cell_txt = "matched with family 'shared' haplotype"
        ws.write(0, 0, cell_txt, cell_format = self.plain_fmts[GT0_COLOR])
        cell_txt = "matched with family 'unshared' haplotype"
        ws.write(1, 0, cell_txt, cell_format = self.plain_fmts[GT1_COLOR])
        cell_txt = "haplotye with significant odds ratio (more than "
        cell_txt += str(self.rpt_params.cutoff_ors) + ")"
        ws.write(2, 0, cell_txt, cell_format = self.plain_fmts[SIG_ORS_COLOR])
        cell_txt = "haplotye that match with the family with the best odds ratio"
        ws.write(3, 0, cell_txt, cell_format = self.plain_fmts[BEST_ORS_COLOR])
        cell_txt = "haplotye with significant p-value (less than "
        cell_txt += str(self.rpt_params.cutoff_pvalue) + ")"
        ws.write(4, 0, cell_txt, cell_format = self.plain_fmts[SIG_PVALUE_COLOR])
        cell_txt = "haplotye that match with the family with the best odds ratio"
        ws.write(5, 0, cell_txt, cell_format = self.plain_fmts[BEST_PVALUE_COLOR])
        ws.set_column(0, 0, 40)

    def gen_report(self, out_file):
        self.__wb = Workbook(filename=out_file)
        self.info(" >>>> generating haplotype association study report")
        self.info(" >>>> report file: " + out_file)
        self.__set_report_format()
        self.__add_fams_haps_sheet()
        self.__add_legend_sheet()
        self.__wb.close()
        return out_file
