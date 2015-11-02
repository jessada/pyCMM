import os

"""

This module is all about constant value 

"""

#FULL_SYSTEM_TEST = True
FULL_SYSTEM_TEST = False
SLURM_TEST = False

FAST_PROJECT_CODE = 'b2011097'
SLOW_PROJECT_CODE = 'b2011158'

ENV_TEST_DIR = "PYTHON_TEST_DIR"

## > > > > > > > > > > > > > GATK Best practice (DNA seq) < < < < < < < < < <
DNASEQ_PIPELINE_DESCRIPTION = "A flow to control a pipeline to process DNA sequencing data"
DNASEQ_PIPELINE_DFLT_LOG_FILE = "dnaseq_pipeline"
DFLT_GATKBP_ALLOC_TIME = "4-00:00:00"

DNASEQ_CREATE_JOB_SETUP_FILE_DESCRIPTION = "An applicationt to generate job setup file for DNASEQ_PIPELINE"

## > > > > > > > > > > > > > CMM DB < < < < < < < < < <
CMMDB_PIPELINE_DESCRIPTION = "A flow to control a pipeline to process CMM database"
CMMDB_MUTSTAT_DESCRIPTION = "A flow to control a pipeline to process mutation statistics database"
CMMDB_TABLEANNOVAR_DESCRIPTION = "A flow to control a pipeline to run tableannovar"
CMMDB_MUTATIONREPORTS_DESCRIPTION = "A flow to control a pipeline to run generate mutation reports"
CMMDB_PIPELINE_DFLT_LOG_FILE = "cmmdb_pipeline"
DFLT_CMMDB_ALLOC_TIME = "1-00:00:00"

CMMDB_CREATE_JOB_SETUP_FILE_DESCRIPTION = "An applicationt to generate job setup file to process CMM database"

## > > > > > > > > > > > > > Mutation report < < < < < < < < < <
MUTREP_FAMILY_REPORT_BIN = 'pyCMM-mutrep-family-report'
MUTREP_SUMMARY_REPORT_BIN = 'pyCMM-mutrep-summary-report'

MUTREP_PIPELINE_DESCRIPTION = "A flow to control a pipeline to process mutations report"
MUTREP_PIPELINE_DFLT_LOG_FILE = "mutrep_pipeline"
DFLT_MUTREP_ALLOC_TIME = "1-00:00:00"

MUTREP_FAMILY_REPORT_DESCRIPTION = "An appliation to generate mutation report for a given family at given regions"
MUTREP_SUMMARY_REPORT_DESCRIPTION = "An appliation to generate summary report at given regions"
MUTREP_CREATE_JOB_SETUP_FILE_DESCRIPTION = "An application to generate job setup file to process mutations report"

DFLT_MAF_VAR = "1000g2014oct_all"
DFLT_MUTREP_FREQ_RATIOS = DFLT_MAF_VAR + ":0.2"
FUNC_REFGENE_VAR = "Func.refGene"

PREDICTION_LIST = []
LJB_SIFT_PREDICTION_COL_NAME = "SIFT_pred"
PREDICTION_LIST.append(LJB_SIFT_PREDICTION_COL_NAME)
LJB_POLYPHEN2_HDIV_PREDICTION_COL_NAME = "Polyphen2_HDIV_pred"
PREDICTION_LIST.append(LJB_POLYPHEN2_HDIV_PREDICTION_COL_NAME)
LJB_POLYPHEN2_HVAR_PREDICTION_COL_NAME = "Polyphen2_HVAR_pred"
PREDICTION_LIST.append(LJB_POLYPHEN2_HVAR_PREDICTION_COL_NAME)
LJB_LRT_PREDICTION_COL_NAME = "LRT_pred"
PREDICTION_LIST.append(LJB_LRT_PREDICTION_COL_NAME)
LJB_MUTATIONTASTER_PREDICTION_COL_NAME = "MutationTaster_pred"
PREDICTION_LIST.append(LJB_MUTATIONTASTER_PREDICTION_COL_NAME)
LJB_MUTATIONASSESSOR_PREDICTION_COL_NAME = "MutationAssessor_pred"
PREDICTION_LIST.append(LJB_MUTATIONASSESSOR_PREDICTION_COL_NAME)
LJB_FATHMM_PREDICTION_COL_NAME = "FATHMM_pred"
PREDICTION_LIST.append(LJB_FATHMM_PREDICTION_COL_NAME)
LJB_RADIALSVM_PREDICTION_COL_NAME = "RadialSVM_pred"
PREDICTION_LIST.append(LJB_RADIALSVM_PREDICTION_COL_NAME)
LJB_LR_PREDICTION_COL_NAME = "LR_pred"
PREDICTION_LIST.append(LJB_LR_PREDICTION_COL_NAME)

# This list will show all possible annotation columns (based on the dbs
# annotated by table_annovar.pl), in the INFO fields.
#
# Note. I use variable names instead of hard-code because I want to emphasize 
# that these columns are used somewhere else.
MT_ANNO_COLS = []
MT_ANNO_COLS.append("Gene_symbol")
MT_ANNO_COLS.append("OXPHOS_Complex")
MT_ANNO_COLS.append("Ensembl_Gene_ID")
MT_ANNO_COLS.append("Ensembl_Protein_ID")
MT_ANNO_COLS.append("Uniprot_Name")
MT_ANNO_COLS.append("Uniprot_ID")
MT_ANNO_COLS.append("NCBI_Gene_ID")
MT_ANNO_COLS.append("NCBI_Protein_ID")
MT_ANNO_COLS.append("Gene_pos")
MT_ANNO_COLS.append("AA_pos")
MT_ANNO_COLS.append("AA_sub")
MT_ANNO_COLS.append("Codon_sub")
MT_ANNO_COLS.append("dbSNP_ID")
MT_ANNO_COLS.append("PhyloP_46V")
MT_ANNO_COLS.append("PhastCons_46V")
MT_ANNO_COLS.append("PhyloP_100V")
MT_ANNO_COLS.append("PhastCons_100V")
MT_ANNO_COLS.append("SiteVar")
MT_ANNO_COLS.append("PolyPhen2_prediction")
MT_ANNO_COLS.append("PolyPhen2_score")
MT_ANNO_COLS.append("SIFT_prediction")
MT_ANNO_COLS.append("SIFT_score")
MT_ANNO_COLS.append("FatHmm_prediction")
MT_ANNO_COLS.append("FatHmm_score")
MT_ANNO_COLS.append("PROVEAN_prediction")
MT_ANNO_COLS.append("PROVEAN_score")
MT_ANNO_COLS.append("MutAss_prediction")
MT_ANNO_COLS.append("MutAss_score")
MT_ANNO_COLS.append("EFIN_Swiss_Prot_Score")
MT_ANNO_COLS.append("EFIN_Swiss_Prot_Prediction")
MT_ANNO_COLS.append("EFIN_HumDiv_Score")
MT_ANNO_COLS.append("EFIN_HumDiv_Prediction")
MT_ANNO_COLS.append("CADD_score")
MT_ANNO_COLS.append("CADD_Phred_score")
MT_ANNO_COLS.append("CADD_prediction")
MT_ANNO_COLS.append("Carol_prediction")
MT_ANNO_COLS.append("Carol_score")
MT_ANNO_COLS.append("Condel_score")
MT_ANNO_COLS.append("Condel_pred")
MT_ANNO_COLS.append("COVEC_WMV")
MT_ANNO_COLS.append("COVEC_WMV_prediction")
MT_ANNO_COLS.append("PolyPhen2_score_transf")
MT_ANNO_COLS.append("PolyPhen2_pred_transf")
MT_ANNO_COLS.append("SIFT_score_transf")
MT_ANNO_COLS.append("SIFT_pred_transf")
MT_ANNO_COLS.append("MutAss_score_transf")
MT_ANNO_COLS.append("MutAss_pred_transf")
MT_ANNO_COLS.append("Perc_coevo_Sites")
MT_ANNO_COLS.append("Mean_MI_score")
MT_ANNO_COLS.append("COSMIC_ID")
MT_ANNO_COLS.append("Tumor_site")
MT_ANNO_COLS.append("Examined_samples")
MT_ANNO_COLS.append("Mutation_frequency")
MT_ANNO_COLS.append("US")
MT_ANNO_COLS.append("Status")
MT_ANNO_COLS.append("Associated_disease")
MT_ANNO_COLS.append("Presence_in_TD")
MT_ANNO_COLS.append("Class_predicted")
MT_ANNO_COLS.append("Prob_N")
MT_ANNO_COLS.append("Prob_P")

DFLT_MUTREP_ANNO_COLS = []
DFLT_MUTREP_ANNO_COLS.append(FUNC_REFGENE_VAR)
DFLT_MUTREP_ANNO_COLS.append("ExonicFunc.refGene")
DFLT_MUTREP_ANNO_COLS.append("Gene.refGene")
DFLT_MUTREP_ANNO_COLS.append("GeneDetail.refGene")
DFLT_MUTREP_ANNO_COLS.append("snp138")
DFLT_MUTREP_ANNO_COLS.append("1000g2014oct_eur")
DFLT_MUTREP_ANNO_COLS.append(DFLT_MAF_VAR)
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_PF")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_GF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_PF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_GF")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_PF")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_GF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_PF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_GF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_PF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_GF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_PF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_GF")
DFLT_MUTREP_ANNO_COLS.append("249_SWEDES")
DFLT_MUTREP_ANNO_COLS.append("200_DANES")
DFLT_MUTREP_ANNO_COLS.append("ExAC_ALL")
DFLT_MUTREP_ANNO_COLS.append("ExAC_AFR")
DFLT_MUTREP_ANNO_COLS.append("ExAC_AMR")
DFLT_MUTREP_ANNO_COLS.append("ExAC_EAS")
DFLT_MUTREP_ANNO_COLS.append("ExAC_FIN")
DFLT_MUTREP_ANNO_COLS.append("ExAC_NFE")
DFLT_MUTREP_ANNO_COLS.append("ExAC_OTH")
DFLT_MUTREP_ANNO_COLS.append("ExAC_SAS")
DFLT_MUTREP_ANNO_COLS.append("AAChange.refGene")
DFLT_MUTREP_ANNO_COLS.append("cytoBand")
DFLT_MUTREP_ANNO_COLS.append("genomicSuperDups")
DFLT_MUTREP_ANNO_COLS.append("esp6500siv2_all")
DFLT_MUTREP_ANNO_COLS.append("cg69")
DFLT_MUTREP_ANNO_COLS.append("cosmic70")
DFLT_MUTREP_ANNO_COLS.append("nci60")
DFLT_MUTREP_ANNO_COLS.append("clinvar_20150330")
DFLT_MUTREP_ANNO_COLS += MT_ANNO_COLS
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_WT")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_HET")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_HOM")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_OTH")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_NA")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_GT")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_AF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_WT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_HET")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_HOM")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_OTH")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_NA")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_GT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_ALL_AF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_WT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_HET")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_HOM")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_OTH")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_NA")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_GT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_PF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_AF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CAFAM_GF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_WT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_HET")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_HOM")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_OTH")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_NA")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_GT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_PF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_AF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_COLON_GF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_WT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_HET")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_HOM")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_OTH")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_NA")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_GT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_PF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_AF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_CRC_GF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_WT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_HET")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_HOM")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_OTH")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_NA")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_GT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_PF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_AF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_RECTAL_GF")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_WT")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_HET")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_HOM")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_OTH")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_NA")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_GT")
DFLT_MUTREP_ANNO_COLS.append("PMS2_SPORADIC_FAM24_AF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_WT")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_HET")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_HOM")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_OTH")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_NA")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_GT")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR3_6_14_18_AF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_WT")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_HET")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_HOM")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_OTH")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_NA")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_GT")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR5_19_AF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_WT")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_HET")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_HOM")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_OTH")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_NA")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_GT")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_CHR9_AF")
DFLT_MUTREP_ANNO_COLS.append("SIFT_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_SIFT_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("Polyphen2_HDIV_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_POLYPHEN2_HDIV_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("Polyphen2_HVAR_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_POLYPHEN2_HVAR_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("LRT_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_LRT_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("MutationTaster_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_MUTATIONTASTER_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("MutationAssessor_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_MUTATIONASSESSOR_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("FATHMM_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_FATHMM_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("RadialSVM_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_RADIALSVM_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("LR_score")
DFLT_MUTREP_ANNO_COLS.append(LJB_LR_PREDICTION_COL_NAME)
DFLT_MUTREP_ANNO_COLS.append("VEST3_score")
DFLT_MUTREP_ANNO_COLS.append("CADD_raw")
DFLT_MUTREP_ANNO_COLS.append("CADD_phred")
DFLT_MUTREP_ANNO_COLS.append("GERP++_RS")
DFLT_MUTREP_ANNO_COLS.append("phyloP46way_placental")
DFLT_MUTREP_ANNO_COLS.append("phyloP100way_vertebrate")
DFLT_MUTREP_ANNO_COLS.append("SiPhy_29way_logOdds")

# This part will remove some columns that are not necessary to be shown in
# normal reports, for example, some STAT columns.
#DFLT_MUTREP_ANNO_COLS.remove("ILL_BR_WT")
#DFLT_MUTREP_ANNO_COLS.remove("ILL_BR_HET")
#DFLT_MUTREP_ANNO_COLS.remove("ILL_BR_HOM")
#DFLT_MUTREP_ANNO_COLS.remove("ILL_BR_OTH")
#DFLT_MUTREP_ANNO_COLS.remove("ILL_BR_NA")
#DFLT_MUTREP_ANNO_COLS.remove("ILL_BR_GT")
#DFLT_MUTREP_ANNO_COLS.remove("ILL_BR_AF")
#DFLT_MUTREP_ANNO_COLS.remove("AXEQ_BR_WT")
#DFLT_MUTREP_ANNO_COLS.remove("AXEQ_BR_HET")
#DFLT_MUTREP_ANNO_COLS.remove("AXEQ_BR_HOM")
#DFLT_MUTREP_ANNO_COLS.remove("AXEQ_BR_OTH")
#DFLT_MUTREP_ANNO_COLS.remove("AXEQ_BR_NA")
#DFLT_MUTREP_ANNO_COLS.remove("AXEQ_BR_GT")
#DFLT_MUTREP_ANNO_COLS.remove("AXEQ_BR_AF")
#DFLT_MUTREP_ANNO_COLS.remove("101CRC_WT")
#DFLT_MUTREP_ANNO_COLS.remove("101CRC_HET")
#DFLT_MUTREP_ANNO_COLS.remove("101CRC_HOM")
#DFLT_MUTREP_ANNO_COLS.remove("101CRC_OTH")
#DFLT_MUTREP_ANNO_COLS.remove("101CRC_NA")
#DFLT_MUTREP_ANNO_COLS.remove("101CRC_GT")
#DFLT_MUTREP_ANNO_COLS.remove("101CRC_AF")


## > > > > > > > > > > > > > ANNOVAR < < < < < < < < < <
DFLT_ANNOVAR_DB_FOLDER = "/glob/jessada/tools/annovar/humandb"
DFLT_ANNOVAR_DB_NAMES = "refGene"
DFLT_ANNOVAR_DB_OPS = "g"
DFLT_ANNOVAR_DB_NAMES += ",cytoBand"
DFLT_ANNOVAR_DB_OPS += ",r"
DFLT_ANNOVAR_DB_NAMES += ",genomicSuperDups"
DFLT_ANNOVAR_DB_OPS += ",r"
DFLT_ANNOVAR_DB_NAMES += ",esp6500siv2_all"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",1000g2014oct_eur"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",1000g2014oct_all"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",snp138"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",exac03"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",cg69"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",cosmic70"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",nci60"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",clinvar_20150330"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",mitimpact2"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_Ill_Br"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_101CRC_all"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_101CRC_CAFAM"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_101CRC_COLON"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_101CRC_CRC"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_101CRC_RECTAL"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_PMS2_sporadic_fam24"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_Axeq_chr3_6_14_18"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_Axeq_chr5_19"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_Axeq_chr9"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",200_Danes"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",249_Swedes"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",ljb26_all"
DFLT_ANNOVAR_DB_OPS += ",f"
