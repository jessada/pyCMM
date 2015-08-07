import os

"""

This module is all about constant value 

"""

#FULL_SYSTEM_TEST = True
FULL_SYSTEM_TEST = False
SLURM_TEST = False

FAST_PROJECT_CODE = 'b2011158'
SLOW_PROJECT_CODE = 'b2011097'

## > > > > > > > > > > > > > GATK Best practice (DNA seq) < < < < < < < < < <
DNASEQ_PIPELINE_DESCRIPTION = "A flow to control a pipeline to process DNA sequencing data"
DNASEQ_PIPELINE_DFLT_LOG_FILE = "dnaseq_pipeline"
GATK_ALLOC_TIME = "10:00:00"
#GATK_ALLOC_TIME = "4-00:00:00"

DNASEQ_CREATE_JOB_SETUP_FILE_DESCRIPTION = "An applicationt to generate job setup file for DNASEQ_PIPELINE"

## > > > > > > > > > > > > > CMM DB < < < < < < < < < <
CMMDB_PIPELINE_DESCRIPTION = "A flow to control a pipeline to process CMM database"
CMMDB_PIPELINE_DFLT_LOG_FILE = "cmmdb_pipeline"
CMMDB_ALLOC_TIME = "6:00:00"

CMMDB_CREATE_JOB_SETUP_FILE_DESCRIPTION = "An applicationt to generate job setup file to process CMM database"

## > > > > > > > > > > > > > Mutation report < < < < < < < < < <
MUTREP_PIPELINE_DESCRIPTION = "A flow to control a pipeline to process mutations report"
MuTREP_PIPELINE_DFLT_LOG_FILE = "mutrep_pipeline"
MUTREP_ALLOC_TIME = "1-00:00:00"

MUTREP_CREATE_JOB_SETUP_FILE_DESCRIPTION = "An applicationt to generate job setup file to process mutations report"

DFLT_MUTREP_FREQ_RATIOS = "esp6500siv2_all:0.2"

DFLT_MUTREP_ANNO_COLS = []
#DFLT_MUTREP_ANNO_COLS.append("Chr")
#DFLT_MUTREP_ANNO_COLS.append("Start")
#DFLT_MUTREP_ANNO_COLS.append("End")
#DFLT_MUTREP_ANNO_COLS.append("Ref")
#DFLT_MUTREP_ANNO_COLS.append("Alt")
DFLT_MUTREP_ANNO_COLS.append("Func.refGene")
DFLT_MUTREP_ANNO_COLS.append("Gene.refGene")
DFLT_MUTREP_ANNO_COLS.append("GeneDetail.refGene")
DFLT_MUTREP_ANNO_COLS.append("ExonicFunc.refGene")
DFLT_MUTREP_ANNO_COLS.append("1000g2014oct_eur")
DFLT_MUTREP_ANNO_COLS.append("1000g2014oct_all")
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
DFLT_MUTREP_ANNO_COLS.append("snp138")
DFLT_MUTREP_ANNO_COLS.append("cg69")
DFLT_MUTREP_ANNO_COLS.append("cosmic70")
DFLT_MUTREP_ANNO_COLS.append("nci60")
DFLT_MUTREP_ANNO_COLS.append("clinvar_20150330")
DFLT_MUTREP_ANNO_COLS.append("Gene_symbol")
DFLT_MUTREP_ANNO_COLS.append("OXPHOS_Complex")
DFLT_MUTREP_ANNO_COLS.append("Ensembl_Gene_ID")
DFLT_MUTREP_ANNO_COLS.append("Ensembl_Protein_ID")
DFLT_MUTREP_ANNO_COLS.append("Uniprot_Name")
DFLT_MUTREP_ANNO_COLS.append("Uniprot_ID")
DFLT_MUTREP_ANNO_COLS.append("NCBI_Gene_ID")
DFLT_MUTREP_ANNO_COLS.append("NCBI_Protein_ID")
DFLT_MUTREP_ANNO_COLS.append("Gene_pos")
DFLT_MUTREP_ANNO_COLS.append("AA_pos")
DFLT_MUTREP_ANNO_COLS.append("AA_sub")
DFLT_MUTREP_ANNO_COLS.append("Codon_sub")
DFLT_MUTREP_ANNO_COLS.append("dbSNP_ID")
DFLT_MUTREP_ANNO_COLS.append("PhyloP_46V")
DFLT_MUTREP_ANNO_COLS.append("PhastCons_46V")
DFLT_MUTREP_ANNO_COLS.append("PhyloP_100V")
DFLT_MUTREP_ANNO_COLS.append("PhastCons_100V")
DFLT_MUTREP_ANNO_COLS.append("SiteVar")
DFLT_MUTREP_ANNO_COLS.append("PolyPhen2_prediction")
DFLT_MUTREP_ANNO_COLS.append("PolyPhen2_score")
DFLT_MUTREP_ANNO_COLS.append("SIFT_prediction")
DFLT_MUTREP_ANNO_COLS.append("SIFT_score")
DFLT_MUTREP_ANNO_COLS.append("FatHmm_prediction")
DFLT_MUTREP_ANNO_COLS.append("FatHmm_score")
DFLT_MUTREP_ANNO_COLS.append("PROVEAN_prediction")
DFLT_MUTREP_ANNO_COLS.append("PROVEAN_score")
DFLT_MUTREP_ANNO_COLS.append("MutAss_prediction")
DFLT_MUTREP_ANNO_COLS.append("MutAss_score")
DFLT_MUTREP_ANNO_COLS.append("EFIN_Swiss_Prot_Score")
DFLT_MUTREP_ANNO_COLS.append("EFIN_Swiss_Prot_Prediction")
DFLT_MUTREP_ANNO_COLS.append("EFIN_HumDiv_Score")
DFLT_MUTREP_ANNO_COLS.append("EFIN_HumDiv_Prediction")
DFLT_MUTREP_ANNO_COLS.append("CADD_score")
DFLT_MUTREP_ANNO_COLS.append("CADD_Phred_score")
DFLT_MUTREP_ANNO_COLS.append("CADD_prediction")
DFLT_MUTREP_ANNO_COLS.append("Carol_prediction")
DFLT_MUTREP_ANNO_COLS.append("Carol_score")
DFLT_MUTREP_ANNO_COLS.append("Condel_score")
DFLT_MUTREP_ANNO_COLS.append("Condel_pred")
DFLT_MUTREP_ANNO_COLS.append("COVEC_WMV")
DFLT_MUTREP_ANNO_COLS.append("COVEC_WMV_prediction")
DFLT_MUTREP_ANNO_COLS.append("PolyPhen2_score_transf")
DFLT_MUTREP_ANNO_COLS.append("PolyPhen2_pred_transf")
DFLT_MUTREP_ANNO_COLS.append("SIFT_score_transf")
DFLT_MUTREP_ANNO_COLS.append("SIFT_pred_transf")
DFLT_MUTREP_ANNO_COLS.append("MutAss_score_transf")
DFLT_MUTREP_ANNO_COLS.append("MutAss_pred_transf")
DFLT_MUTREP_ANNO_COLS.append("Perc_coevo_Sites")
DFLT_MUTREP_ANNO_COLS.append("Mean_MI_score")
DFLT_MUTREP_ANNO_COLS.append("COSMIC_ID")
DFLT_MUTREP_ANNO_COLS.append("Tumor_site")
DFLT_MUTREP_ANNO_COLS.append("Examined_samples")
DFLT_MUTREP_ANNO_COLS.append("Mutation_frequency")
DFLT_MUTREP_ANNO_COLS.append("US")
DFLT_MUTREP_ANNO_COLS.append("Status")
DFLT_MUTREP_ANNO_COLS.append("Associated_disease")
DFLT_MUTREP_ANNO_COLS.append("Presence_in_TD")
DFLT_MUTREP_ANNO_COLS.append("Class_predicted")
DFLT_MUTREP_ANNO_COLS.append("Prob_N")
DFLT_MUTREP_ANNO_COLS.append("Prob_P")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_WT")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_HET")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_HOM")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_OTH")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_NA")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_GT")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_GF")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_AF")
DFLT_MUTREP_ANNO_COLS.append("ILL_BR_PF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_WT")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_HET")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_HOM")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_OTH")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_NA")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_GT")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_GF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_AF")
DFLT_MUTREP_ANNO_COLS.append("AXEQ_BR_PF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_WT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_HET")
DFLT_MUTREP_ANNO_COLS.append("101CRC_HOM")
DFLT_MUTREP_ANNO_COLS.append("101CRC_OTH")
DFLT_MUTREP_ANNO_COLS.append("101CRC_NA")
DFLT_MUTREP_ANNO_COLS.append("101CRC_GT")
DFLT_MUTREP_ANNO_COLS.append("101CRC_GF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_AF")
DFLT_MUTREP_ANNO_COLS.append("101CRC_PF")
DFLT_MUTREP_ANNO_COLS.append("SIFT_score")
DFLT_MUTREP_ANNO_COLS.append("SIFT_pred")
DFLT_MUTREP_ANNO_COLS.append("Polyphen2_HDIV_score")
DFLT_MUTREP_ANNO_COLS.append("Polyphen2_HDIV_pred")
DFLT_MUTREP_ANNO_COLS.append("Polyphen2_HVAR_score")
DFLT_MUTREP_ANNO_COLS.append("Polyphen2_HVAR_pred")
DFLT_MUTREP_ANNO_COLS.append("LRT_score")
DFLT_MUTREP_ANNO_COLS.append("LRT_pred")
DFLT_MUTREP_ANNO_COLS.append("MutationTaster_score")
DFLT_MUTREP_ANNO_COLS.append("MutationTaster_pred")
DFLT_MUTREP_ANNO_COLS.append("MutationAssessor_score")
DFLT_MUTREP_ANNO_COLS.append("MutationAssessor_pred")
DFLT_MUTREP_ANNO_COLS.append("FATHMM_score")
DFLT_MUTREP_ANNO_COLS.append("FATHMM_pred")
DFLT_MUTREP_ANNO_COLS.append("RadialSVM_score")
DFLT_MUTREP_ANNO_COLS.append("RadialSVM_pred")
DFLT_MUTREP_ANNO_COLS.append("LR_score")
DFLT_MUTREP_ANNO_COLS.append("LR_pred")
DFLT_MUTREP_ANNO_COLS.append("VEST3_score")
DFLT_MUTREP_ANNO_COLS.append("CADD_raw")
DFLT_MUTREP_ANNO_COLS.append("CADD_phred")
DFLT_MUTREP_ANNO_COLS.append("GERP++_RS")
DFLT_MUTREP_ANNO_COLS.append("phyloP46way_placental")
DFLT_MUTREP_ANNO_COLS.append("phyloP100way_vertebrate")
DFLT_MUTREP_ANNO_COLS.append("SiPhy_29way_logOdds")

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
DFLT_ANNOVAR_DB_NAMES += ",CMM_Axeq_Br"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",CMM_101CRC"
DFLT_ANNOVAR_DB_OPS += ",f"
DFLT_ANNOVAR_DB_NAMES += ",ljb26_all"
DFLT_ANNOVAR_DB_OPS += ",f"
