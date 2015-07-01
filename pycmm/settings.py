import os

"""

This module is all about constant value 

"""

#FULL_SYSTEM_TEST = True
FULL_SYSTEM_TEST = False
SLURM_TEST = False

FAST_PROJECT_CODE = 'b2011158'
SLOW_PROJECT_CODE = 'b2011097'

#DFLT_WORKING_DIR = "."
#VCF_MUT_KEY_INT_FMT = "{chrom:02d}_{pos:012d}_{ref:s}_{alt:s}"
#VCF_MUT_KEY_STR_FMT = "{chrom:s}_{pos:012d}_{ref:s}_{alt:s}"

## > > > > > > > > > > > > > vcf2avdb < < < < < < < < < <
#VCF2AVDB_DESCRIPTION = "A script to parse all point mutations from vcf format into avdb format"
#VCF2AVDB_DFLT_LOG_FILE = "vcf2avdb"

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
