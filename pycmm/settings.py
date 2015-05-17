#import os

"""

This module is all about constant value 

"""

#FULL_SYSTEM_TEST = True
FULL_SYSTEM_TEST = False
SLURM_TEST = False

DFLT_WORKING_DIR = "."
VCF_MUT_KEY_INT_FMT = "{chrom:02d}_{pos:012d}_{ref:s}_{alt:s}"
VCF_MUT_KEY_STR_FMT = "{chrom:s}_{pos:012d}_{ref:s}_{alt:s}"

## > > > > > > > > > > > > > vcf2avdb < < < < < < < < < <
VCF2AVDB_DESCRIPTION = "A script to parse all point mutations from vcf format into avdb format"
VCF2AVDB_DFLT_LOG_FILE = "vcf2avdb"

## > > > > > > > > > > > > > GATK Best practice (DNA seq) < < < < < < < < < <
GATK_ALLOC_TIME = "72:00:00"

## > > > > > > > > > > > > > random values < < < < < < < < < <
#DFLT_SEED = None
#DEMO_SEED = 20
#TEST_SEED = DEMO_SEED
#
## > > > > > > > > > > > > > development files & folders < < < < < < < < < <
#PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))
#
#OUT_DIR = '/home/jessada/development/projects/COMPAS/out/'
#FIGS_TMP_OUT_DIR = os.path.join(OUT_DIR, 'figs_tmp')
#TERM_TMP_OUT_DIR = os.path.join(OUT_DIR, 'term_tmp')
#
## > > > > > > > > > > > > > model parameters < < < < < < < < < <
#DFLT_MAP_SIZE = 100
#DFLT_WEIGHT_STEP_SIZE = 0.2
#DFLT_NBH_STEP_SIZE = 1
#DFLT_MAX_NBH_SIZE = 10
#DFLT_MAP_ROWS = 17
#DFLT_MAP_COLS = 17
#
## > > > > > > > > > > > > > samples type < < < < < < < < < <
#TYPE_TRAINING_SAMPLE = 'training'
#TYPE_TEST_SAMPLE = 'test'
