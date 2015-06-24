import sys
import re
import datetime
import pyaml
import yaml
from os import listdir
from os.path import join as join_path
from os.path import isdir
from os.path import isfile
#from pycmm.template import pyCMMBase
from pycmm.utils import mylogger
from pycmm.utils.jobman import JobManager
from pycmm.proc.annovar import Annovar
from pycmm.proc.annovar import ANNOVAR_PARAMS_INPUT_FILE_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_FOLDER_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_BUILDVER_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_OUT_PREFIX_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_LIST_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_NAME_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_OP_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_NASTRING_KEY
#from pycmm.utils.jobman import JOB_STATUS_COMPLETED
#from pycmm.utils import exec_sh
#from pycmm.settings import GATK_ALLOC_TIME

# *************** jobs metadata section ***************
JOBS_SETUP_DATASET_NAME_KEY = "DATASET_NAME"
JOBS_SETUP_PROJECT_CODE_KEY = "PROJECT_CODE"
JOBS_SETUP_OUTPUT_DIR_KEY = "OUTPUT_DIR"
JOBS_SETUP_JOBS_REPORT_FILE_KEY = "JOBS_REPORT_FILE"
JOBS_SETUP_VCF_TABIX_FILE_KEY = "VCF_TABIX_FILE"

# *************** ANNOVAR section ***************
JOBS_SETUP_ANNOVAR_SECTION = "ANNOVAR"
JOBS_SETUP_ANNOVAR_DB_FOLDER_KEY = ANNOVAR_PARAMS_DB_FOLDER_KEY
JOBS_SETUP_ANNOVAR_BUILDVER_KEY = ANNOVAR_PARAMS_BUILDVER_KEY
# >> table_annovar section
#JOBS_SETUP_TA_ANNOVAR_SECTION = "TABLE_ANNOVAR"
#JOBS_SETUP_TA_INPUT_FILE_KEY = ANNOVAR_PARAMS_INPUT_FILE_KEY
JOBS_SETUP_TA_DB_LIST_KEY = ANNOVAR_PARAMS_DB_LIST_KEY
JOBS_SETUP_TA_DB_NAME_KEY = ANNOVAR_PARAMS_DB_NAME_KEY
JOBS_SETUP_TA_DB_OP_KEY = ANNOVAR_PARAMS_DB_OP_KEY
JOBS_SETUP_TA_NASTRING_KEY = ANNOVAR_PARAMS_NASTRING_KEY
#JOBS_SETUP_TA_OUT_PREFIX = ANNOVAR_PARAMS_OUT_PREFIX_KEY

## jobs sample info section
#JOBS_SETUP_SAMPLES_SECTION = "SAMPLES"
#JOBS_SETUP_SAMPLE_NAME_KEY = "sample_name"

class MutRepPipeline(JobManager):
    """ A class to control GATK best practice pipeline """

    def __init__(self,
                 jobs_setup_file,
                 ):
        mylogger.getLogger(__name__)
        self.__load_jobs_info(jobs_setup_file)
        JobManager.__init__(self,
                            jobs_report_file=self.jobs_report_file)
        self.__create_directories()

    def get_raw_repr(self):
        return {"dataset name": self.dataset_name,
                "project code": self.project_code,
                "jobs report file": self.jobs_report_file,
                }

    @property
    def dataset_name(self):
        return self.__jobs_info[JOBS_SETUP_DATASET_NAME_KEY]

    @property
    def project_code(self):
        return self.__jobs_info[JOBS_SETUP_PROJECT_CODE_KEY]

    @property
    def input_vcf_tabix(self):
        return self.__jobs_info[JOBS_SETUP_VCF_TABIX_FILE_KEY]

    @property
    def jobs_report_file(self):
        return join_path(self.output_dir,
                         self.dataset_name+"_rpt.txt")

    @property
    def output_dir(self):
        return self.__jobs_info[JOBS_SETUP_OUTPUT_DIR_KEY]

    @property
    def rpts_out_dir(self):
        return join_path(self.output_dir,
                         "rpts")

    @property
    def slurm_log_dir(self):
        return join_path(self.output_dir,
                         "slurm_log")

    @property
    def working_dir(self):
        return join_path(self.output_dir,
                         "tmp")

    @property
    def time_stamp(self):
        return self.__time_stamp

    @property
    def annovar_config(self):
        cfg = self.__jobs_info[JOBS_SETUP_ANNOVAR_SECTION]
        cfg[ANNOVAR_PARAMS_OUT_PREFIX_KEY] = join_path(self.rpts_out_dir,
                                                       self.dataset_name)
        cfg[ANNOVAR_PARAMS_INPUT_FILE_KEY] = self.input_vcf_tabix
        return cfg

    def __load_jobs_info(self, jobs_setup_file):
        self.__meta_data = {}
        stream = file(jobs_setup_file, "r")
        self.__jobs_info = yaml.safe_load(stream)

    def __create_directories(self):
        self.create_dir(self.rpts_out_dir)
        self.create_dir(self.slurm_log_dir)
        self.create_dir(self.working_dir)

    def __garbage_collecting(self):
        pass

    def monitor_action(self):
        JobManager.monitor_action(self)
        pass
