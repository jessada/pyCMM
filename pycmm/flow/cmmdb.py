import sys
import re
import datetime
import pyaml
import yaml
from os import listdir
from os.path import join as join_path
from os.path import isdir
from os.path import isfile
from pycmm.settings import DFLT_ANNOVAR_DB_FOLDER
from pycmm.settings import DFLT_ANNOVAR_DB_NAMES
from pycmm.settings import DFLT_ANNOVAR_DB_OPS
from pycmm.utils import exec_sh
from pycmm.utils import mylogger
from pycmm.utils.jobman import JobManager
from pycmm.proc.annovar import Annovar
from pycmm.proc.annovar import ANNOVAR_PARAMS_INPUT_FILE_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_FOLDER_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_BUILDVER_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_OUT_PREFIX_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_NAMES_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_DB_OPS_KEY
from pycmm.proc.annovar import ANNOVAR_PARAMS_NASTRING_KEY
from pycmm.settings import CMMDB_ALLOC_TIME
from pycmm.settings import DFLT_MUTREP_ANNO_COLS
from pycmm.settings import DFLT_MUTREP_FREQ_RATIOS

ALL_CHROMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "Y"]
CAL_MUTATIONS_STAT_SCRIPT = "$PYCMM/bash/cal_mutations_stat.sh"

# *************** jobs metadata section ***************
JOBS_SETUP_DATASET_NAME_KEY = "DATASET_NAME"
JOBS_SETUP_PROJECT_CODE_KEY = "PROJECT_CODE"
JOBS_SETUP_OUTPUT_DIR_KEY = "OUTPUT_DIR"
JOBS_SETUP_JOBS_REPORT_FILE_KEY = "JOBS_REPORT_FILE"
JOBS_SETUP_VCF_TABIX_FILE_KEY = "VCF_TABIX_FILE"
JOBS_SETUP_VCF_REGION_KEY = "VCF_REGION"

# *************** jobs ANNOVAR section ***************
JOBS_SETUP_ANNOVAR_SECTION = "ANNOVAR"
JOBS_SETUP_ANNOVAR_DB_FOLDER_KEY = ANNOVAR_PARAMS_DB_FOLDER_KEY
JOBS_SETUP_ANNOVAR_BUILDVER_KEY = ANNOVAR_PARAMS_BUILDVER_KEY
# >> table_annovar section
JOBS_SETUP_TA_DB_NAMES_KEY = ANNOVAR_PARAMS_DB_NAMES_KEY
JOBS_SETUP_TA_DB_OPS_KEY = ANNOVAR_PARAMS_DB_OPS_KEY
JOBS_SETUP_TA_NASTRING_KEY = ANNOVAR_PARAMS_NASTRING_KEY

# *************** jobs samples section ***************
JOBS_SETUP_SAMPLES_INFO_KEY = "SAMPLES_INFO"
JOBS_SETUP_FAMILY_ID_KEY = "FAMILY_ID"
JOBS_SETUP_MEMBERS_LIST_KEY = "MEMBERS"
JOBS_SETUP_MEMBER_NAME_KEY = "NAME"

# *************** report layout section ***************
JOBS_SETUP_REPORT_LAYOUT_SECTION = "REPORT_LAYOUT"
JOBS_SETUP_REPORT_ANNO_COLS_KEY = "COLUMNS"
JOBS_SETUP_REPORT_REGIONS_KEY = "REGIONS"
JOBS_SETUP_REPORT_CALL_INFO_KEY = "CALL_INFO"
JOBS_SETUP_REPORT_FREQ_RATIOS_KEY = "FREQUENCY_RATIOS"
JOBS_SETUP_REPORT_FREQ_RATIOS_COL_KEY = "COLUMN"
JOBS_SETUP_REPORT_FREQ_RATIOS_FREQ_KEY = "FREQUENCY"

class CMMDBPipeline(JobManager):
    """ A class to control CMM data(base) workflow """

    def __init__(self,
                 jobs_setup_file,
                 ):
        mylogger.getLogger(__name__)
        self.__load_jobs_info(jobs_setup_file)
        JobManager.__init__(self,
                            jobs_report_file=self.jobs_report_file)
        self.__parse_samples_info()
        self.__time_stamp = datetime.datetime.now()
        self.__create_directories()

    def get_raw_repr(self):
        return {"dataset name": self.dataset_name,
                "project code": self.project_code,
                "jobs report file": self.jobs_report_file,
                }

    @property
    def dataset_name(self):
        return self._jobs_info[JOBS_SETUP_DATASET_NAME_KEY]

    @property
    def project_code(self):
        if JOBS_SETUP_PROJECT_CODE_KEY in self._jobs_info:
            return self._jobs_info[JOBS_SETUP_PROJECT_CODE_KEY]
        else:
            return None

    @property
    def input_vcf_tabix(self):
        return self._jobs_info[JOBS_SETUP_VCF_TABIX_FILE_KEY]

    @property
    def vcf_region(self):
        if JOBS_SETUP_VCF_REGION_KEY in self._jobs_info:
            return str(self._jobs_info[JOBS_SETUP_VCF_REGION_KEY])
        else:
            return None

    @property
    def out_stat_file(self):
        return self.__out_stat_file

    @property
    def samples_list(self):
        return self.__samples_list

    @property
    def families_list(self):
        return self.__families_list

    @property
    def jobs_report_file(self):
        return join_path(self.output_dir,
                         self.dataset_name+"_rpt.txt")

    @property
    def output_dir(self):
        return self._jobs_info[JOBS_SETUP_OUTPUT_DIR_KEY]

    @property
    def rpts_out_dir(self):
        return join_path(self.output_dir,
                         "rpts")

    @property
    def data_out_dir(self):
        return join_path(self.output_dir,
                         "data_out")

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
        cfg = self._jobs_info[JOBS_SETUP_ANNOVAR_SECTION]
        cfg[ANNOVAR_PARAMS_OUT_PREFIX_KEY] = join_path(self.rpts_out_dir,
                                                       self.dataset_name)
        cfg[ANNOVAR_PARAMS_INPUT_FILE_KEY] = self.input_vcf_tabix
        return cfg

    @property
    def report_layout(self):
        return self._jobs_info[JOBS_SETUP_REPORT_LAYOUT_SECTION]

    def __load_jobs_info(self, jobs_setup_file):
        self.__meta_data = {}
        stream = file(jobs_setup_file, "r")
        self._jobs_info = yaml.safe_load(stream)

    def __parse_samples_info(self):
        if JOBS_SETUP_SAMPLES_INFO_KEY not in self._jobs_info:
            self.__samples_list = None
            self.__families_list = None
            return
        samples_info = self._jobs_info[JOBS_SETUP_SAMPLES_INFO_KEY]
        if type(samples_info) is str:
            self.__samples_list = samples_info.split(",")
            self.__families_list = None
        else:
            self.__samples_list = []
            for family_info in samples_info:
                family_info[JOBS_SETUP_FAMILY_ID_KEY] = str(family_info[JOBS_SETUP_FAMILY_ID_KEY])
                for member in family_info[JOBS_SETUP_MEMBERS_LIST_KEY]:
                    self.__samples_list.append(member[JOBS_SETUP_MEMBER_NAME_KEY])
            self.__families_list = samples_info

    def __create_directories(self):
        self.create_dir(self.rpts_out_dir)
        self.create_dir(self.data_out_dir)
        self.create_dir(self.slurm_log_dir)
        self.create_dir(self.working_dir)

    def __get_cal_mut_stat_params(self,
                                  dataset_name,
                                  out_stat_file,
                                  vcf_region=None,
                                  samples_list=None,
                                  ):
        params = " -k " + dataset_name
        params += " -i " + self.input_vcf_tabix
        if vcf_region is not None:
            params += " -r " + vcf_region
        if samples_list is not None:
            params += " -c " + samples_list
        params += " -o " + out_stat_file
        return params

    def cal_mut_stat(self):
        if self.project_code is None:
            self.__out_stat_file = join_path(self.data_out_dir,
                                             self.dataset_name + ".stat")
            params = self.__get_cal_mut_stat_params(dataset_name=self.dataset_name,
                                                    out_stat_file=self.out_stat_file,
                                                    vcf_region=self.vcf_region,
                                                    samples_list=self.samples_list,
                                                    )
            cmd = CAL_MUTATIONS_STAT_SCRIPT + params
            exec_sh(cmd)
        elif self.vcf_region is None:
            self.__out_stat_file = []
            for chrom in ALL_CHROMS:
                dataset_name = self.dataset_name + "_" + chrom
                out_stat_file = join_path(self.data_out_dir,
                                          dataset_name + ".stat")
                params = self.__get_cal_mut_stat_params(dataset_name=dataset_name,
                                                        out_stat_file=out_stat_file,
                                                        vcf_region=chrom,
                                                        samples_list=self.samples_list,
                                                        )
                self.__out_stat_file.append(out_stat_file)
                job_name = dataset_name + "_cal_stat"
                slurm_log_file = join_path(self.slurm_log_dir,
                                           job_name)
                slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
                slurm_log_file += ".log"
                self.submit_job(job_name,
                                self.project_code,
                                "core",
                                "1",
                                CMMDB_ALLOC_TIME,
                                slurm_log_file,
                                CAL_MUTATIONS_STAT_SCRIPT,
                                params,
                                )
        else:
            self.__out_stat_file = join_path(self.data_out_dir,
                                             self.dataset_name + ".stat")
            params = self.__get_cal_mut_stat_params(dataset_name=self.dataset_name,
                                                    out_stat_file=self.out_stat_file,
                                                    vcf_region=self.vcf_region,
                                                    samples_list=self.samples_list,
                                                    )
            job_name = self.dataset_name + "_cal_stat"
            slurm_log_file = join_path(self.slurm_log_dir,
                                       job_name)
            slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
            slurm_log_file += ".log"
            self.submit_job(job_name,
                            self.project_code,
                            "core",
                            "1",
                            CMMDB_ALLOC_TIME,
                            slurm_log_file,
                            CAL_MUTATIONS_STAT_SCRIPT,
                            params,
                            )

    def table_annovar(self):
        annovar = Annovar(self.annovar_config)
        if self.project_code is not None:
            job_name = self.dataset_name + "_ta"
            slurm_log_file = join_path(self.slurm_log_dir,
                                       job_name)
            slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
            slurm_log_file += ".log"
            self.submit_job(job_name,
                            self.project_code,
                            "core",
                            "1",
                            CMMDB_ALLOC_TIME,
                            slurm_log_file,
                            annovar.table_annovar_cmd,
                            "",
                            )
        else:
            exec_sh(annovar.table_annovar_cmd)

    def __garbage_collecting(self):
        pass

    def monitor_action(self):
        JobManager.monitor_action(self)
        self.__garbage_collecting()

def create_jobs_setup_file(dataset_name,
                           project_out_dir,
                           vcf_tabix_file,
                           vcf_region=None,
                           samples_info=None,
                           project_code=None,
                           jobs_report_file=None,
                           annovar_human_db_dir=DFLT_ANNOVAR_DB_FOLDER,
                           annovar_buildver="hg19",
                           annovar_db_names=DFLT_ANNOVAR_DB_NAMES,
                           annovar_db_ops=DFLT_ANNOVAR_DB_OPS,
                           annovar_nastring=".",
                           anno_cols=",".join(DFLT_MUTREP_ANNO_COLS),
                           report_regions=None,
                           call_info="NO",
                           frequency_ratios=DFLT_MUTREP_FREQ_RATIOS,
                           out_jobs_setup_file=None,
                           ):
    mylogger.getLogger(__name__)
    if jobs_report_file is None:
        jobs_report_file = join_path(project_out_dir,
                                     dataset_name+"_rpt.txt")
    if out_jobs_setup_file is None:
        out_jobs_setup_file = join_path(samples_root_dir,
                                       dataset_name+"_job_setup.txt")
    annovar_config = {}
    annovar_config[JOBS_SETUP_ANNOVAR_DB_FOLDER_KEY] = annovar_human_db_dir
    annovar_config[JOBS_SETUP_ANNOVAR_BUILDVER_KEY] = annovar_buildver
    annovar_config[JOBS_SETUP_TA_DB_NAMES_KEY] = annovar_db_names
    annovar_config[JOBS_SETUP_TA_DB_OPS_KEY] = annovar_db_ops
    annovar_config[JOBS_SETUP_TA_NASTRING_KEY] = annovar_nastring
    stream = file(out_jobs_setup_file, 'w')
    job_setup_document = {}
    job_setup_document[JOBS_SETUP_DATASET_NAME_KEY] = dataset_name
    if project_code is not None:
        job_setup_document[JOBS_SETUP_PROJECT_CODE_KEY] = project_code
    if vcf_region is not None:
        job_setup_document[JOBS_SETUP_VCF_REGION_KEY] = vcf_region
    if (samples_info is not None) and (samples_info.find(":") == -1):
        job_setup_document[JOBS_SETUP_SAMPLES_INFO_KEY] = samples_info
    if (samples_info is not None) and (samples_info.find(":") > -1):
        families_info = []
        for family_item in samples_info.split(","):
            family_info = {}
            family_infos = family_item.split(":")
            family_info[JOBS_SETUP_FAMILY_ID_KEY] = family_infos[0]
            members = []
            for member_infos in family_infos[1:]:
                member = {}
                member[JOBS_SETUP_MEMBER_NAME_KEY] = member_infos
                members.append(member)
            family_info[JOBS_SETUP_MEMBERS_LIST_KEY] = members
            families_info.append(family_info)
        job_setup_document[JOBS_SETUP_SAMPLES_INFO_KEY] = families_info
    job_setup_document[JOBS_SETUP_OUTPUT_DIR_KEY] = project_out_dir
    job_setup_document[JOBS_SETUP_VCF_TABIX_FILE_KEY] = vcf_tabix_file
    job_setup_document[JOBS_SETUP_JOBS_REPORT_FILE_KEY] = jobs_report_file
    job_setup_document[JOBS_SETUP_ANNOVAR_SECTION] = annovar_config
    report_layout_config = {}
    report_layout_config[JOBS_SETUP_REPORT_ANNO_COLS_KEY] = anno_cols.split(",")
    if report_regions is not None:
        report_layout_config[JOBS_SETUP_REPORT_REGIONS_KEY] = report_regions.split(",")
    report_layout_config[JOBS_SETUP_REPORT_CALL_INFO_KEY] = call_info
    if frequency_ratios is not None:
        job_freq_ratios = []
        for frequency_ratio in frequency_ratios.split(","):
            (col, freq) = frequency_ratio.split(":")
            job_freq_ratios.append({JOBS_SETUP_REPORT_FREQ_RATIOS_COL_KEY: col,
                                    JOBS_SETUP_REPORT_FREQ_RATIOS_FREQ_KEY: freq,
                                    })
        report_layout_config[JOBS_SETUP_REPORT_FREQ_RATIOS_KEY] = job_freq_ratios
    job_setup_document[JOBS_SETUP_REPORT_LAYOUT_SECTION] = report_layout_config

    pyaml.dump(job_setup_document, stream)
