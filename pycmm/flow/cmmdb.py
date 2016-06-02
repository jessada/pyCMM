import pyaml
#import sys
#import re
#import datetime
#import pysam
#import yaml
#from os import listdir
from os.path import join as join_path
#from os.path import isdir
#from os.path import isfile
#from collections import defaultdict
from pycmm.utils import exec_sh
from pycmm.cmmlib import CMMParams
from pycmm.flow import CMMPipeline
from pycmm.flow import init_jobs_setup_file
from pycmm.flow import get_func_arg
from pycmm.cmmlib.annovarlib import AnnovarParams
from pycmm.cmmlib.annovarlib import get_annovar_params_sections
#from pycmm.settings import DFLT_CMMDB_ALLOC_TIME
#from pycmm.settings import DFLT_MUTREP_ALLOC_TIME
#from pycmm.settings import DUMMY_TABLE_ANNOVAR_BIN
#from pycmm.settings import ALL_MUTREP_ANNO_COLS
#from pycmm.template import pyCMMBase
#from pycmm.utils import mylogger
#from pycmm.utils.jobman import JobManager
#from pycmm.cmmlib.dnalib import ALL_CHROMS
#from pycmm.cmmlib.annovarlib import Annovar
#from pycmm.cmmlib.annovarlib import AnnovarParams
#from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_INPUT_FILE_KEY
#from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_DB_DIR_KEY
#from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_BUILDVER_KEY
#from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_OUT_PREFIX_KEY
#from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_DB_NAMES_KEY
#from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_DB_OPS_KEY
#from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_NASTRING_KEY
#
CAL_MUTATIONS_STAT_SCRIPT = "$PYCMM/bash/cal_mutations_stat.sh"
#
## *************** jobs metadata section ***************
#JOBS_SETUP_DATASET_NAME_KEY = "DATASET_NAME"
#JOBS_SETUP_PROJECT_CODE_KEY = "PROJECT_CODE"
#JOBS_SETUP_DB_ALLOC_TIME_KEY = "DB_ALLOC_TIME"
#JOBS_SETUP_RPT_ALLOC_TIME_KEY = "RPT_ALLOC_TIME"
#JOBS_SETUP_OUTPUT_DIR_KEY = "OUTPUT_DIR"
#JOBS_SETUP_JOBS_REPORT_FILE_KEY = "JOBS_REPORT_FILE"
#
## *************** ANNOVAR section ***************
#JOBS_SETUP_ANV_SECTION = "ANNOVAR"
#JOBS_SETUP_ANV_DB_DIR_KEY = ANNOVAR_PARAMS_DB_DIR_KEY
#JOBS_SETUP_ANV_BUILDVER_KEY = ANNOVAR_PARAMS_BUILDVER_KEY
## >> table_annovar section
#JOBS_SETUP_ANV_DB_NAMES_KEY = ANNOVAR_PARAMS_DB_NAMES_KEY
#JOBS_SETUP_ANV_DB_OPS_KEY = ANNOVAR_PARAMS_DB_OPS_KEY
#JOBS_SETUP_ANV_NASTRING_KEY = ANNOVAR_PARAMS_NASTRING_KEY
#
## *************** jobs samples section ***************
#JOBS_SETUP_SAMPLE_INFOS_KEY = "SAMPLE_INFOS"
#JOBS_SETUP_FAMILY_ID_KEY = "FAMILY_ID"
#JOBS_SETUP_MEMBERS_LIST_KEY = "MEMBERS"
#JOBS_SETUP_SAMPLE_ID_KEY = "SAMPLE_ID"
#

# *************** mustat db section ***************
JOBS_SETUP_MUTSTAT_PARAMS_SECTION = "MUTSTAT_PARAMS"
JOBS_SETUP_VCF_TABIX_FILE_KEY = "VCF_TABIX_FILE"
JOBS_SETUP_DB_REGION_KEY = "DB_REGION"

# *************** table_annovar db section ***************
JOBS_SETUP_ANV_PARAMS_SECTION = "ANNOVAR_PARAMS"

#class MemberInfo(pyCMMBase):
#    """  To encapsulate family information so that it is readable """
#
#    def __init__(self, info):
#        pyCMMBase.__init__(self)
#        self.__info = info
#
#    def get_raw_repr(self):
#        return {"sample id": self.sample_id,
#                }
#
#    @property
#    def sample_id(self):
#        return self.__info[JOBS_SETUP_SAMPLE_ID_KEY]
#
#class FamilyInfo(pyCMMBase):
#    """  To encapsulate family information so that it is readable """
#
#    def __init__(self, info):
#        pyCMMBase.__init__(self)
#        self.__info = info
#        self.__shared_mutations = None
#
#    def get_raw_repr(self):
#        out_txt = {}
#        out_txt["family id"] = self.fam_id
#        out_txt["members"] = self.members
#        if self.shared_mutations is not None:
#            out_txt["shared mutations"] = self.shared_mutations
#        return out_txt
#
#    @property
#    def fam_id(self):
#        return str(self.__info[JOBS_SETUP_FAMILY_ID_KEY])
#
#    @property
#    def members(self):
#        return map(lambda x: MemberInfo(x),
#                   self.__info[JOBS_SETUP_MEMBERS_LIST_KEY])
#
#    @property
#    def shared_mutations(self):
#        return self.__shared_mutations
#
#    @shared_mutations.setter
#    def shared_mutations(self, val):
#        self.__shared_mutations = val
#
class MutStatParams(CMMParams):
    """ To handle and parse mutation statistics parameters """

    def __init__(self, **kwargs):
        super(MutStatParams, self).__init__(**kwargs)

    def get_raw_repr(self, **kwargs):
        raw_repr = super(MutStatParams, self).get_raw_repr(**kwargs)
        raw_repr["input vcf tabix file"] = self.input_vcf_tabix
        return raw_repr

    @property
    def input_vcf_tabix(self):
        return self._get_job_config(JOBS_SETUP_VCF_TABIX_FILE_KEY,
                                    required=True)

    @property
    def db_region(self):
        return self._get_job_config(JOBS_SETUP_DB_REGION_KEY)

class CMMDBPipeline(CMMPipeline):
    """ A class to control CMMDB best practice pipeline """

    def __init__(self, **kwargs):
        super(CMMDBPipeline, self).__init__(**kwargs)
        self.__init_properties()

    def get_raw_repr(self, **kwargs):
        raw_repr = super(CMMDBPipeline, self).get_raw_repr(**kwargs)
        return raw_repr

    def __init_properties(self):
        pass

    @property
    def mutstat_params(self):
        return MutStatParams(entries=self._get_job_config(JOBS_SETUP_MUTSTAT_PARAMS_SECTION,
                                                          required=True))

    @property
    def annovar_params(self):
        return AnnovarParams(entries=self._get_job_config(JOBS_SETUP_ANV_PARAMS_SECTION,
                                                          required=True))

    @property
    def out_stat_files(self):
        return self.__out_stat_files

    def get_third_party_software_version(self):
        vm = VersionManager()
        versions = OrderedDict()
        versions['vcftools'] = vm.vcftools
        return versions

    def __get_cal_mut_stat_params(self,
                                  dataset_name,
                                  out_stat_file,
                                  db_region=None,
                                  samples_id=None,
                                  ):
        params = " -k " + dataset_name
        params += " -i " + self.mutstat_params.input_vcf_tabix
        if db_region is not None:
            params += " -r " + db_region
        if samples_id is not None:
            if type(samples_id) is list:
                params += " -c " + ",".join(samples_id)
            else:
                params += " -c " + samples_id
        params += " -o " + out_stat_file
        return params

    def cal_mut_stat(self):
        if self.project_code is None:
            self.__out_stat_files = join_path(self.data_out_dir,
                                              self.dataset_name + ".stat")
            params = self.__get_cal_mut_stat_params(dataset_name=self.dataset_name,
                                                    out_stat_file=self.out_stat_files,
                                                    db_region=self.mutstat_params.db_region,
                                                    samples_id=self.samples_id,
                                                    )
            cmd = CAL_MUTATIONS_STAT_SCRIPT + params
            exec_sh(cmd)
        elif self.mutstat_params.db_region is None:
            self.__out_stat_files = []
            for chrom in ALL_CHROMS:
                dataset_name = self.dataset_name + "_" + chrom
                out_stat_file = join_path(self.data_out_dir,
                                          dataset_name + ".stat")
                params = self.__get_cal_mut_stat_params(dataset_name=dataset_name,
                                                        out_stat_file=out_stat_file,
                                                        db_region=chrom,
                                                        samples_id=self.samples_id,
                                                        )
                self.__out_stat_files.append(out_stat_file)
                job_name = dataset_name + "_cal_stat"
                slurm_log_file = join_path(self.slurm_log_dir,
                                           job_name)
                slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
                slurm_log_file += ".log"
                self._submit_slurm_job(job_name,
                                       "1",
                                       CAL_MUTATIONS_STAT_SCRIPT,
                                       params,
                                       )
        else:
            self.__out_stat_files = join_path(self.data_out_dir,
                                              self.dataset_name + ".stat")
            params = self.__get_cal_mut_stat_params(dataset_name=self.dataset_name,
                                                    out_stat_file=self.out_stat_files,
                                                    db_region=self.mutstat_params.db_region,
                                                    samples_id=self.samples_id,
                                                    )
            job_name = self.dataset_name + "_cal_stat"
            self._submit_slurm_job(job_name,
                                   "1",
                                   CAL_MUTATIONS_STAT_SCRIPT,
                                   params,
                                   )

#class CMMDBPipeline(JobManager):
#    """ A class to control CMM data(base) workflow """
#
#    def __init__(self,
#                 jobs_setup_file,
#                 ):
#        self.__load_jobs_info(jobs_setup_file)
#        self.__jobs_setup_file = jobs_setup_file
#        JobManager.__init__(self,
#                            jobs_report_file=self.jobs_report_file)
#        self.__parse_sample_infos()
#        self.__annovar_configs = None
#        self.__time_stamp = datetime.datetime.now()
#        self.__create_directories()
#
#    def get_raw_repr(self):
#        return {"dataset name": self.dataset_name,
#                "project code": self.project_code,
#                "jobs report file": self.jobs_report_file,
#                }
#
#    @property
#    def samples_list_w_fam_pref(self):
#        return self.__samples_list_w_fam_pref
#
#    @property
#    def data_out_dir(self):
#        return join_path(self.output_dir,
#                         "data_out")
#
#    @property
#    def annovar_config(self):
#        if self.__annovar_configs is None:
#            self.__parse_annovar_configs()
#        return self.__annovar_configs
#
#    @property
#    def report_layout(self):
#        return self._jobs_info[JOBS_SETUP_RPT_LAYOUT_SECTION]
#
#    def table_annovar(self):
#        cfg = self.annovar_config
#        table_annovar_cmd = DUMMY_TABLE_ANNOVAR_BIN
#        table_annovar_cmd += " --dataset_name " + self.dataset_name
#        table_annovar_cmd += " --input_file " + cfg.input_file
#        table_annovar_cmd += " --db_folder " + cfg.db_folder
#        table_annovar_cmd += " --buildver " + cfg.buildver
#        table_annovar_cmd += " --protocols " + cfg.protocols
#        table_annovar_cmd += " --operations " + cfg.operations
#        table_annovar_cmd += " --nastring " + cfg.nastring
#        table_annovar_cmd += " --data_out_folder " + self.output_dir
#        if self.project_code is not None:
#            job_name = self.dataset_name + "_ta"
#            slurm_log_file = join_path(self.slurm_log_dir,
#                                       job_name)
#            slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
#            slurm_log_file += ".log"
## *********************************************************************************************** Need refactoring ***********************************************************************************************
#            self.submit_job(job_name,
#                            self.project_code,
#                            "core",
#                            "2",
#                            self.job_alloc_time,
#                            slurm_log_file,
#                            table_annovar_cmd,
#                            "",
#                            )
## *********************************************************************************************** Need refactoring ***********************************************************************************************
#        else:
#            exec_sh(table_annovar_cmd)
#
#    def __garbage_collecting(self):
#        pass
#
#    def monitor_action(self):
#        JobManager.monitor_action(self)
#        self.__garbage_collecting()

#def create_jobs_setup_file(project_name,
#                           project_out_dir,
#                           sample_info=None,
#                           project_code=None,
#                           flow_alloc_time=None,
#                           rpt_alloc_time=None,
#                           job_alloc_time=None,
#                           jobs_report_file=None,
#                           out_jobs_setup_file=None,
#                           vcf_tabix_file=None,
#                           db_region=None,
#                           annovar_human_db_dir=DFLT_ANV_DB_DIR,
#                           annovar_buildver="hg19",
#                           annovar_db_names=DFLT_ANV_DB_NAMES,
#                           annovar_db_ops=DFLT_ANV_DB_OPS,
#                           annovar_nastring=".",
#                           anno_cols=None,
#                           anno_excl_tags=None,
#                           rows_filter_actions=None,
#                           annotated_vcf_tabix=None,
#                           report_regions=None,
#                           frequency_ratios=None,
#                           expression_patterns=None,
#                           expression_usages=None,
#                           split_chrom=False,
#                           summary_families_sheet=False,
#                           call_detail=False,
#                           only_summary=False,
#                           only_families=False,
#                           ):

def create_cmmdb_jobs_setup_file(*args, **kwargs):
    job_setup_document, stream = init_jobs_setup_file(*args, **kwargs)
    mutstat_params = {}
    vcf_tabix_file = get_func_arg('vcf_tabix_file', kwargs)
    if vcf_tabix_file is not None:
        mutstat_params[JOBS_SETUP_VCF_TABIX_FILE_KEY] = vcf_tabix_file
    db_region = get_func_arg('db_region', kwargs)
    if db_region is not None:
        mutstat_params[JOBS_SETUP_DB_REGION_KEY] = '"' + db_region + '"'
    job_setup_document[JOBS_SETUP_MUTSTAT_PARAMS_SECTION] = mutstat_params

    anv_params = get_annovar_params_sections(*args, **kwargs)
    job_setup_document[JOBS_SETUP_ANV_PARAMS_SECTION] = anv_params
    return job_setup_document, stream

# Right now no one actually want to use create_cmmdb_jobs_setup_file
# But reserve it here in case it need to be a complete report and database pipeline
def create_jobs_setup_file(*args, **kwargs):
    job_setup_document, stream = create_cmmdb_jobs_setup_file(*args, **kwargs)
    pyaml.dump(job_setup_document, stream)
#    annovar_config = {}
#    annovar_config[JOBS_SETUP_ANV_DB_DIR_KEY] = annovar_human_db_dir
#    annovar_config[JOBS_SETUP_ANV_BUILDVER_KEY] = annovar_buildver
#    annovar_config[JOBS_SETUP_ANV_DB_NAMES_KEY] = annovar_db_names
#    annovar_config[JOBS_SETUP_ANV_DB_OPS_KEY] = annovar_db_ops
#    annovar_config[JOBS_SETUP_ANV_NASTRING_KEY] = annovar_nastring
#    stream = file(out_jobs_setup_file, 'w')
#    job_setup_document = {}
#    job_setup_document[JOBS_SETUP_DATASET_NAME_KEY] = dataset_name
#    if project_code is not None:
#        job_setup_document[JOBS_SETUP_PROJECT_CODE_KEY] = project_code
#    if job_alloc_time is None:
#        job_alloc_time = DFLT_CMMDB_ALLOC_TIME
#    job_setup_document[JOBS_SETUP_DB_ALLOC_TIME_KEY] = '"' + job_alloc_time + '"'
#    if rpt_alloc_time is None:
#        rpt_alloc_time = DFLT_MUTREP_ALLOC_TIME
#    job_setup_document[JOBS_SETUP_RPT_ALLOC_TIME_KEY] = '"' + rpt_alloc_time + '"'
#    if db_region is not None:
#        job_setup_document[JOBS_SETUP_DB_REGION_KEY] = db_region
#    if (sample_infos is not None) and isfile(sample_infos):
#        s_stream = file(sample_infos, "r")
#        document = yaml.safe_load(s_stream)
#        job_setup_document[JOBS_SETUP_SAMPLE_INFOS_KEY] = document[JOBS_SETUP_SAMPLE_INFOS_KEY]
#    elif (sample_infos is not None) and (sample_infos.find(":") == -1):
#        job_setup_document[JOBS_SETUP_SAMPLE_INFOS_KEY] = sample_infos
#    elif (sample_infos is not None) and (sample_infos.find(":") > -1):
#        families_info = []
#        for family_item in sample_infos.split(","):
#            family_info = {}
#            family_infos = family_item.split(":")
#            family_info[JOBS_SETUP_FAMILY_ID_KEY] = family_infos[0]
#            members = []
#            for member_infos in family_infos[1:]:
#                member = {}
#                member[JOBS_SETUP_SAMPLE_ID_KEY] = '"' + member_infos + '"'
#                members.append(member)
#            family_info[JOBS_SETUP_MEMBERS_LIST_KEY] = members
#            families_info.append(family_info)
#        job_setup_document[JOBS_SETUP_SAMPLE_INFOS_KEY] = families_info
#    job_setup_document[JOBS_SETUP_OUTPUT_DIR_KEY] = project_out_dir
#    if vcf_tabix_file is not None:
#        job_setup_document[JOBS_SETUP_VCF_TABIX_FILE_KEY] = vcf_tabix_file
#    job_setup_document[JOBS_SETUP_JOBS_REPORT_FILE_KEY] = jobs_report_file
#    job_setup_document[JOBS_SETUP_ANV_SECTION] = annovar_config
#

