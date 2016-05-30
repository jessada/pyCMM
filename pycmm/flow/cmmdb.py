import sys
import re
import datetime
import pyaml
import pysam
import yaml
from os import listdir
from os.path import join as join_path
from os.path import isdir
from os.path import isfile
from collections import defaultdict
from pycmm.settings import DFLT_ANNOVAR_DB_FOLDER
from pycmm.settings import DFLT_ANNOVAR_DB_NAMES
from pycmm.settings import DFLT_ANNOVAR_DB_OPS
from pycmm.settings import DFLT_CMMDB_ALLOC_TIME
from pycmm.settings import DFLT_MUTREP_ALLOC_TIME
from pycmm.settings import DUMMY_TABLE_ANNOVAR_BIN
from pycmm.template import pyCMMBase
from pycmm.utils import exec_sh
from pycmm.utils import mylogger
from pycmm.utils.jobman import JobManager
from pycmm.cmmlib.dnalib import ALL_CHROMS
from pycmm.cmmlib.annovarlib import Annovar
from pycmm.cmmlib.annovarlib import AnnovarParams
from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_INPUT_FILE_KEY
from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_DB_FOLDER_KEY
from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_BUILDVER_KEY
from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_OUT_PREFIX_KEY
from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_DB_NAMES_KEY
from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_DB_OPS_KEY
from pycmm.cmmlib.annovarlib import ANNOVAR_PARAMS_NASTRING_KEY
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.settings import DFLT_MUTREP_FREQ_RATIOS

CAL_MUTATIONS_STAT_SCRIPT = "$PYCMM/bash/cal_mutations_stat.sh"

# *************** jobs metadata section ***************
JOBS_SETUP_DATASET_NAME_KEY = "DATASET_NAME"
JOBS_SETUP_PROJECT_CODE_KEY = "PROJECT_CODE"
JOBS_SETUP_DB_ALLOC_TIME_KEY = "DB_ALLOC_TIME"
JOBS_SETUP_RPT_ALLOC_TIME_KEY = "RPT_ALLOC_TIME"
JOBS_SETUP_OUTPUT_DIR_KEY = "OUTPUT_DIR"
JOBS_SETUP_JOBS_REPORT_FILE_KEY = "JOBS_REPORT_FILE"
JOBS_SETUP_VCF_TABIX_FILE_KEY = "VCF_TABIX_FILE"
JOBS_SETUP_DB_REGION_KEY = "DB_REGION"

# *************** jobs ANNOVAR section ***************
JOBS_SETUP_ANNOVAR_SECTION = "ANNOVAR"
JOBS_SETUP_ANNOVAR_DB_FOLDER_KEY = ANNOVAR_PARAMS_DB_FOLDER_KEY
JOBS_SETUP_ANNOVAR_BUILDVER_KEY = ANNOVAR_PARAMS_BUILDVER_KEY
# >> table_annovar section
JOBS_SETUP_TA_DB_NAMES_KEY = ANNOVAR_PARAMS_DB_NAMES_KEY
JOBS_SETUP_TA_DB_OPS_KEY = ANNOVAR_PARAMS_DB_OPS_KEY
JOBS_SETUP_TA_NASTRING_KEY = ANNOVAR_PARAMS_NASTRING_KEY

# *************** jobs samples section ***************
JOBS_SETUP_SAMPLE_INFOS_KEY = "SAMPLE_INFOS"
JOBS_SETUP_FAMILY_ID_KEY = "FAMILY_ID"
JOBS_SETUP_MEMBERS_LIST_KEY = "MEMBERS"
JOBS_SETUP_SAMPLE_ID_KEY = "SAMPLE_ID"

# *************** report layout section ***************
JOBS_SETUP_RPT_LAYOUT_SECTION = "REPORT_LAYOUT"
JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX = "ANNOTATED_VCF_TABIX"
JOBS_SETUP_RPT_ANNO_COLS_KEY = "COLUMNS"
JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY = "ANNOTATION_EXCLUSION_TAGS"
JOBS_SETUP_RPT_REGIONS_KEY = "REGIONS"
JOBS_SETUP_RPT_FREQ_RATIOS_KEY = "FREQUENCY_RATIOS"
JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY = "COLUMN"
JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY = "FREQUENCY"
JOBS_SETUP_RPT_EXPRESSIONS_KEY = "EXPRESSIONS"
JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY = "NAME"
JOBS_SETUP_RPT_EXPRESSIONS_PATTERN_KEY = "PATTERN"
JOBS_SETUP_RPT_EXPRESSIONS_USAGES_KEY = "USAGES"
JOBS_SETUP_RPT_EXPRESSIONS_ACTION_KEY = "ACTION"
JOBS_SETUP_RPT_EXPRESSIONS_INFO_KEY = "INFO"
JOBS_SETUP_RPT_SPLIT_CHROM_KEY = "SPLIT_CHROM"
JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY = "SUMMARY_FAMILIES"
JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY = "EXTRA_ANNOTATION_COLUMNS"
JOBS_SETUP_RPT_CALL_DETAIL_KEY = "Calling_detail"
JOBS_SETUP_RPT_MT_KEY = "Mitochondria"
JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY = "ROWS_FILTER_ACTIONS_CRITERIA"
JOBS_SETUP_RPT_FILTER_RARE = "Rare"
JOBS_SETUP_RPT_FILTER_NON_INTERGENIC = "Non-Intergenic"
JOBS_SETUP_RPT_FILTER_NON_INTRONIC = "Non-Intronic"
JOBS_SETUP_RPT_FILTER_NON_UPSTREAM = "Non-Upstream"
JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM = "Non-Downstream"
JOBS_SETUP_RPT_FILTER_NON_UTR = "Non-UTR"
JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS = "Non-Synonymous"
JOBS_SETUP_RPT_FILTER_HAS_MUTATION = "Has-Mutation"
JOBS_SETUP_RPT_FILTER_HAS_SHARED = "Has-Shared"
JOBS_SETUP_RPT_ONLY_SUMMARY_KEY = "ONLY_SUMMARY"
JOBS_SETUP_RPT_ONLY_FAMILIES_KEY = "ONLY_FAMILIES"

class MemberInfo(pyCMMBase):
    """  To encapsulate family information so that it is readable """

    def __init__(self, info):
        pyCMMBase.__init__(self)
        self.__info = info

    def get_raw_repr(self):
        return {"sample id": self.sample_id,
                }

    @property
    def sample_id(self):
        return self.__info[JOBS_SETUP_SAMPLE_ID_KEY]

class FamilyInfo(pyCMMBase):
    """  To encapsulate family information so that it is readable """

    def __init__(self, info):
        pyCMMBase.__init__(self)
        self.__info = info
        self.__shared_mutations = None

    def get_raw_repr(self):
        out_txt = {}
        out_txt["family id"] = self.fam_id
        out_txt["members"] = self.members
        if self.shared_mutations is not None:
            out_txt["shared mutations"] = self.shared_mutations
        return out_txt

    @property
    def fam_id(self):
        return str(self.__info[JOBS_SETUP_FAMILY_ID_KEY])

    @property
    def members(self):
        return map(lambda x: MemberInfo(x),
                   self.__info[JOBS_SETUP_MEMBERS_LIST_KEY])

    @property
    def shared_mutations(self):
        return self.__shared_mutations

    @shared_mutations.setter
    def shared_mutations(self, val):
        self.__shared_mutations = val

class CMMDBPipeline(JobManager):
    """ A class to control CMM data(base) workflow """

    def __init__(self,
                 jobs_setup_file,
                 ):
        self.__load_jobs_info(jobs_setup_file)
        self.__jobs_setup_file = jobs_setup_file
        JobManager.__init__(self,
                            jobs_report_file=self.jobs_report_file)
        self.__parse_sample_infos()
        self.__annovar_configs = None
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
    def db_alloc_time(self):
        return self._jobs_info[JOBS_SETUP_DB_ALLOC_TIME_KEY]

    @property
    def input_vcf_tabix(self):
        if JOBS_SETUP_VCF_TABIX_FILE_KEY in self._jobs_info:
            return self._jobs_info[JOBS_SETUP_VCF_TABIX_FILE_KEY]

    @property
    def db_region(self):
        if JOBS_SETUP_DB_REGION_KEY in self._jobs_info:
            return str(self._jobs_info[JOBS_SETUP_DB_REGION_KEY])
        else:
            return None

    @property
    def out_stat_file(self):
        return self.__out_stat_file

    @property
    def samples_list(self):
        return self.__samples_list

    @property
    def samples_list_w_fam_pref(self):
        return self.__samples_list_w_fam_pref

    @property
    def family_infos(self):
        return self.__family_infos

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
        if self.__annovar_configs is None:
            self.__parse_annovar_configs()
        return self.__annovar_configs

    @property
    def report_layout(self):
        return self._jobs_info[JOBS_SETUP_RPT_LAYOUT_SECTION]

    @property
    def jobs_setup_file(self):
        return self.__jobs_setup_file

    def __load_jobs_info(self, jobs_setup_file):
        stream = file(jobs_setup_file, "r")
        self._jobs_info = yaml.safe_load(stream)

    def __parse_sample_infos(self):
        if JOBS_SETUP_SAMPLE_INFOS_KEY not in self._jobs_info:
            self.__samples_list = None
            self.__family_infos = None
            self.__samples_list_w_fam_pref = None
            return
        sample_infos = self._jobs_info[JOBS_SETUP_SAMPLE_INFOS_KEY]
        if type(sample_infos) is str:
            self.__samples_list = sample_infos.split(",")
            self.__family_infos = None
            self.__samples_list_w_fam_pref = None
        else:
            self.__samples_list = []
            self.__family_infos = {}
            self.__samples_list_w_fam_pref = []
            for entry in sample_infos:
                family_info = FamilyInfo(entry)
                for member in family_info.members:
                    sample_id = member.sample_id
                    sample_id_w_pref = family_info.fam_id + '-' + sample_id
                    self.__samples_list.append(sample_id)
                    self.__samples_list_w_fam_pref.append(sample_id_w_pref)
                self.__family_infos[family_info.fam_id] = family_info

    def __parse_annovar_configs(self):
        cfg = self._jobs_info[JOBS_SETUP_ANNOVAR_SECTION]
        cfg[ANNOVAR_PARAMS_OUT_PREFIX_KEY] = join_path(self.data_out_dir,
                                                       self.dataset_name)
        cfg[ANNOVAR_PARAMS_INPUT_FILE_KEY] = self.input_vcf_tabix
        self.__annovar_configs = AnnovarParams(cfg)

    def __create_directories(self):
        self.create_dir(self.rpts_out_dir)
        self.create_dir(self.data_out_dir)
        self.create_dir(self.slurm_log_dir)
        self.create_dir(self.working_dir)

    def __get_cal_mut_stat_params(self,
                                  dataset_name,
                                  out_stat_file,
                                  db_region=None,
                                  samples_list=None,
                                  ):
        params = " -k " + dataset_name
        params += " -i " + self.input_vcf_tabix
        if db_region is not None:
            params += " -r " + db_region
        if samples_list is not None:
            if type(samples_list) is list:
                params += " -c " + ",".join(samples_list)
            else:
                params += " -c " + samples_list
        params += " -o " + out_stat_file
        return params

    def cal_mut_stat(self):
        if self.project_code is None:
            self.__out_stat_file = join_path(self.data_out_dir,
                                             self.dataset_name + ".stat")
            params = self.__get_cal_mut_stat_params(dataset_name=self.dataset_name,
                                                    out_stat_file=self.out_stat_file,
                                                    db_region=self.db_region,
                                                    samples_list=self.samples_list,
                                                    )
            cmd = CAL_MUTATIONS_STAT_SCRIPT + params
            exec_sh(cmd)
        elif self.db_region is None:
            self.__out_stat_file = []
            for chrom in ALL_CHROMS:
                dataset_name = self.dataset_name + "_" + chrom
                out_stat_file = join_path(self.data_out_dir,
                                          dataset_name + ".stat")
                params = self.__get_cal_mut_stat_params(dataset_name=dataset_name,
                                                        out_stat_file=out_stat_file,
                                                        db_region=chrom,
                                                        samples_list=self.samples_list,
                                                        )
                self.__out_stat_file.append(out_stat_file)
                job_name = dataset_name + "_cal_stat"
                slurm_log_file = join_path(self.slurm_log_dir,
                                           job_name)
                slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
                slurm_log_file += ".log"
# *********************************************************************************************** Need refactoring ***********************************************************************************************
                self.submit_job(job_name,
                                self.project_code,
                                "core",
                                "1",
                                self.db_alloc_time,
                                slurm_log_file,
                                CAL_MUTATIONS_STAT_SCRIPT,
                                params,
                                )
# *********************************************************************************************** Need refactoring ***********************************************************************************************
        else:
            self.__out_stat_file = join_path(self.data_out_dir,
                                             self.dataset_name + ".stat")
            params = self.__get_cal_mut_stat_params(dataset_name=self.dataset_name,
                                                    out_stat_file=self.out_stat_file,
                                                    db_region=self.db_region,
                                                    samples_list=self.samples_list,
                                                    )
            job_name = self.dataset_name + "_cal_stat"
            slurm_log_file = join_path(self.slurm_log_dir,
                                       job_name)
            slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
            slurm_log_file += ".log"
# *********************************************************************************************** Need refactoring ***********************************************************************************************
            self.submit_job(job_name,
                            self.project_code,
                            "core",
                            "1",
                            self.db_alloc_time,
                            slurm_log_file,
                            CAL_MUTATIONS_STAT_SCRIPT,
                            params,
                            )
# *********************************************************************************************** Need refactoring ***********************************************************************************************

    def table_annovar(self):
        cfg = self.annovar_config
        table_annovar_cmd = DUMMY_TABLE_ANNOVAR_BIN
        table_annovar_cmd += " --dataset_name " + self.dataset_name
        table_annovar_cmd += " --input_file " + cfg.input_file
        table_annovar_cmd += " --db_folder " + cfg.db_folder
        table_annovar_cmd += " --buildver " + cfg.buildver
        table_annovar_cmd += " --protocols " + cfg.protocols
        table_annovar_cmd += " --operations " + cfg.operations
        table_annovar_cmd += " --nastring " + cfg.nastring
        table_annovar_cmd += " --data_out_folder " + self.output_dir
        if self.project_code is not None:
            job_name = self.dataset_name + "_ta"
            slurm_log_file = join_path(self.slurm_log_dir,
                                       job_name)
            slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
            slurm_log_file += ".log"
# *********************************************************************************************** Need refactoring ***********************************************************************************************
            self.submit_job(job_name,
                            self.project_code,
                            "core",
                            "2",
                            self.db_alloc_time,
                            slurm_log_file,
                            table_annovar_cmd,
                            "",
                            )
# *********************************************************************************************** Need refactoring ***********************************************************************************************
        else:
            exec_sh(table_annovar_cmd)

    def __garbage_collecting(self):
        pass

    def monitor_action(self):
        JobManager.monitor_action(self)
        self.__garbage_collecting()

def create_jobs_setup_file(dataset_name,
                           project_out_dir,
                           vcf_tabix_file=None,
                           db_region=None,
                           sample_infos=None,
                           project_code=None,
                           db_alloc_time=DFLT_CMMDB_ALLOC_TIME,
                           rpt_alloc_time=DFLT_MUTREP_ALLOC_TIME,
                           jobs_report_file=None,
                           annovar_human_db_dir=DFLT_ANNOVAR_DB_FOLDER,
                           annovar_buildver="hg19",
                           annovar_db_names=DFLT_ANNOVAR_DB_NAMES,
                           annovar_db_ops=DFLT_ANNOVAR_DB_OPS,
                           annovar_nastring=".",
                           anno_cols=None,
                           anno_excl_tags=None,
                           rows_filter_actions=None,
                           annotated_vcf_tabix=None,
                           report_regions=None,
                           frequency_ratios=None,
                           expression_patterns=None,
                           expression_usages=None,
                           split_chrom=False,
                           summary_families_sheet=False,
                           call_detail=False,
                           only_summary=False,
                           only_families=False,
                           out_jobs_setup_file=None,
                           ):
    if jobs_report_file is None:
        jobs_report_file = join_path(project_out_dir,
                                     dataset_name+"_rpt.txt")
    if out_jobs_setup_file is None:
        out_jobs_setup_file = join_path(project_out_dir,
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
    if db_alloc_time is None:
        db_alloc_time = DFLT_CMMDB_ALLOC_TIME
    job_setup_document[JOBS_SETUP_DB_ALLOC_TIME_KEY] = '"' + db_alloc_time + '"'
    if rpt_alloc_time is None:
        rpt_alloc_time = DFLT_MUTREP_ALLOC_TIME
    job_setup_document[JOBS_SETUP_RPT_ALLOC_TIME_KEY] = '"' + rpt_alloc_time + '"'
    if db_region is not None:
        job_setup_document[JOBS_SETUP_DB_REGION_KEY] = db_region
    if (sample_infos is not None) and isfile(sample_infos):
        s_stream = file(sample_infos, "r")
        document = yaml.safe_load(s_stream)
        job_setup_document[JOBS_SETUP_SAMPLE_INFOS_KEY] = document[JOBS_SETUP_SAMPLE_INFOS_KEY]
    elif (sample_infos is not None) and (sample_infos.find(":") == -1):
        job_setup_document[JOBS_SETUP_SAMPLE_INFOS_KEY] = sample_infos
    elif (sample_infos is not None) and (sample_infos.find(":") > -1):
        families_info = []
        for family_item in sample_infos.split(","):
            family_info = {}
            family_infos = family_item.split(":")
            family_info[JOBS_SETUP_FAMILY_ID_KEY] = family_infos[0]
            members = []
            for member_infos in family_infos[1:]:
                member = {}
                member[JOBS_SETUP_SAMPLE_ID_KEY] = '"' + member_infos + '"'
                members.append(member)
            family_info[JOBS_SETUP_MEMBERS_LIST_KEY] = members
            families_info.append(family_info)
        job_setup_document[JOBS_SETUP_SAMPLE_INFOS_KEY] = families_info
    job_setup_document[JOBS_SETUP_OUTPUT_DIR_KEY] = project_out_dir
    if vcf_tabix_file is not None:
        job_setup_document[JOBS_SETUP_VCF_TABIX_FILE_KEY] = vcf_tabix_file
    job_setup_document[JOBS_SETUP_JOBS_REPORT_FILE_KEY] = jobs_report_file
    job_setup_document[JOBS_SETUP_ANNOVAR_SECTION] = annovar_config
    report_layout_config = {}
    if anno_cols is None:
        report_layout_config[JOBS_SETUP_RPT_ANNO_COLS_KEY] = ALL_MUTREP_ANNO_COLS.keys()
    else:
        report_layout_config[JOBS_SETUP_RPT_ANNO_COLS_KEY] = anno_cols.split(",")
    if anno_excl_tags is not None:
        report_layout_config[JOBS_SETUP_RPT_ANNO_EXCL_TAGS_KEY] = anno_excl_tags.split(",")
    if annotated_vcf_tabix is not None:
        report_layout_config[JOBS_SETUP_RPT_ANNOTATED_VCF_TABIX] = annotated_vcf_tabix
    if report_regions is not None:
        report_layout_config[JOBS_SETUP_RPT_REGIONS_KEY] = report_regions.split(",")
    if frequency_ratios is None:
        frequency_ratios = DFLT_MUTREP_FREQ_RATIOS
    job_freq_ratios = []
    for frequency_ratio in frequency_ratios.split(","):
        (col, freq) = frequency_ratio.split(":")
        job_freq_ratios.append({JOBS_SETUP_RPT_FREQ_RATIOS_COL_KEY: col,
                                JOBS_SETUP_RPT_FREQ_RATIOS_FREQ_KEY: freq,
                                })
    report_layout_config[JOBS_SETUP_RPT_FREQ_RATIOS_KEY] = job_freq_ratios
    if expression_patterns is not None:
        exprs = []
        for raw_pattern in  expression_patterns.split(";"):
            name, pattern = raw_pattern.split(":")
            expr = defaultdict(list)
            expr[JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY] = name.strip()
            expr[JOBS_SETUP_RPT_EXPRESSIONS_PATTERN_KEY] = pattern.strip()
            exprs.append(expr)
        if expression_usages is not None:
            for raw_usage in expression_usages.split(";"):
                usage = {}
                usage_expr_name = raw_usage.split(":")[0]
                usage[JOBS_SETUP_RPT_EXPRESSIONS_ACTION_KEY] = raw_usage.split(":")[1]
                if len(raw_usage.split(":")) > 2:
                    usage[JOBS_SETUP_RPT_EXPRESSIONS_INFO_KEY] = ":".join(raw_usage.split(":")[2:])
                for expr in exprs:
                    if expr[JOBS_SETUP_RPT_EXPRESSIONS_NAME_KEY] == usage_expr_name:
                        expr[JOBS_SETUP_RPT_EXPRESSIONS_USAGES_KEY].append(usage)
        report_layout_config[JOBS_SETUP_RPT_EXPRESSIONS_KEY] = exprs
    report_layout_config[JOBS_SETUP_RPT_SPLIT_CHROM_KEY] = split_chrom
    report_layout_config[JOBS_SETUP_RPT_SUMMARY_FAMILIES_KEY] = summary_families_sheet
    if call_detail:
        extra_anno_cols = []
        if call_detail:
            extra_anno_cols.append(JOBS_SETUP_RPT_CALL_DETAIL_KEY)
        report_layout_config[JOBS_SETUP_RPT_EXTRA_ANNO_COLS_KEY] = extra_anno_cols
    if rows_filter_actions is not None:
        filter_criterias = []
        for filter_criteria in rows_filter_actions.split(","):
            if filter_criteria == JOBS_SETUP_RPT_FILTER_RARE:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_RARE)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_INTERGENIC:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_INTERGENIC)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_INTRONIC:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_INTRONIC)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_UPSTREAM:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_UPSTREAM)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_DOWNSTREAM)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_UTR:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_UTR)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_NON_SYNONYMOUS)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_HAS_MUTATION:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_HAS_MUTATION)
            if filter_criteria == JOBS_SETUP_RPT_FILTER_HAS_SHARED:
                filter_criterias.append(JOBS_SETUP_RPT_FILTER_HAS_SHARED)
        report_layout_config[JOBS_SETUP_RPT_ROWS_FILTER_ACTIONS_CRITERIA_KEY] = filter_criterias
    report_layout_config[JOBS_SETUP_RPT_ONLY_SUMMARY_KEY] = only_summary
    report_layout_config[JOBS_SETUP_RPT_ONLY_FAMILIES_KEY] = only_families
    job_setup_document[JOBS_SETUP_RPT_LAYOUT_SECTION] = report_layout_config


    pyaml.dump(job_setup_document, stream)
