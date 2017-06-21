import yaml
import pyaml
from os.path import isfile
from collections import OrderedDict
from pycmm.settings import DBMS_EXECUTE_DB_JOBS_BIN
from pycmm.template import pyCMMBase
from pycmm.utils.ver import VersionManager
from pycmm.cmmlib import CMMParams
from pycmm.proc import CMMPipeline
from pycmm.proc import init_jobs_setup_file
from pycmm.proc import get_func_arg
from pycmm.proc.db.connector import SQLiteDB
from pycmm.proc.db.connector import DATA_TYPE_AVDB_INFO
from pycmm.proc.db.connector import DATA_TYPE_ANNOVAR_INFO
from pycmm.proc.db.connector import DATA_TYPE_GTZ
from pycmm.proc.db.dbinput import AVDBReader
from pycmm.proc.db.dbinput import TAVcfInfoReader as VcfInfoReader
from pycmm.proc.db.dbinput import TAVcfGTZReader as VcfGTZReader

## *************** DB params section ***************
JOBS_SETUP_DB_PARAMS_SECTION = "DB_PARAMS"
JOBS_SETUP_DB_FILE_KEY = "DB_FIlE"
JOBS_SETUP_DB_JOBS = "DB_JOBS"
JOBS_SETUP_DB_JOB_TBL_NAME_KEY = "Table_name"
JOBS_SETUP_DB_JOB_DROP_TABLE_KEY = "Drop_table"
JOBS_SETUP_DB_JOB_DATA_FILE_KEY = "Data_file"
JOBS_SETUP_DB_JOB_OPERATION_KEY = "Operation"
JOBS_SETUP_DB_JOB_HEADER_EXIST_KEY = "Header_exist"
## *************** table_annovar db section ***************

OPERATION_LOAD_AVDB = 'LOAD_AVDB'
OPERATION_LOAD_ANNOVAR_VCF = 'LOAD_ANNOVAR_VCF'
OPERATION_LOAD_GTZ_VCF = 'LOAD_GTZ_VCF'


class SQLiteDBWriter(SQLiteDB):
    """
    A class for loading and updating data to CMM SQLite database
    """

    def __init__(self, *args, **kwargs):
        super(SQLiteDBWriter, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(SQLiteDB, self).get_raw_obj_str(*args, **kwargs)
        return raw_repr

    def __load_data_file(self,
                         reader,
                         tbl_name,
                         ):
        rec_len = len(reader.header_cols)
        sql = 'INSERT OR IGNORE INTO ' + tbl_name
        sql += ' VALUES (' + ('?,'*(rec_len*2))[:rec_len*2-1] + ')'
        row_count = self._exec_many_sql(sql, reader).rowcount
        return row_count

    def load_avdb_data(self, data_file, header_exist, *args, **kwargs):
        reader = AVDBReader(file_name=data_file, header_exist=header_exist)
        tbl_name = kwargs['tbl_name']
        self._create_table(tbl_name=tbl_name,
                            col_names=reader.header_cols,
                            pkeys=reader.pkeys,
                            )
        self._update_sys_info(data_type=DATA_TYPE_AVDB_INFO,
                              updating_info_cols=reader.avdb_cols,
                              updating_tbl_name=tbl_name,
                              )
        return self.__load_data_file(reader=reader, *args, **kwargs)

    def load_annovar_vcf_data(self, data_file, *args, **kwargs):
        reader = VcfInfoReader(file_name=data_file)
        tbl_name = kwargs['tbl_name']
        self._create_table(tbl_name=tbl_name,
                            col_names=reader.header_cols,
                            pkeys=reader.pkeys,
                            )
        self._update_sys_info(data_type=DATA_TYPE_ANNOVAR_INFO,
                              updating_info_cols=reader.info_cols,
                              updating_tbl_name=tbl_name,
                              )
        return self.__load_data_file(reader=reader, *args, **kwargs)

    def load_gtz_vcf_data(self, data_file, *args, **kwargs):
        reader = VcfGTZReader(file_name=data_file)
        tbl_name = kwargs['tbl_name']
        self._create_table(tbl_name=tbl_name,
                            col_names=reader.header_cols,
                            pkeys=reader.pkeys,
                            )
        self._update_sys_info(data_type=DATA_TYPE_GTZ,
                              updating_info_cols=reader.samples_id,
                              updating_tbl_name=tbl_name,
                              )
        return self.__load_data_file(reader=reader, *args, **kwargs)

    def __add_column(self, tbl_name, col_name, col_def=""):
        sql = "ALTER TABLE " + tbl_name
        sql += " ADD " + col_name
        sql += " " + col_def
        self._exec_sql(sql)

    def cal_hw(self, tbl_name):
        # look for dataset name
        col_names = self.get_col_names(tbl_name)
        af_col_names = filter(lambda x: x.endswith("AF"),
                              col_names)
        if len(af_col_names) != 1:
            raise Exception('cannot perform Hardy-Weinberg calculation for table ' + tbl_name)
        dataset_name = af_col_names[0].replace("_AF","")
        hw_col_name = dataset_name.strip("OAF_") + "_hw_chisq"
        # perform calculation
        self.__add_column(tbl_name, "p")
        self.__add_column(tbl_name, "q")
        self.__add_column(tbl_name, "exp_wt")
        self.__add_column(tbl_name, "exp_het")
        self.__add_column(tbl_name, "exp_hom")
        self.__add_column(tbl_name, hw_col_name)
        sql = "UPDATE " + tbl_name
        sql += " SET p = 1-" + dataset_name + "_AF"
        sql += ", q = " + dataset_name + "_AF"
        self._exec_sql(sql)
        sql = "UPDATE " + tbl_name
        sql += " SET exp_wt = p*p*" + dataset_name + "_GT"
        sql += ", exp_het = 2*p*q*" + dataset_name + "_GT"
        sql += ", exp_hom = q*q*" + dataset_name + "_GT"
        self._exec_sql(sql)
        sql = "UPDATE " + tbl_name
        sql += " SET " + hw_col_name + " = 0.0"
        sql += " WHERE exp_wt = 0"
        sql += " OR exp_het = 0"
        sql += " OR exp_hom = 0"
        self._exec_sql(sql)
        sql = "UPDATE " + tbl_name
        sql += " SET " + hw_col_name + " = (" + dataset_name + "_WT-exp_wt)*(" + dataset_name + "_WT-exp_wt)/exp_wt"
        sql += " + (" + dataset_name + "_HET-exp_het)*(" + dataset_name + "_HET-exp_het)/exp_het"
        sql += " + (" + dataset_name + "_HOM-exp_hOM)*(" + dataset_name + "_HOM-exp_hom)/exp_hom"
        sql += " WHERE exp_wt != 0"
        sql += " AND exp_het != 0"
        sql += " AND exp_hom != 0"
        self._exec_sql(sql)
        # update system information
        avdb_info = self.get_avdb_info()
        info_cols = avdb_info[tbl_name]
        if hw_col_name not in info_cols:
            info_cols.append(hw_col_name)
        self._update_sys_info(DATA_TYPE_AVDB_INFO, info_cols, tbl_name)

class DBJob(CMMParams):
    """ A class to keep DB job information """

    def __init__(self, *args, **kwargs):
        super(DBJob, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_str = super(DBJob, self).get_raw_obj_str(*args, **kwargs)
        raw_str["table name"] = self.tbl_name
        raw_str["drop table"] = self.drop_table
        raw_str["data file"] = self.data_file
        raw_str["operation"] = self.operation
        raw_str["header exist"] = self.header_exist
        return raw_str

    @property
    def tbl_name(self):
        return self._get_job_config(JOBS_SETUP_DB_JOB_TBL_NAME_KEY)

    @property
    def drop_table(self):
        return self._get_job_config(JOBS_SETUP_DB_JOB_DROP_TABLE_KEY) is True

    @property
    def data_file(self):
        return self._get_job_config(JOBS_SETUP_DB_JOB_DATA_FILE_KEY)

    @property
    def operation(self):
        return self._get_job_config(JOBS_SETUP_DB_JOB_OPERATION_KEY,
                                    required=True)

    @property
    def header_exist(self):
        return self._get_job_config(JOBS_SETUP_DB_JOB_HEADER_EXIST_KEY)

class DBParams(CMMParams):
    """ A class to keep all the job params related to db controller """

    def __init__(self, *args, **kwargs):
        super(DBParams, self).__init__(*args, **kwargs)
        self.__db_jobs = self.__parse_db_jobs(self._get_job_config(JOBS_SETUP_DB_JOBS))

    def get_raw_obj_str(self, *args, **kwargs):
        raw_str = super(DBParams, self).get_raw_obj_str(*args, **kwargs)
        raw_str["db file"] = self.db_file
        raw_str["db jobs"] = self.db_jobs
        return raw_str

    def __parse_db_jobs(self, db_jobs_params):
        return map(lambda x: DBJob(x), db_jobs_params)

    @property
    def db_file(self):
        return self._get_job_config(JOBS_SETUP_DB_FILE_KEY,
                                    required=True)

    @property
    def db_jobs(self):
        return self.__db_jobs

class SQLiteDBController(CMMPipeline):
    """ A class to control all writing processes related to SQLiteDB """

    def __init__(self, jobs_setup_file, verbose=True, *args, **kwargs):
        kwargs['jobs_setup_file'] = jobs_setup_file
        super(SQLiteDBController, self).__init__(*args, **kwargs)
        self.__db_params = DBParams(self._get_job_config(JOBS_SETUP_DB_PARAMS_SECTION,
                                                         required=True))
        self.__db = SQLiteDBWriter(self.db_params.db_file, verbose=verbose)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_str = super(SQLiteDBController, self).get_raw_obj_str(*args, **kwargs)
        return raw_str

    def get_third_party_software_version(self):
        vm = VersionManager()
        versions = OrderedDict()
        return versions

    @property
    def db(self):
        return self.__db

    @property
    def db_params(self):
        return self.__db_params

    def __execute_db_job(self, db_job):
        info_msg = "executing db job ("
        info_msg += " table_name = {tbl_name}"
        info_msg += ", drop_table = {drop_table}"
        info_msg += ", data_file = {data_file}"
        info_msg += ", operation = {operation}"
        info_msg += ", header_exist = {header_exist}"
        info_msg += " )"
        self.info()
        self.info(info_msg.format(tbl_name=db_job.tbl_name,
                                  drop_table=db_job.drop_table,
                                  data_file=db_job.data_file,
                                  operation=db_job.operation,
                                  header_exist=db_job.header_exist,
                                  ))
        if db_job.drop_table:
            self.db.drop_table(db_job.tbl_name)
        if db_job.operation == OPERATION_LOAD_AVDB:
            return self.db.load_avdb_data(data_file=db_job.data_file,
                                          tbl_name=db_job.tbl_name,
                                          header_exist=db_job.header_exist,
                                          )
        elif db_job.operation == OPERATION_LOAD_ANNOVAR_VCF:
            return self.db.load_annovar_vcf_data(data_file=db_job.data_file,
                                                 tbl_name=db_job.tbl_name,
                                                 )
        elif db_job.operation == OPERATION_LOAD_GTZ_VCF:
            return self.db.load_gtz_vcf_data(data_file=db_job.data_file,
                                             tbl_name=db_job.tbl_name,
                                             )
        else:
            raise Exception('Unknown operation: ' + db_job.operation)

    def __execute_db_jobs(self):
        for db_job in self.db_params.db_jobs:
            n_records = self.__execute_db_job(db_job)
            self.info("Finished: " + str(n_records) + " records were inserted")

    def run_offline_pipeline(self):
        self.__execute_db_jobs()

    def run_controller(self):
        if self.project_code is None:
            self.run_offline_pipeline()
        else:
            job_name = self.dataset_name + "_exec_db"
            params = " -j " + self.jobs_setup_file
            self._submit_slurm_job(job_name,
                                   "1",
                                   DBMS_EXECUTE_DB_JOBS_BIN,
                                   params,
                                   )


def parse_db_jobs(db_jobs_params):
    db_jobs_yaml = []
    if db_jobs_params is None or db_jobs_params == "":
        return db_jobs_yaml
    if isfile(db_jobs_params):
        s_stream = file(db_jobs_params, "r")
        document = yaml.safe_load(s_stream)
        return document[JOBS_SETUP_DB_JOBS]
    # so the params are in text
    db_jobs = db_jobs_params.split(",")
    for db_job in db_jobs:
        job_info = {}
        job_details = db_job.split(":")
        # tbl_name:drop_table:data_file:opertion
        job_info[JOBS_SETUP_DB_JOB_TBL_NAME_KEY] = job_details[0]
        if len(job_details[1]) > 0:
            job_info[JOBS_SETUP_DB_JOB_DROP_TABLE_KEY] = job_details[1]
        job_info[JOBS_SETUP_DB_JOB_DATA_FILE_KEY] = job_details[2]
        job_info[JOBS_SETUP_DB_JOB_OPERATION_KEY] = job_details[3]
        if len(job_details[4]) > 0:
            job_info[JOBS_SETUP_DB_JOB_HEADER_EXIST_KEY] = job_details[4]
        db_jobs_yaml.append(job_info)
    return db_jobs_yaml

def create_dbms_jobs_setup_file(*args, **kwargs):
    job_setup_document, stream = init_jobs_setup_file(*args, **kwargs)

    db_job_infos = {}
    db_job_infos[JOBS_SETUP_DB_FILE_KEY] = get_func_arg('db_file', kwargs)
    db_job_infos[JOBS_SETUP_DB_JOBS] = parse_db_jobs(get_func_arg('db_jobs', kwargs))

    job_setup_document[JOBS_SETUP_DB_PARAMS_SECTION] = db_job_infos
    pyaml.dump(job_setup_document, stream)
