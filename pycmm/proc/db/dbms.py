import yaml
import pyaml
import sqlite3
from os.path import isfile
from collections import OrderedDict
from pycmm.settings import DBMS_EXECUTE_DB_JOBS_BIN
from pycmm.template import pyCMMBase
from pycmm.utils.ver import VersionManager
from pycmm.cmmlib import CMMParams
from pycmm.proc import CMMPipeline
from pycmm.proc import init_jobs_setup_file
from pycmm.proc import get_func_arg
from pycmm.proc.db.dbinput import AVDBReader
from pycmm.proc.db.dbinput import TAVcfInfoReader as VcfInfoReader
from pycmm.proc.db.dbinput import TAVcfGTZReader as VcfGTZReader

## *************** DB params section ***************
JOBS_SETUP_DB_PARAMS_SECTION = "DB_PARAMS"
JOBS_SETUP_DB_FILE_KEY = "DB_FIlE"
JOBS_SETUP_DB_JOBS = "DB_JOBS"
JOBS_SETUP_DB_JOB_TABLE_NAME_KEY = "Table_name"
JOBS_SETUP_DB_JOB_DROP_TABLE_KEY = "Drop_table"
JOBS_SETUP_DB_JOB_DATA_FILE_KEY = "Data_file"
JOBS_SETUP_DB_JOB_OPERATION_KEY = "Operation"
JOBS_SETUP_DB_JOB_HEADER_EXIST_KEY = "Header_exist"
## *************** table_annovar db section ***************

OPERATION_LOAD_AVDB = 'LOAD_AVDB'
OPERATION_LOAD_ANNOVAR_VCF = 'LOAD_ANNOVAR_VCF'
OPERATION_LOAD_GTZ_VCF = 'LOAD_GTZ_VCF'


class SQLiteDB(pyCMMBase):
    """ A class to provide APIs for conntecting to CMM custrom SQLite DB  """

    def __init__(self, db_path, **kwargs):
        super(SQLiteDB, self).__init__(**kwargs)
        self.__db_path = db_path
        self.info("connecting to database: " + db_path)

    def get_raw_obj_str(self, **kwargs):
        raw_repr = super(SQLiteDB, self).get_raw_obj_str(**kwargs)
        return raw_repr

    @property
    def db_path(self):
        return self.__db_path

    def __exec_sql(self, sql):
        # to encapsulate sqlite sql execution
        # so far just for logging
        conn = self.__connect_db() 
        self.dbg("executing sql: " + sql)
        result = conn.execute(sql)
        self.__close_connection()
        return result

    def __exec_many_sql(self, sql, records):
        # to encapsulate sqlite sql execution
        # so far just for logging
        conn = self.__connect_db() 
        self.dbg("executing many sql: " + sql)
        result = conn.executemany(sql, records)
        self.__close_connection()
        return result

    def __connect_db(self):
        self.__conn = sqlite3.connect(self.db_path)
        return self.__conn
#        return self.__conn.cursor()

    def __close_connection(self):
        self.__conn.commit()
        self.__conn.close()

    def __load_data_file(self,
                         data_file,
                         reader,
                         table_name,
                         header_exist=None,
                         drop_table=True,
                         ):
        if header_exist is not None:
            rd = reader(file_name=data_file, header_exist=header_exist)
        else:
            rd = reader(file_name=data_file)
        rec_len = len(rd.header_cols)
        if drop_table:
            sql = 'DROP TABLE if exists ' + table_name
            self.__exec_sql(sql)
        sql = 'CREATE TABLE if not exists ' + table_name
        sql += ' (' + ",".join(rd.header_cols) + ')' 
        self.__exec_sql(sql)
        sql = 'INSERT INTO ' + table_name
        sql += ' VALUES (' + ('?,'*(rec_len*2))[:rec_len*2-1] + ')'
        row_count = self.__exec_many_sql(sql, rd).rowcount
        return row_count

    def load_avdb_data(self, *args, **kwargs):
        kwargs['reader'] = AVDBReader
        return self.__load_data_file(*args, **kwargs)

    def load_annovar_vcf_data(self, *args, **kwargs):
        kwargs['table_name'] = 'annovar_vcf'
        kwargs['reader'] = VcfInfoReader
        return self.__load_data_file(*args, **kwargs)

    def load_gtz_vcf_data(self, *args, **kwargs):
        kwargs['reader'] = VcfGTZReader
        return self.__load_data_file(*args, **kwargs)

    def __add_column(self, table_name, col_name, col_def=""):
        sql = "ALTER TABLE " + table_name
        sql += " ADD " + col_name
        sql += " " + col_def
        self.__exec_sql(sql)

    def cal_hw(self, table_name):
        # look for dataset name
        cursor = self.__exec_sql('select * from ' + table_name + " LIMIT 1")
        col_names = map(lambda x: x[0], cursor.description)
        af_col_names = filter(lambda x: x.endswith("AF"),
                              col_names)
        if len(af_col_names) != 1:
            raise Exception('cannot perform Hardy-Weinberg calculation for table ' + table_name)
        dataset_name = af_col_names[0].replace("_AF","")
        self.__add_column(table_name, "p")
        self.__add_column(table_name, "q")
        self.__add_column(table_name, "exp_wt")
        self.__add_column(table_name, "exp_het")
        self.__add_column(table_name, "exp_hom")
        self.__add_column(table_name, "hw_chisq")
        sql = "update " + table_name
        sql += " set p = 1-" + dataset_name + "_AF"
        sql += ", q = " + dataset_name + "_AF"
        self.__exec_sql(sql)
        sql = "update " + table_name
        sql += " set exp_wt = p*p*" + dataset_name + "_GT"
        sql += ", exp_het = 2*p*q*" + dataset_name + "_GT"
        sql += ", exp_hom = q*q*" + dataset_name + "_GT"
        self.__exec_sql(sql)
        sql = "update " + table_name
        sql += " set hw_chisq = 0.0"
        sql += " where exp_wt = 0"
        sql += " or exp_het = 0"
        sql += " or exp_hom = 0"
        self.__exec_sql(sql)
        sql = "update " + table_name
        sql += " set hw_chisq = (" + dataset_name + "_WT-exp_wt)*(" + dataset_name + "_WT-exp_wt)/exp_wt"
        sql += " + (" + dataset_name + "_HET-exp_het)*(" + dataset_name + "_HET-exp_het)/exp_het"
        sql += " + (" + dataset_name + "_HOM-exp_hOM)*(" + dataset_name + "_HOM-exp_hom)/exp_hom"
        sql += " where exp_wt != 0"
        sql += " and exp_het != 0"
        sql += " and exp_hom != 0"
        self.__exec_sql(sql)

    def table_exist(self, table_name):
        c = self.__connect_db()
        sql = "SELECT name FROM sqlite_master WHERE type='table'"
        sql += "AND name='" + table_name + "'"
        tbl_exist = c.execute(sql).fetchone()
        self.__close_connection() 
        if tbl_exist:
            return True
        return False

class DBJob(CMMParams):
    """ A class to keep DB job information """

    def __init__(self, *args, **kwargs):
        super(DBJob, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_str = super(DBJob, self).get_raw_obj_str(*args, **kwargs)
        raw_str["table name"] = self.table_name
        raw_str["drop table"] = self.drop_table
        raw_str["data file"] = self.data_file
        raw_str["operation"] = self.operation
        raw_str["header exist"] = self.header_exist
        return raw_str

    @property
    def table_name(self):
        return self._get_job_config(JOBS_SETUP_DB_JOB_TABLE_NAME_KEY)

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
    """ A class to control all processes related to SQLiteDB """

    def __init__(self, *args, **kwargs):
        super(SQLiteDBController, self).__init__(*args, **kwargs)
        self.__db_params = DBParams(self._get_job_config(JOBS_SETUP_DB_PARAMS_SECTION,
                                                         required=True))
        self.__db = SQLiteDB(self.db_params.db_file)

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
        info_msg += " table_name = {table_name}"
        info_msg += ", drop_table = {drop_table}"
        info_msg += ", data_file = {data_file}"
        info_msg += ", operation = {operation}"
        info_msg += ", header_exist = {header_exist}"
        info_msg += " )"
        self.info()
        self.info(info_msg.format(table_name=db_job.table_name,
                                  drop_table=db_job.drop_table,
                                  data_file=db_job.data_file,
                                  operation=db_job.operation,
                                  header_exist=db_job.header_exist,
                                  ))
        if db_job.operation == OPERATION_LOAD_AVDB:
            return self.db.load_avdb_data(data_file=db_job.data_file,
                                          table_name=db_job.table_name,
                                          header_exist=db_job.header_exist,
                                          drop_table=db_job.drop_table,
                                          )
        elif db_job.operation == OPERATION_LOAD_ANNOVAR_VCF:
            return self.db.load_annovar_vcf_data(data_file=db_job.data_file,
                                                 drop_table=db_job.drop_table,
                                                 )
        elif db_job.operation == OPERATION_LOAD_GTZ_VCF:
            return self.db.load_gtz_vcf_data(data_file=db_job.data_file,
                                             table_name=db_job.table_name,
                                             drop_table=db_job.drop_table,
                                             )
        else:
            raise Exception('Unknown operation: ' + db_job.operation)

    def __execute_db_jobs(self):
        for db_job in self.db_params.db_jobs:
            n_records = self.__execute_db_job(db_job)
            self.info("Finished: " + str(n_records) + " records are inserted")

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
        job_info[JOBS_SETUP_DB_JOB_TABLE_NAME_KEY] = job_details[0]
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
