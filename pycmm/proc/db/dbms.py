import yaml
import pyaml
import sqlite3
from os.path import isfile
from collections import OrderedDict
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
JOBS_SETUP_DB_JOB_DATA_TYPE_KEY = "Data_type"
JOBS_SETUP_DB_JOB_HEADER_EXIST_KEY = "Header_exist"
## *************** table_annovar db section ***************

DATA_TYPE_AVDB = 'AVDB'
DATA_TYPE_ANNOVAR_VCF = 'ANNOVAR_VCF'
DATA_TYPE_GTZ_VCF = 'GTZ_VCF'


class SQLiteDB(pyCMMBase):
    """ A class to provide APIs for conntecting to CMM custrom SQLite DB  """

    def __init__(self, db_path, **kwargs):
        super(SQLiteDB, self).__init__(**kwargs)
        self.__db_path = db_path

    def get_raw_obj_str(self, **kwargs):
        raw_repr = super(SQLiteDB, self).get_raw_obj_str(**kwargs)
        return raw_repr

    @property
    def db_path(self):
        return self.__db_path

    def table_exist(self, table_name):
        c = self.__connect_db()
        qry = "SELECT name FROM sqlite_master WHERE type='table'"
        qry += "AND name='" + table_name + "'"
        tbl_exist = c.execute(qry).fetchone()
        self.__close_connection() 
        if tbl_exist:
            return True
        return False

    def __connect_db(self):
        self.__conn = sqlite3.connect(self.db_path)
        return self.__conn.cursor()

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
        c = self.__connect_db()
        if drop_table:
            sql_script = 'DROP TABLE if exists ' + table_name
            c.execute(sql_script)
        sql_script = 'CREATE TABLE if not exists ' + table_name
        sql_script += ' (' + ",".join(rd.header_cols) + ')' 
        c.execute(sql_script)
        sql_script = 'INSERT INTO ' + table_name
        sql_script += ' VALUES (' + ('?,'*(rec_len*2))[:rec_len*2-1] + ')'
        c.executemany(sql_script, rd)
        row_count = c.rowcount
        self.__close_connection() 
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

class DBJob(CMMParams):
    """ A class to keep DB job information """

    def __init__(self, *args, **kwargs):
        super(DBJob, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_str = super(DBJob, self).get_raw_obj_str(*args, **kwargs)
        raw_str["table name"] = self.table_name
        raw_str["drop table"] = self.drop_table
        raw_str["data file"] = self.data_file
        raw_str["data type"] = self.data_type
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
        return self._get_job_config(JOBS_SETUP_DB_JOB_DATA_FILE_KEY,
                                    required=True)

    @property
    def data_type(self):
        return self._get_job_config(JOBS_SETUP_DB_JOB_DATA_TYPE_KEY,
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
#        versions['SQLite'] = vm.sqlite_version
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
        info_msg += ", data_type = {data_type}"
        info_msg += ", header_exist = {header_exist}"
        info_msg += " )"
        self.info()
        self.info(info_msg.format(table_name=db_job.table_name,
                                  drop_table=db_job.drop_table,
                                  data_file=db_job.data_file,
                                  data_type=db_job.data_type,
                                  header_exist=db_job.header_exist,
                                  ))
        if db_job.data_type == DATA_TYPE_AVDB:
            return self.db.load_avdb_data(data_file=db_job.data_file,
                                          table_name=db_job.table_name,
                                          header_exist=db_job.header_exist,
                                          drop_table=db_job.drop_table,
                                          )
        elif db_job.data_type == DATA_TYPE_ANNOVAR_VCF:
            return self.db.load_annovar_vcf_data(data_file=db_job.data_file,
                                                 drop_table=db_job.drop_table,
                                                 )
        elif db_job.data_type == DATA_TYPE_GTZ_VCF:
            return self.db.load_gtz_vcf_data(data_file=db_job.data_file,
                                             table_name=db_job.table_name,
                                             drop_table=db_job.drop_table,
                                             )
        else:
            raise Exception('Unknown data type: ' + db_job.data_type)

    def __execute_db_jobs(self):
        for db_job in self.db_params.db_jobs:
            n_records = self.__execute_db_job(db_job)
            self.info("Finished: " + str(n_records) + " records are inserted")

    def run_offline_pipeline(self):
        self.__execute_db_jobs()

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
        job_info[JOBS_SETUP_DB_JOB_DATA_TYPE_KEY] = job_details[3]
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
