import yaml
import pyaml
import sqlite3
from os.path import isfile
from collections import OrderedDict
from collections import defaultdict
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
JOBS_SETUP_DB_JOB_TBL_NAME_KEY = "Table_name"
JOBS_SETUP_DB_JOB_DROP_TABLE_KEY = "Drop_table"
JOBS_SETUP_DB_JOB_DATA_FILE_KEY = "Data_file"
JOBS_SETUP_DB_JOB_OPERATION_KEY = "Operation"
JOBS_SETUP_DB_JOB_HEADER_EXIST_KEY = "Header_exist"
## *************** table_annovar db section ***************

OPERATION_LOAD_AVDB = 'LOAD_AVDB'
OPERATION_LOAD_ANNOVAR_VCF = 'LOAD_ANNOVAR_VCF'
OPERATION_LOAD_GTZ_VCF = 'LOAD_GTZ_VCF'

TBL_NAME_AVDB_INFO = "avdb_info"
TBL_AVDBS_COL_NAME_COL_NAME = "column_name"
TBL_AVDBS_TBL_NAME_COL_NAME = "tbl_name"

TBL_NAME_ANNOVAR_INFO = "annovar_info"
TBL_ANNOVAR_COL_NAME_COL_NAME = "column_name"
TBL_ANNOVAR_TBL_NAME_COL_NAME = "tbl_name"

TBL_NAME_SAMPLES = "samples"
TBL_SAMPLES_SAMPLE_ID_COL_NAME = "sample_id"
TBL_SAMPLES_GTZ_TBL_NAME_COL_NAME = "gtz_tbl_name"


class SQLiteDB(pyCMMBase):
    """ A class to provide APIs for conntecting to CMM custrom SQLite DB  """

    def __init__(self, db_path, verbose=True, **kwargs):
        super(SQLiteDB, self).__init__(**kwargs)
        self.__db_path = db_path
        self.__verbose = verbose
        self.__init_db()

    def get_raw_obj_str(self, **kwargs):
        raw_repr = super(SQLiteDB, self).get_raw_obj_str(**kwargs)
        return raw_repr

    @property
    def db_path(self):
        return self.__db_path

    @property
    def verbose(self):
        return self.__verbose

    def __init_db(self):
        self.info()
        self.info("connecting to database: " + self.db_path)
        self.__create_table(tbl_name=TBL_NAME_AVDB_INFO,
                            col_names=[TBL_AVDBS_COL_NAME_COL_NAME,
                                       TBL_AVDBS_TBL_NAME_COL_NAME],
                            pkeys=[TBL_AVDBS_COL_NAME_COL_NAME,
                                   TBL_AVDBS_TBL_NAME_COL_NAME],
                            )
        self.__create_table(tbl_name=TBL_NAME_ANNOVAR_INFO,
                            col_names=[TBL_ANNOVAR_COL_NAME_COL_NAME,
                                       TBL_ANNOVAR_TBL_NAME_COL_NAME],
                            pkeys=[TBL_ANNOVAR_COL_NAME_COL_NAME,
                                   TBL_ANNOVAR_TBL_NAME_COL_NAME],
                            )
        self.__create_table(tbl_name=TBL_NAME_SAMPLES,
                            col_names=[TBL_SAMPLES_SAMPLE_ID_COL_NAME,
                                       TBL_SAMPLES_GTZ_TBL_NAME_COL_NAME],
                            pkeys=[TBL_SAMPLES_SAMPLE_ID_COL_NAME,
                                   TBL_SAMPLES_GTZ_TBL_NAME_COL_NAME],
                            )
    def __log_sql(self, sql):
        if self.verbose:
            self.info("executing sql: " + sql)

    def __fetch_all(self, sql):
        # to encapsulate sqlite sql fetchall execution
        # so far just for logging
        conn = self.__connect_db()
        self.__log_sql(sql)
        result = conn.execute(sql).fetchall()
        self.__close_connection()
        return result

    def __fetch_none(self, sql):
        # to encapsulate sqlite sql fetchnone execution
        # so far just for logging
        conn = self.__connect_db()
        self.__log_sql(sql)
        result = conn.execute(sql).fetchone()
        self.__close_connection()
        return result

    def __exec_sql(self, sql):
        # to encapsulate sqlite sql execution
        # so far just for logging
        conn = self.__connect_db() 
        self.__log_sql(sql)
        result = conn.execute(sql)
        self.__close_connection()
        return result

    def __exec_many_sql(self, sql, records):
        # to encapsulate sqlite sql execution
        # so far just for logging
        conn = self.__connect_db() 
        self.__log_sql(sql)
        result = conn.executemany(sql, records)
        self.__close_connection()
        return result

    def __connect_db(self):
        self.__conn = sqlite3.connect(self.db_path)
        return self.__conn

    def __close_connection(self):
        self.__conn.commit()
        self.__conn.close()

    def drop_table(self, tbl_name):
        sql = 'DROP TABLE IF EXISTS ' + tbl_name
        self.__exec_sql(sql)

    def __drop_view(self, view_name):
        sql = 'DROP VIEW IF EXISTS ' + view_name
        self.__exec_sql(sql)

    def __create_table(self, tbl_name, col_names, pkeys=None):
        sql = 'CREATE TABLE IF NOT EXISTS ' + tbl_name
        sql += ' ( ' + ", ".join(col_names) 
        if pkeys is not None:
            sql += ', PRIMARY KEY (' + ",".join(pkeys) + ')' 
        sql += ' )'
        self.__exec_sql(sql)

    def __add_avdb_column(self, col_name, tbl_name):
        sql = 'INSERT INTO ' + TBL_NAME_AVDB_INFO
        sql += " ( " + TBL_AVDBS_COL_NAME_COL_NAME + ", " + TBL_AVDBS_TBL_NAME_COL_NAME + " ) "
        sql += " select '" + col_name + "', '" + tbl_name + "' "
        sql += " WHERE NOT EXISTS ( "
        sql += " SELECT 1 FROM " + TBL_NAME_AVDB_INFO
        sql += " WHERE " + TBL_AVDBS_COL_NAME_COL_NAME + "='" + col_name + "'"
        sql += " AND " + TBL_AVDBS_TBL_NAME_COL_NAME + "='" + tbl_name + "' )"
        self.__exec_sql(sql)

    def __add_annovar_column(self, col_name, tbl_name):
        sql = 'INSERT INTO ' + TBL_NAME_ANNOVAR_INFO
        sql += " ( " + TBL_ANNOVAR_COL_NAME_COL_NAME + ", " + TBL_ANNOVAR_TBL_NAME_COL_NAME + " ) "
        sql += " select '" + col_name + "', '" + tbl_name + "' "
        sql += " WHERE NOT EXISTS ( "
        sql += " SELECT 1 FROM " + TBL_NAME_ANNOVAR_INFO
        sql += " WHERE " + TBL_ANNOVAR_COL_NAME_COL_NAME + "='" + col_name + "'"
        sql += " AND " + TBL_ANNOVAR_TBL_NAME_COL_NAME + "='" + tbl_name + "' )"
        self.__exec_sql(sql)

    def __add_sample(self, sample_id, tbl_name):
        sql = 'INSERT INTO ' + TBL_NAME_SAMPLES
        sql += " ( " + TBL_SAMPLES_SAMPLE_ID_COL_NAME + ", " + TBL_SAMPLES_GTZ_TBL_NAME_COL_NAME + " ) "
        sql += " select '" + sample_id + "', '" + tbl_name + "' "
        sql += " WHERE NOT EXISTS ( "
        sql += " SELECT 1 FROM " + TBL_NAME_SAMPLES
        sql += " WHERE " + TBL_SAMPLES_SAMPLE_ID_COL_NAME + "='" + sample_id + "'"
        sql += " AND " + TBL_SAMPLES_GTZ_TBL_NAME_COL_NAME + "='" + tbl_name + "' )"
        self.__exec_sql(sql)

    def __load_data_file(self,
                         reader,
                         tbl_name,
                         ):
        rec_len = len(reader.header_cols)
        sql = 'INSERT OR IGNORE INTO ' + tbl_name
        sql += ' VALUES (' + ('?,'*(rec_len*2))[:rec_len*2-1] + ')'
        row_count = self.__exec_many_sql(sql, reader).rowcount
        return row_count

    def load_avdb_data(self, data_file, header_exist, *args, **kwargs):
        reader = AVDBReader(file_name=data_file, header_exist=header_exist)
        tbl_name = kwargs['tbl_name']
        self.__create_table(tbl_name=tbl_name,
                            col_names=reader.header_cols,
                            pkeys=reader.pkeys,
                            )
        for avdb_col in reader.avdb_cols:
            self.__add_avdb_column(avdb_col, tbl_name)
        return self.__load_data_file(reader=reader, *args, **kwargs)

    def load_annovar_vcf_data(self, data_file, *args, **kwargs):
        reader = VcfInfoReader(file_name=data_file)
        tbl_name = kwargs['tbl_name']
        self.__create_table(tbl_name=tbl_name,
                            col_names=reader.header_cols,
                            )
        for info_col in reader.info_cols:
            self.__add_annovar_column(info_col, tbl_name)
        return self.__load_data_file(reader=reader, *args, **kwargs)

    def load_gtz_vcf_data(self, data_file, *args, **kwargs):
        reader = VcfGTZReader(file_name=data_file)
        tbl_name = kwargs['tbl_name']
        self.__create_table(tbl_name=tbl_name,
                            col_names=reader.header_cols,
                            )
        for sample_id in reader.samples_id:
            self.__add_sample(sample_id, tbl_name)
        return self.__load_data_file(reader=reader, *args, **kwargs)

    def __add_column(self, tbl_name, col_name, col_def=""):
        sql = "ALTER TABLE " + tbl_name
        sql += " ADD " + col_name
        sql += " " + col_def
        self.__exec_sql(sql)

    def cal_hw(self, tbl_name):
        # look for dataset name
        col_names = self.get_col_names(tbl_name)
        af_col_names = filter(lambda x: x.endswith("AF"),
                              col_names)
        if len(af_col_names) != 1:
            raise Exception('cannot perform Hardy-Weinberg calculation for table ' + tbl_name)
        dataset_name = af_col_names[0].replace("_AF","")
        self.__add_column(tbl_name, "p")
        self.__add_column(tbl_name, "q")
        self.__add_column(tbl_name, "exp_wt")
        self.__add_column(tbl_name, "exp_het")
        self.__add_column(tbl_name, "exp_hom")
        self.__add_column(tbl_name, "hw_chisq")
        sql = "UPDATE " + tbl_name
        sql += " SET p = 1-" + dataset_name + "_AF"
        sql += ", q = " + dataset_name + "_AF"
        self.__exec_sql(sql)
        sql = "UPDATE " + tbl_name
        sql += " SET exp_wt = p*p*" + dataset_name + "_GT"
        sql += ", exp_het = 2*p*q*" + dataset_name + "_GT"
        sql += ", exp_hom = q*q*" + dataset_name + "_GT"
        self.__exec_sql(sql)
        sql = "UPDATE " + tbl_name
        sql += " SET hw_chisq = 0.0"
        sql += " WHERE exp_wt = 0"
        sql += " OR exp_het = 0"
        sql += " OR exp_hom = 0"
        self.__exec_sql(sql)
        sql = "UPDATE " + tbl_name
        sql += " SET hw_chisq = (" + dataset_name + "_WT-exp_wt)*(" + dataset_name + "_WT-exp_wt)/exp_wt"
        sql += " + (" + dataset_name + "_HET-exp_het)*(" + dataset_name + "_HET-exp_het)/exp_het"
        sql += " + (" + dataset_name + "_HOM-exp_hOM)*(" + dataset_name + "_HOM-exp_hom)/exp_hom"
        sql += " WHERE exp_wt != 0"
        sql += " AND exp_het != 0"
        sql += " AND exp_hom != 0"
        self.__exec_sql(sql)

    def create_mutrep_view(self, view_name, gtz_tbl_name):
        # create view for further used in mutation report
        # it's based on an assumption that all avdb columns in the database
        # will be included into the view
        samples_id = self.get_samples_id(gtz_tbl_name=gtz_tbl_name)
        avdb_info = self.get_avdb_info()
        annovar_info = self.get_annovar_info()
        self.__drop_view(view_name)
        sql = 'CREATE VIEW IF NOT EXISTS ' + view_name + " AS "
        sql += " SELECT gtz.CHROM, gtz.POS, gtz.REF, gtz.ALT, gtz.FILTER"
        for annovar_tbl_name in annovar_info:
            col_names = annovar_info[annovar_tbl_name]
            for col_name in col_names:
                sql += ", " + annovar_tbl_name + "." + col_name
        for avdb_tbl_name in avdb_info:
            col_names = avdb_info[avdb_tbl_name]
            for col_name in col_names:
                sql += ", " + avdb_tbl_name + "." + col_name
        sql += ", " + ", ".join(map(lambda x: "gtz."+x, samples_id))
        sql += " FROM " + gtz_tbl_name + " AS gtz"
        for annovar_tbl_name in annovar_info:
            sql += " LEFT OUTER JOIN " + annovar_tbl_name + " ON "
            join_conditions = []
            join_conditions.append("gtz.CHROM=" + annovar_tbl_name + ".CHROM")
            join_conditions.append("gtz.POS=" + annovar_tbl_name + ".POS")
            join_conditions.append("gtz.REF=" + annovar_tbl_name + ".REF")
            join_conditions.append("gtz.ALT=" + annovar_tbl_name + ".ALT")
            sql += " AND ".join(join_conditions)
        for avdb_tbl_name in avdb_info:
            sql += " LEFT OUTER JOIN " + avdb_tbl_name + " ON "
            join_conditions = []
            join_conditions.append("gtz.CHROM=" + avdb_tbl_name + ".Chr")
            join_conditions.append("gtz.POS=" + avdb_tbl_name + ".Start")
            join_conditions.append("gtz.REF=" + avdb_tbl_name + ".Ref")
            join_conditions.append("gtz.ALT=" + avdb_tbl_name + ".ALt")
            sql += " AND ".join(join_conditions)
        self.__exec_sql(sql)

    def table_exist(self, tbl_name):
        sql = "SELECT name FROM sqlite_master WHERE type='table'"
        sql += "AND name='" + tbl_name + "'"
        tbl_exist = self.__fetch_none(sql)
        if tbl_exist:
            return True
        return False

    def count_rows(self, tbl_name):
        sql = "SELECT COUNT(*) FROM " + tbl_name
        (rows_count,) = self.__fetch_none(sql)
        return rows_count

    def get_avdb_info(self):
        avdb_info = defaultdict(list)
        sql = "SELECT " + TBL_AVDBS_COL_NAME_COL_NAME + ", " + TBL_AVDBS_TBL_NAME_COL_NAME
        sql += " FROM " + TBL_NAME_AVDB_INFO
        rows = self.__fetch_all(sql)
        for (col_name, tbl_name) in rows:
            avdb_info[tbl_name].append(col_name)
        return avdb_info

    def get_annovar_info(self):
        annovar_info = defaultdict(list)
        sql = "SELECT " + TBL_ANNOVAR_COL_NAME_COL_NAME + ", " + TBL_ANNOVAR_TBL_NAME_COL_NAME
        sql += " FROM " + TBL_NAME_ANNOVAR_INFO
        rows = self.__fetch_all(sql)
        for (col_name, tbl_name) in rows:
            annovar_info[tbl_name].append(col_name)
        return annovar_info

    def get_samples_id(self, gtz_tbl_name=None):
        sql = "SELECT " + TBL_SAMPLES_SAMPLE_ID_COL_NAME
        sql += " FROM " + TBL_NAME_SAMPLES
        if gtz_tbl_name is not None:
            sql += " WHERE " + TBL_SAMPLES_GTZ_TBL_NAME_COL_NAME + "='" + gtz_tbl_name + "'"
        rows = self.__fetch_all(sql)
        return map(lambda x: x[0], rows)

    def get_col_names(self, tbl_name):
        cursor = self.__exec_sql('SELECT * FROM ' + tbl_name + " LIMIT 1")
        col_names = map(lambda x: x[0], cursor.description)
        return col_names

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
    """ A class to control all processes related to SQLiteDB """

    def __init__(self, jobs_setup_file, verbose=True, *args, **kwargs):
        kwargs['jobs_setup_file'] = jobs_setup_file
        super(SQLiteDBController, self).__init__(*args, **kwargs)
        self.__db_params = DBParams(self._get_job_config(JOBS_SETUP_DB_PARAMS_SECTION,
                                                         required=True))
        self.__db = SQLiteDB(self.db_params.db_file, verbose=verbose)

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
