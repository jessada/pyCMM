#import pyaml
#from os.path import join as join_path
#from collections import OrderedDict
import sqlite3
from pycmm.template import pyCMMBase
#from pycmm.utils import exec_sh
#from pycmm.utils.ver import VersionManager
#from pycmm.cmmlib import CMMParams
from pycmm.proc import CMMPipeline
from pycmm.proc import init_jobs_setup_file
from pycmm.proc import get_func_arg
from pycmm.proc.db.dbinput import AVDBReader
from pycmm.proc.db.dbinput import TAVcfInfoReader as VcfInfoReader
from pycmm.proc.db.dbinput import TAVcfGTZReader as VcfGTZReader
#from pycmm.cmmlib.annovarlib import AnnovarParams
#from pycmm.cmmlib.annovarlib import get_annovar_params_sections
#from pycmm.cmmlib.dnalib import ALL_CHROMS
#from pycmm.settings import DUMMY_TABLE_ANNOVAR_BIN

## *************** mustat db section ***************
#JOBS_SETUP_MUTSTAT_PARAMS_SECTION = "MUTSTAT_PARAMS"
#JOBS_SETUP_VCF_TABIX_FILE_KEY = "VCF_TABIX_FILE"
#JOBS_SETUP_DB_REGION_KEY = "DB_REGION"
#JOBS_SETUP_VCF2AVDB_KEY_TABLE_KEY = "VCF2AVDB_KEY_TABLE"
#
## *************** table_annovar db section ***************
#JOBS_SETUP_ANV_PARAMS_SECTION = "ANNOVAR_PARAMS"

#class MutStatParams(CMMParams):
#    """ To handle and parse mutation statistics parameters """
#
#    def __init__(self, **kwargs):
#        super(MutStatParams, self).__init__(**kwargs)
#
#    def get_raw_obj_str(self, **kwargs):
#        raw_repr = super(MutStatParams, self).get_raw_obj_str(**kwargs)
#        raw_repr["input vcf tabix file"] = self.input_vcf_tabix
#        db_region = self.db_region
#        if db_region is None:
#            db_region = "ALL"
#        raw_repr["input vcf database region"] = db_region
#        return raw_repr
#
#    @property
#    def input_vcf_tabix(self):
#        return self._get_job_config(JOBS_SETUP_VCF_TABIX_FILE_KEY,
#                                    required=True)
#
#    @property
#    def db_region(self):
#        return self._get_job_config(JOBS_SETUP_DB_REGION_KEY)
#
#    @property
#    def vcf2avdb_key_table(self):
#        return self._get_job_config(JOBS_SETUP_VCF2AVDB_KEY_TABLE_KEY)

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

#    def load_avdb_data(self, data_file, table_name, header_exist=False, drop_table=True):
#        rd = AVDBReader(file_name=data_file, header_exist=header_exist)
#        rec_len = len(rd.header_cols)
#        c = self.__connect_db()
#        if drop_table:
#            sql_script = 'DROP TABLE if exists ' + table_name
#            c.execute(sql_script)
#        sql_script = 'CREATE TABLE if not exists ' + table_name
#        sql_script += ' (' + ",".join(rd.header_cols) + ')' 
#        c.execute(sql_script)
#        sql_script = 'INSERT INTO ' + table_name
#        sql_script += ' VALUES (' + ('?,'*(rec_len*2))[:rec_len*2-1] + ')'
#        c.executemany(sql_script, rd)
#        row_count = c.rowcount
#        self.__close_connection() 
#        return row_count
#
#    def load_annovar_vcf_data(self,
#                              data_file,
#                              drop_table=True,
#                              ):
#        table_name = 'annovar_vcf'
#        rd = VcfInfoReader(file_name=data_file)
#        rec_len = len(rd.header_cols)
#        c = self.__connect_db()
#        if drop_table:
#            sql_script = 'DROP TABLE if exists ' + table_name
#            c.execute(sql_script)
#        sql_script = 'CREATE TABLE if not exists ' + table_name
#        sql_script += ' (' + ",".join(rd.header_cols) + ')' 
#        c.execute(sql_script)
#        sql_script = 'INSERT INTO ' + table_name
#        sql_script += ' VALUES (' + ('?,'*(rec_len*2))[:rec_len*2-1] + ')'
#        c.executemany(sql_script, rd)
#        row_count = c.rowcount
#        self.__close_connection() 
#        return row_count

#    def load_avdb_data(self, data_file, table_name, header_exist=False,):
#        return self.load_data(data_file, AVDBReader, table_name, header_exist)

#    def count_row(self, table_name):
#        # it works but people said that it has to go through the whole table
#        # might be very poor performance
#        c = self.__connect_db()
#        sql_script = 'SELECT Count(*) FROM ' + table_name
#        c.execute(sql_script)
#        row_count = c.fetchone()[0]
#        self.__close_connection() 
#        return row_count
#    @property
#    def mutstat_params(self):
#        return MutStatParams(entries=self._get_job_config(JOBS_SETUP_MUTSTAT_PARAMS_SECTION,
#                                                          required=True))
#
#    @property
#    def annovar_params(self):
#        return AnnovarParams(entries=self._get_job_config(JOBS_SETUP_ANV_PARAMS_SECTION,
#                                                          required=True))
#
#    @property
#    def out_files(self):
#        return self.__out_files
#
#    def get_third_party_software_version(self):
#        vm = VersionManager()
#        versions = OrderedDict()
#        versions['vcftools'] = vm.vcftools_version
#        versions['table_annovar.pl'] = vm.table_annovar_version
#        return versions
#
#    def __get_cmmdb_script_params(self,
#                                  dataset_name,
#                                  out_file,
#                                  db_region=None,
#                                  samples_id=None,
#                                  ):
#        params = " -k " + dataset_name
#        params += " -i " + self.mutstat_params.input_vcf_tabix
#        if self.mutstat_params.vcf2avdb_key_table is not None:
#            params += " -t " + self.mutstat_params.vcf2avdb_key_table
#        if db_region is not None:
#            params += " -r " + db_region
#        if samples_id is not None:
#            if type(samples_id) is list:
#                params += " -c " + ",".join(samples_id)
#            else:
#                params += " -c " + samples_id
#        params += " -o " + out_file
#        return params
#
#    def __run_cmm_script(self,
#                         script_name,
#                         out_suffix=".stat",
#                         job_name_suffix="_cal_stat",
#                         ):
#        if self.project_code is None:
#            self.__out_files = join_path(self.data_out_dir,
#                                              self.dataset_name + out_suffix)
#            params = self.__get_cmmdb_script_params(dataset_name=self.dataset_name,
#                                                    out_file=self.out_files,
#                                                    db_region=self.mutstat_params.db_region,
#                                                    samples_id=self.samples_id,
#                                                    )
#            cmd = script_name + params
#            exec_sh(cmd)
#        elif self.mutstat_params.db_region is None:
#            self.__out_files = []
#            for chrom in ALL_CHROMS:
#                dataset_name = self.dataset_name + "_" + chrom
#                out_file = join_path(self.data_out_dir,
#                                          dataset_name + out_suffix)
#                params = self.__get_cmmdb_script_params(dataset_name=dataset_name,
#                                                        out_file=out_file,
#                                                        db_region=chrom,
#                                                        samples_id=self.samples_id,
#                                                        )
#                self.__out_files.append(out_file)
#                job_name = dataset_name + job_name_suffix
#                self._submit_slurm_job(job_name,
#                                       "1",
#                                       script_name,
#                                       params,
#                                       )
#        else:
#            self.__out_files = join_path(self.data_out_dir,
#                                              self.dataset_name + out_suffix)
#            params = self.__get_cmmdb_script_params(dataset_name=self.dataset_name,
#                                                    out_file=self.out_files,
#                                                    db_region=self.mutstat_params.db_region,
#                                                    samples_id=self.samples_id,
#                                                    )
#            job_name = self.dataset_name + job_name_suffix
#            self._submit_slurm_job(job_name,
#                                   "1",
#                                   script_name,
#                                   params,
#                                   )
#
#    def cal_mut_stat(self):
#        self.__run_cmm_script(CAL_MUTATIONS_STAT_CMMDB_SCRIPT)
#
#    def vcfaf_to_annovar(self):
#        self.__run_cmm_script(VCF_AF_TO_ANNOVAR_CMMDB_SCRIPT)
#
#    def vcf2annovar_key(self):
#        self.__run_cmm_script(VCF2AVDB_KEY_CMMDB_SCRIPT,
#                              out_suffix=".key",
#                              job_name_suffix="_key",
#                              )
#
#    def table_annovar(self):
#        # calling dummy table_annovar mainly because to clariy the configurations and for logging
#        annovar_params = self.annovar_params
#        table_annovar_cmd = DUMMY_TABLE_ANNOVAR_BIN
#        table_annovar_params = " --dataset_name " + self.dataset_name
#        table_annovar_params += " --input_file " + annovar_params.input_vcf_tabix
#        table_annovar_params += " --db_folder " + annovar_params.db_dir
#        table_annovar_params += " --buildver " + annovar_params.buildver
#        table_annovar_params += " --protocols " + annovar_params.protocols
#        table_annovar_params += " --operations " + annovar_params.operations
#        table_annovar_params += " --nastring " + annovar_params.nastring
#        table_annovar_params += " --data_out_folder " + self.data_out_dir
#        if self.project_code is not None:
#            job_name = self.dataset_name + "_ta"
#            self._submit_slurm_job(job_name,
#                                   "8",
#                                   table_annovar_cmd,
#                                   table_annovar_params,
#                                   )
#        else:
#            exec_sh(table_annovar_cmd + table_annovar_params)


def create_dbms_jobs_setup_file(*args, **kwargs):
    job_setup_document, stream = init_jobs_setup_file(*args, **kwargs)

#    mutstat_params = {}
#    vcf_tabix_file = get_func_arg('vcf_tabix_file', kwargs)
#    if vcf_tabix_file is not None:
#        mutstat_params[JOBS_SETUP_VCF_TABIX_FILE_KEY] = vcf_tabix_file
#    db_region = get_func_arg('db_region', kwargs)
#    if db_region is not None:
#        mutstat_params[JOBS_SETUP_DB_REGION_KEY] = '"' + db_region + '"'
#    job_setup_document[JOBS_SETUP_MUTSTAT_PARAMS_SECTION] = mutstat_params
#
#    vcf2avdb_key_table = get_func_arg('vcf2avdb_key_table', kwargs)
#    if vcf2avdb_key_table is not None:
#        mutstat_params[JOBS_SETUP_VCF2AVDB_KEY_TABLE_KEY] = vcf2avdb_key_table
#
#    anv_params = get_annovar_params_sections(*args, **kwargs)
    job_setup_document[JOBS_SETUP_ANV_PARAMS_SECTION] = anv_params
    return job_setup_document, stream

## Right now no one actually want to use create_cmmdb_jobs_setup_file
## But reserve it here in case it need to be a complete report and database pipeline
#def create_jobs_setup_file(*args, **kwargs):
#    job_setup_document, stream = create_cmmdb_jobs_setup_file(*args, **kwargs)
#    pyaml.dump(job_setup_document, stream)
