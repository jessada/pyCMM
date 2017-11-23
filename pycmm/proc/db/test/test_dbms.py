# This Python file uses the following encoding: utf-8
import unittest
import filecmp
from os.path import join as join_path
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.settings import DB_TEST
from pycmm.template import SafeTester
from pycmm.utils import file_size
from pycmm.proc.db.dbinput import AVDBReader
from pycmm.proc.db.dbinput import TAVcfInfoReader as VcfInfoReader
from pycmm.proc.db.dbinput import TAVcfGTZReader as VcfGTZReader
from pycmm.proc.db.dbinput import DATA_TYPE_AVDB_CMM_AF
from pycmm.proc.db.dbinput import DATA_TYPE_AVDB_PUBLIC_AF
from pycmm.proc.db.dbms import SQLiteDBWriter
from pycmm.proc.db.dbms import SQLiteDBController
from pycmm.proc.db.dbms import create_dbms_jobs_setup_file as create_jobs_setup_file
from pycmm.proc.db.connector import TBL_NAME_CMM_TBLS
from pycmm.proc.db.connector import TBL_CMM_TBLS_TBL_NAME_COL_NAME
from pycmm.proc.db.connector import TBL_NAME_GTZ_COORS
from pycmm.proc.db.connector import TBL_NAME_ALL_GTZ_ANNOS
from pycmm.proc.db.connector import SQLiteDB

DBMS_TEST = False


class TestSQLiteDBWriter(SafeTester):

    def __init__(self, methodName):
        super(TestSQLiteDBWriter, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_load_avdb_data_1(self):
        """ test loading avdb data with header and multiple annotations """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "db_input.txt")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        tbl_name = self.current_func_name
        db = SQLiteDBWriter(db_file, verbose=False)
        row_count = db.load_avdb_data(data_file=input_file,
                                      tbl_name=tbl_name,
                                      header_exist=True,
                                      )
        db.drop_table(tbl_name)
        self.assertEqual(row_count,
                         11,
                         "SQLiteDB cannot correctly load avdb data")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_load_avdb_data_2(self):
        """ test loading avdb data with header and multiple annotations """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "db_input.txt")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        tbl_name = self.current_func_name
        db = SQLiteDBWriter(db_file, verbose=False)
        row_count = db.load_avdb_data(data_file=input_file,
                                      tbl_name=tbl_name,
                                      header_exist=True,
                                      data_type=DATA_TYPE_AVDB_CMM_AF,
                                      )
        db.drop_table(tbl_name)
        self.assertEqual(row_count,
                         23,
                         "SQLiteDB cannot correctly load avdb data")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_load_avdb_data_3(self):
        """ test unicode """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "db_input.txt")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        tbl_name = self.current_func_name
        db = SQLiteDBWriter(db_file, verbose=False)
        row_count = db.load_avdb_data(data_file=input_file,
                                      tbl_name=tbl_name,
                                      header_exist=False,
                                      )
        db.drop_table(tbl_name)
        self.assertEqual(row_count,
                         4,
                         "SQLiteDB cannot correctly load avdb data")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_load_avdb_data_4(self):
        """ test avdb w/o header and name start with number """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "hg19_1234_test_avdb.txt")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        tbl_name = self.current_func_name
        db = SQLiteDBWriter(db_file, verbose=False)
        row_count = db.load_avdb_data(data_file=input_file,
                                      tbl_name=tbl_name,
                                      header_exist=True,
                                      data_type=DATA_TYPE_AVDB_PUBLIC_AF,
                                      )
        db.drop_table(tbl_name)
        self.assertEqual(row_count,
                         10,
                         "SQLiteDB cannot correctly load avdb data")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_load_annovar_vcf_data_1(self):
        """ test loading annotation file from ANNOVAR """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        tbl_name = self.current_func_name
        db = SQLiteDBWriter(db_file, verbose=False)
        row_count = db.load_annovar_vcf_data(data_file=input_file,
                                             tbl_name=tbl_name,
                                             )
        db.drop_table(tbl_name)
        self.assertEqual(row_count,
                         12,
                         "SQLiteDB cannot correctly load avdb data")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_load_gtz_vcf_data_1(self):
        """ test loading genotyping zygosities from vcf file """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        tbl_name = self.current_func_name
        db = SQLiteDBWriter(db_file, verbose=False)
        row_count = db.load_gtz_vcf_data(data_file=input_file,
                                         tbl_name=tbl_name,
                                         )
        db.drop_table(tbl_name)
        self.assertEqual(row_count,
                         6,
                         "SQLiteDB cannot correctly load avdb data")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_calculate_hardy_weinberg_1(self):
        """
        test calculating Hardy-Weinberg equilibrium in avdb table
        from MAF OAF data
        """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "db_input.txt")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        tbl_name = self.current_func_name
        db = SQLiteDBWriter(db_file, verbose=False)
        db.drop_table(tbl_name)
        row_count = db.load_avdb_data(data_file=input_file,
                                      tbl_name=tbl_name,
                                      header_exist=True,
                                      )
        self.assertEqual(row_count,
                         11,
                         "SQLiteDB cannot correctly load avdb data")
        avdb_info = db.get_avdb_info()
        self.assertEqual(len(avdb_info[tbl_name]),
                         9,
                         "SQLiteDB cannot correctly calculating Hardy-Weinberg ")
        db.cal_hw(tbl_name)
        avdb_info = db.get_avdb_info()
        self.assertEqual(len(avdb_info[tbl_name]),
                         10,
                         "SQLiteDB cannot correctly calculating Hardy-Weinberg ")
        row_count = db.count_rows(TBL_NAME_CMM_TBLS,
                                  TBL_CMM_TBLS_TBL_NAME_COL_NAME+"= '"+tbl_name+"'")
        self.assertEqual(row_count,
                         1,
                         "SQLiteDB cannot correctly load avdb data")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_calculate_hardy_weinberg_2(self):
        """
        test calculating Hardy-Weinberg equilibrium in avdb table
        from Swegen data
        """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "db_input.txt")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        tbl_name = self.current_func_name
        db = SQLiteDBWriter(db_file, verbose=False)
        db.drop_table(tbl_name)
        row_count = db.load_avdb_data(data_file=input_file,
                                      tbl_name=tbl_name,
                                      header_exist=True,
                                      )
        self.assertEqual(row_count,
                         125,
                         "SQLiteDB cannot correctly load avdb data")
        avdb_info = db.get_avdb_info()
        self.assertEqual(len(avdb_info[tbl_name]),
                         8,
                         "SQLiteDB cannot correctly calculating Hardy-Weinberg ")
        db.cal_hw(tbl_name)
        avdb_info = db.get_avdb_info()
        self.assertEqual(len(avdb_info[tbl_name]),
                         10,
                         "SQLiteDB cannot correctly calculating Hardy-Weinberg ")
        row_count = db.count_rows(TBL_NAME_CMM_TBLS,
                                  TBL_CMM_TBLS_TBL_NAME_COL_NAME+"= '"+tbl_name+"'")
        self.assertEqual(row_count,
                         1,
                         "SQLiteDB cannot correctly load avdb data")

class TestSQLiteDBController(SafeTester):

    def __init__(self, methodName):
        super(TestSQLiteDBController, self).__init__(methodName=methodName,
                                                     test_module_name=__name__,
                                                     )

    def __create_jobs_setup_file(self, *args, **kwargs):
        if 'project_name' not in kwargs:
            kwargs['project_name'] = self.test_function
        kwargs['project_out_dir'] = self.working_dir
        kwargs['jobs_setup_file'] = join_path(self.working_dir,
                                              kwargs['project_name']+'_jobs_setup.txt')
        create_jobs_setup_file(*args, **kwargs)
        return kwargs['jobs_setup_file']

    def test_load_jobs_info_1(self):
        """ test if default configurations are loaded correctly """

        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            'dummy1.db')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file)
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        self.assertEqual(dbc.db_params.db_file,
                         db_file,
                         "SQLiteDBController cannot correctly read job info 'db file' from jobs setup file")
        self.assertEqual(len(dbc.db_params.db_jobs),
                         0,
                         "SQLiteDBController cannot correctly read job info 'db jobs' from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if configurations with defined parameters are loaded correctly """

        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            'dummy2.db')
        db_jobs = "table1::data_file1:type1:Yes"
        db_jobs += ",TableN::test_file2:typo:Yes"
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        self.assertEqual(dbc.db_params.db_file,
                         db_file,
                         "SQLiteDBController cannot correctly read job info 'db file' from jobs setup file")
        self.assertEqual(len(dbc.db_params.db_jobs),
                         2,
                         "SQLiteDBController cannot correctly read job info 'db jobs' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[0].tbl_name,
                         'table1',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[0].operation,
                         'type1',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[0].drop_table,
                         False,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[1].data_file,
                         'test_file2',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[1].header_exist,
                         True,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")

    def test_load_jobs_info_3(self):
        """
        test configurations exception handling and
        configuration with advanced params
        #1
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            'dummy3.db')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs="",
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        self.assertEqual(dbc.db_params.db_file,
                         db_file,
                         "SQLiteDBController cannot correctly read job info 'db file' from jobs setup file")
        self.assertEqual(len(dbc.db_params.db_jobs),
                         0,
                         "SQLiteDBController cannot correctly read job info 'db jobs' from jobs setup file")

    def test_load_jobs_info_4(self):
        """
        test configurations exception handling and
        configuration with advanced params
        #2
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            'dummy4.db')
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        self.assertEqual(dbc.db_params.db_file,
                         db_file,
                         "SQLiteDBController cannot correctly read job info 'db file' from jobs setup file")
        self.assertEqual(len(dbc.db_params.db_jobs),
                         3,
                         "SQLiteDBController cannot correctly read job info 'db jobs' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[0].tbl_name,
                         'test_avdb',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[0].drop_table,
                         True,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        exp_data_file = join_path(self.data_dir,
                                  "test_avdb.txt")
        self.assertTrue(filecmp.cmp(exp_data_file, dbc.db_params.db_jobs[0].data_file),
                        "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        exp_data_file = join_path(self.data_dir,
                                  "test_annovar_vcf.vcf.gz")
        self.assertTrue(filecmp.cmp(exp_data_file, dbc.db_params.db_jobs[1].data_file),
                        "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[1].data_type,
                         DATA_TYPE_AVDB_CMM_AF,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[1].header_exist,
                         None,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        exp_data_file = join_path(self.data_dir,
                                  "test_gtz_vcf.vcf.gz")
        self.assertTrue(filecmp.cmp(exp_data_file, dbc.db_params.db_jobs[2].data_file),
                        "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[2].data_type,
                         DATA_TYPE_AVDB_PUBLIC_AF,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[2].operation,
                         'LOAD_GTZ_VCF',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(dbc.db_params.db_jobs[2].drop_table,
                         False,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_execute_db_jobs_1(self):
        """ test base-line db_jobs execution """

        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        db = SQLiteDBWriter(db_file, verbose=False)
        self.assertEqual(db.table_exist('test_avdb'),
                         True,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.table_exist('annovar_vcf'),
                         False,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.table_exist('test_gtz_vcf'),
                         True,
                         "SQLiteDBController cannot correctly execute db jobs")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_execute_db_jobs_2(self):
        """
        test db_jobs execution with data that
        are likely to represent real input data
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        db = SQLiteDB(db_file, verbose=False)
        # record counting to ensure if the settings (drop_table for example) are working correct
        self.assertEqual(db.count_rows("test_hg19_gnomad_genome"),
                         1010,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_clinvar_20150330"),
                         8,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_clinvar_20170130"),
                         29,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_Axeq_chr5_19"),
                         0,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_OAF_BrC_CRC_prostate"),
                         21,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_OAF_familial_CRCs"),
                         21,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_OAF_CHEK2"),
                         21,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_Swegen_20161223"),
                         133,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_spidex"),
                         1935,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_intervar"),
                         27,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_dbnsfp30a"),
                         128,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_avsnp144"),
                         583,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_annovar_vcf"),
                         18,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_gtz_THYRCA"),
                         10,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_gtz_WES294"),
                         13,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(reduce(lambda x,y: x+y, db.get_avdb_info().values())),
                         171,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_samples_id()),
                         304,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows(TBL_NAME_GTZ_COORS),
                         30,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_col_names(TBL_NAME_ALL_GTZ_ANNOS)),
                         495,
                         "SQLiteDBController cannot correctly execute db jobs")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_execute_db_jobs_2(self):
        """
        test db_jobs execution with data that
        are likely to represent real input data
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        db = SQLiteDB(db_file, verbose=False)
        # record counting to ensure if the settings (drop_table for example) are working correct
        self.assertEqual(db.count_rows("test_hg19_gnomad_genome"),
                         1010,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_clinvar_20150330"),
                         8,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_clinvar_20170130"),
                         29,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_Axeq_chr5_19"),
                         0,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_OAF_BrC_CRC_prostate"),
                         21,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_OAF_familial_CRCs"),
                         21,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_OAF_CHEK2"),
                         21,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_CMM_Swegen_20161223"),
                         133,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_spidex"),
                         1935,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_intervar"),
                         27,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_dbnsfp30a"),
                         128,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_hg19_avsnp144"),
                         583,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_annovar_vcf"),
                         18,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_gtz_THYRCA"),
                         10,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows("test_gtz_WES294"),
                         13,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(reduce(lambda x,y: x+y, db.get_avdb_info().values())),
                         171,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_samples_id()),
                         304,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows(TBL_NAME_GTZ_COORS),
                         30,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_col_names(TBL_NAME_ALL_GTZ_ANNOS)),
                         496,
                         "SQLiteDBController cannot correctly execute db jobs")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_execute_db_jobs_3(self):
        """
        test db_jobs execution with data that
        are likely to represent real input data.
        This test is designed to validate and generate test data for
        downstream testing.
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(len(reduce(lambda x,y: x+y, db.get_avdb_info().values())),
                         171,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_samples_id()),
                         304,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows(TBL_NAME_GTZ_COORS),
                         16,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_col_names(TBL_NAME_ALL_GTZ_ANNOS)),
                         496,
                         "SQLiteDBController cannot correctly execute db jobs")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_execute_db_jobs_4(self):
        """
        Test db_jobs execution with data that
        are likely to represent real input data.
        This test is designed to validate and generate test data for
        downstream testing.
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(len(reduce(lambda x,y: x+y, db.get_avdb_info().values())),
                         171,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_samples_id()),
                         304,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows(TBL_NAME_GTZ_COORS),
                         44,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_col_names(TBL_NAME_ALL_GTZ_ANNOS)),
                         496,
                         "SQLiteDBController cannot correctly execute db jobs")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_execute_db_jobs_5(self):
        """
        Test db_jobs execution with data that
        are likely to represent real input data.
        This test is designed to validate and generate test data for
        downstream testing.
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(len(reduce(lambda x,y: x+y, db.get_avdb_info().values())),
                         171,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_samples_id()),
                         304,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.count_rows(TBL_NAME_GTZ_COORS),
                         3,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(len(db.get_col_names(TBL_NAME_ALL_GTZ_ANNOS)),
                         496,
                         "SQLiteDBController cannot correctly execute db jobs")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_execute_db_jobs_6(self):
        """
        Test db_jobs execution with data that
        are likely to represent real input data.
        This test is designed to validate and generate test data for
        downstream testing.
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(len(reduce(lambda x,y: x+y, db.get_avdb_info().values())),
                         51,
                         "SQLiteDBController cannot correctly execute db jobs")

    @unittest.skipUnless(FULL_SYSTEM_TEST or DBMS_TEST or DB_TEST, "taking too long time to test")
    def test_vacuum_1(self):
        """ test if vacuum can shrink sqlite wih many inserts and drops"""

        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db_wo_vacuum.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        self.assertTrue(file_size(db_file) > 460000,
                        "SQLiteDBController cannot correctly VACUUM")

        db_jobs = join_path(self.data_dir,
                            'db_w_vacuum.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbc.run_offline_pipeline()
        self.assertTrue(file_size(db_file) < 460000,
                        "SQLiteDBController cannot correctly VACUUM")
