# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.proc.db.dbms import SQLiteDB
from pycmm.proc.db.dbms import SQLiteDBController
from pycmm.proc.db.dbms import create_dbms_jobs_setup_file as create_jobs_setup_file


class TestSQLiteDB(SafeTester):

    def __init__(self, methodName):
        super(TestSQLiteDB, self).__init__(methodName=methodName,
                                           test_module_name=__name__,
                                           )

    def test_load_avdb_data_1(self):
        """ test loading avdb data with header and multiple annotations """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "db_input.txt")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        table_name = self.current_func_name
        db = SQLiteDB(db_file)
        row_count = db.load_avdb_data(data_file=input_file,
                                      table_name=table_name,
                                      header_exist=True,
                                      )
        self.assertEqual(row_count,
                         11,
                         "SQLiteDB cannot correctly load avdb data")

    def test_load_annovar_vcf_data_1(self):
        """ test loading annotation file from ANNOVAR """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db = SQLiteDB(db_file)
        row_count = db.load_annovar_vcf_data(data_file=input_file,
                                             )
        self.assertEqual(row_count,
                         12,
                         "SQLiteDB cannot correctly load avdb data")

    def test_load_gtz_vcf_data_1(self):
        """ test loading genotyping zygosities from vcf file """

        self.init_test(self.current_func_name)
        input_file = join_path(self.data_dir,
                               "input.vcf.gz")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        table_name = self.current_func_name
        db = SQLiteDB(db_file)
        row_count = db.load_gtz_vcf_data(data_file=input_file,
                                         table_name=table_name,
                                         )
        self.assertEqual(row_count,
                         6,
                         "SQLiteDB cannot correctly load avdb data")

    def test_calculate_hardy_weinberg(self):
        """
        test calculating Hardy-Weinberg equilibrium in avdb table
        from MAF OAF data
        """

# *********************************************** Require proper testing **********************************************************
# waiting for a good query engine to be implemented
        self.individual_debug = True
        self.init_test(self.current_func_name)
        db_file_with_raw_maf = join_path(self.data_dir,
                                         "input.db")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        self.copy_file(db_file_with_raw_maf,
                       db_file)
        table_name = "test_avdb"
        db = SQLiteDB(db_file)
        db.cal_hw(table_name)
# *********************************************** Require proper testing **********************************************************

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
                                              self.test_function+'_jobs_setup.txt')
        create_jobs_setup_file(*args, **kwargs)
        return kwargs['jobs_setup_file']

    def test_load_jobs_info_1(self):
        """ test if default configurations are loaded correctly """

        self.init_test(self.current_func_name)
        db_file = 'dummy1.db'
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file)
        pl = SQLiteDBController(jobs_setup_file)
        self.assertEqual(pl.db_params.db_file,
                         "dummy1.db",
                         "SQLiteDBController cannot correctly read job info 'db file' from jobs setup file")
        self.assertEqual(len(pl.db_params.db_jobs),
                         0,
                         "SQLiteDBController cannot correctly read job info 'db jobs' from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if configurations with defined parameters are loaded correctly """

        self.init_test(self.current_func_name)
        db_file = 'dummy2.db'
        db_jobs = "table1::data_file1:type1:Yes"
        db_jobs += ",TableN::test_file2:typo:Yes"
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        pl = SQLiteDBController(jobs_setup_file)
        self.assertEqual(pl.db_params.db_file,
                         "dummy2.db",
                         "SQLiteDBController cannot correctly read job info 'db file' from jobs setup file")
        self.assertEqual(len(pl.db_params.db_jobs),
                         2,
                         "SQLiteDBController cannot correctly read job info 'db jobs' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[0].table_name,
                         'table1',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[0].operation,
                         'type1',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[0].drop_table,
                         False,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[1].data_file,
                         'test_file2',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[1].header_exist,
                         True,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")

    def test_load_jobs_info_3(self):
        """
        test configurations exception handling and
        configuration with advanced params
        #1
        """

        self.init_test(self.current_func_name)
        db_file = 'dummy3.db'
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs="",
                                                        )
        pl = SQLiteDBController(jobs_setup_file)
        self.assertEqual(pl.db_params.db_file,
                         "dummy3.db",
                         "SQLiteDBController cannot correctly read job info 'db file' from jobs setup file")
        self.assertEqual(len(pl.db_params.db_jobs),
                         0,
                         "SQLiteDBController cannot correctly read job info 'db jobs' from jobs setup file")

    def test_load_jobs_info_4(self):
        """
        test configurations exception handling and
        configuration with advanced params
        #2
        """

        self.init_test(self.current_func_name)
        db_file = 'dummy4.db'
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        pl = SQLiteDBController(jobs_setup_file)
        self.assertEqual(pl.db_params.db_file,
                         "dummy4.db",
                         "SQLiteDBController cannot correctly read job info 'db file' from jobs setup file")
        self.assertEqual(len(pl.db_params.db_jobs),
                         3,
                         "SQLiteDBController cannot correctly read job info 'db jobs' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[0].table_name,
                         'table1',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[0].drop_table,
                         False,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[1].data_file,
                         'test_file789',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[1].header_exist,
                         None,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[2].data_file,
                         'sss',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[2].operation,
                         'non_type',
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")
        self.assertEqual(pl.db_params.db_jobs[2].drop_table,
                         False,
                         "SQLiteDBController cannot correctly read job info 'db job' from jobs setup file")

    def test_execute_db_jobs_1(self):
        """ test base-line db_jobs execution"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db_jobs = join_path(self.data_dir,
                            'db.jobs')
        jobs_setup_file = self.__create_jobs_setup_file(db_file=db_file,
                                                        db_jobs=db_jobs,
                                                        )
        pl = SQLiteDBController(jobs_setup_file)
        pl.run_offline_pipeline()
        db = SQLiteDB(db_file)
        self.assertEqual(db.table_exist('test_avdb'),
                         True,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.table_exist('test_annovar_vcf'),
                         False,
                         "SQLiteDBController cannot correctly execute db jobs")
        self.assertEqual(db.table_exist('test_gtz_vcf'),
                         True,
                         "SQLiteDBController cannot correctly execute db jobs")
