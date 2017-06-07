# This Python file uses the following encoding: utf-8
import unittest
#import filecmp
#import datetime
#from os.path import exists as path_exists
from os.path import join as join_path
#from os.path import dirname
#from os.path import isdir
from pycmm.template import SafeTester
from pycmm.proc.db.dbms import SQLiteDB
#from pycmm.settings import FULL_SYSTEM_TEST
#from pycmm.settings import TEST_PROJECT_CODE
#from pycmm.utils import count_lines
#from pycmm.utils import exec_sh
#from pycmm.flow.cmmdb import CMMDBPipeline
#from pycmm.flow.cmmdb import create_jobs_setup_file
#from pycmm.flow.cmmdb import CONCAT_STAT_CMMDB_SCRIPT
#from pycmm.cmmlib.test.test_annovarlib import DFLT_ANV_TEST_DB_DIR
#from pycmm.cmmlib.test.test_annovarlib import DFLT_ANV_TEST_DB_NAMES
#from pycmm.cmmlib.test.test_annovarlib import DFLT_ANV_TEST_DB_OPS
#from pycmm.cmmlib.samplelib import NO_FAMILY
#
#CMMDB_TEST = False
#DFLT_TEST_DB_REGION = "18"

class TestSQLiteDB(SafeTester):

    def __init__(self, methodName):
        super(TestSQLiteDB, self).__init__(methodName=methodName,
                                           test_module_name=__name__,
                                           )

    def setUp(self):
        pass

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

        self.individual_debug = True
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

    def tearDown(self):
        self.remove_working_dir()
