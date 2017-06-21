# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.proc.mutrep.dbreader import SQLiteDBReader


class TestSQLiteDBReader(SafeTester):

    def __init__(self, methodName):
        super(TestSQLiteDBReader, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def test_mutrep_view_1(self):
        """ test creating VIEW for further use in mutation report """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        self.init_test(self.current_func_name)
        db_file_with_raw_maf = join_path(self.data_dir,
                                         "input.db")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        self.copy_file(db_file_with_raw_maf,
                       db_file)
        db = SQLiteDBReader(db_file, verbose=False)
        view_name = self.current_func_name
        gtz_tbl_name = "test_gtz_THYRCA"
        db.create_mutrep_view(view_name=view_name,
                              gtz_tbl_name=gtz_tbl_name)
        self.assertEqual(len(db.get_col_names(view_name)),
                         147,
                         "SQLiteDB cannot correctly create view")
        self.assertEqual(db.count_rows(view_name),
                         9,
                         "SQLiteDB cannot correctly create view")
