# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.proc.db.connector import SQLiteDB
from pycmm.proc.db.connector import DATA_TYPE_AVDB_INFO
from pycmm.proc.db.connector import DATA_TYPE_GTZ


class TestSQLiteDB(SafeTester):

    def __init__(self, methodName):
        super(TestSQLiteDB, self).__init__(methodName=methodName,
                                           test_module_name=__name__,
                                           )

    def test_update_sys_info_1(self):
        """ test if the system info can written and read """

        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db = SQLiteDB(db_file, verbose=False)
        tbl_name1 = self.current_func_name + "_avdb_1"
        cols_name = []
        cols_name.append("OAF1")
        cols_name.append("OAF2")
        db._update_sys_info(data_type=DATA_TYPE_AVDB_INFO,
                            updating_tbl_name=tbl_name1,
                            updating_info_cols=cols_name,
                            )
        tbl_name2 = self.current_func_name + "_avdb_2"
        cols_name = []
        cols_name.append("OAF4")
        cols_name.append("OAF5")
        cols_name.append("OAF6")
        db._update_sys_info(data_type=DATA_TYPE_AVDB_INFO,
                            updating_tbl_name=tbl_name2,
                            updating_info_cols=cols_name,
                            )
        avdb_info = db.get_avdb_info()
        self.assertEqual(len(avdb_info[tbl_name1]),
                         2,
                         "SQLiteDB cannot correctly update system info")
        self.assertEqual(avdb_info[tbl_name2][1],
                         "OAF5",
                         "SQLiteDB cannot correctly update system info")
        tbl_name3 = self.current_func_name + "_samples_1"
        cols_name = []
        cols_name.append("Co-1")
        cols_name.append("Co-4")
        cols_name.append("Br-1")
        db._update_sys_info(data_type=DATA_TYPE_GTZ,
                            updating_tbl_name=tbl_name3,
                            updating_info_cols=cols_name,
                            )
        tbl_name4 = self.current_func_name + "_samples_2"
        cols_name = []
        cols_name.append("Co-6")
        cols_name.append("Co-7")
        cols_name.append("Br-8")
        cols_name.append("Br-9")
        db._update_sys_info(data_type=DATA_TYPE_GTZ,
                            updating_tbl_name=tbl_name4,
                            updating_info_cols=cols_name,
                            )
        samples_id = db.get_samples_id()
        self.assertEqual(len(samples_id),
                         7,
                         "SQLiteDB cannot correctly update system info")
        self.assertEqual(samples_id[5],
                         "Br-8",
                         "SQLiteDB cannot correctly update system info")
        samples_id = db.get_samples_id(gtz_tbl_name=tbl_name4)
        self.assertEqual(len(samples_id),
                         4,
                         "SQLiteDB cannot correctly update system info")
        self.assertEqual(samples_id[3],
                         "Br-9",
                         "SQLiteDB cannot correctly update system info")
