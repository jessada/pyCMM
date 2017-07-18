# This Python file uses the following encoding: utf-8
import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.settings import DB_TEST
from pycmm.settings import MAX_REF_MAF_COL_NAME
from pycmm.settings import GNOMAD_GENOME_ALL_COL_NAME
from pycmm.proc.db.connector import SQLiteDB
from pycmm.proc.db.connector import DATA_TYPE_AVDB_INFO
from pycmm.proc.db.connector import DATA_TYPE_GTZ
from pycmm.proc.db.connector import TBL_NAME_GTZ_COORS
from pycmm.proc.db.connector import TBL_NAME_ALL_GTZ_ANNOS
from pycmm.proc.db.connector import REF_MUTATED_COL_NAME

CONNECTOR_TEST = False
CONNECTOR_TEST = True


class TestSQLiteDB(SafeTester):

    def __init__(self, methodName):
        super(TestSQLiteDB, self).__init__(methodName=methodName,
                                           test_module_name=__name__,
                                           )

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_update_tbl_info_1(self):
        """ test if the system info can written and read """

        self.init_test(self.current_func_name)
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        db = SQLiteDB(db_file, verbose=False)
        tbl_name1 = self.current_func_name + "_avdb_1"
        cols_name = []
        cols_name.append("OAF1")
        cols_name.append("OAF2")
        db._update_tbl_info(data_type=DATA_TYPE_AVDB_INFO,
                            updating_tbl_name=tbl_name1,
                            updating_info_cols=cols_name,
                            )
        tbl_name2 = self.current_func_name + "_avdb_2"
        cols_name = []
        cols_name.append("OAF4")
        cols_name.append("OAF5")
        cols_name.append("OAF6")
        db._update_tbl_info(data_type=DATA_TYPE_AVDB_INFO,
                            updating_tbl_name=tbl_name2,
                            updating_info_cols=cols_name,
                            )
        avdb_info = db.get_avdb_info()
        self.assertEqual(len(avdb_info[tbl_name1]),
                         2,
                         "SQLiteDB cannot correctly update system table info")
        self.assertEqual(avdb_info[tbl_name2][1],
                         "OAF5",
                         "SQLiteDB cannot correctly update system table info")
        tbl_name3 = self.current_func_name + "_samples_1"
        cols_name = []
        cols_name.append("Co-1")
        cols_name.append("Co-4")
        cols_name.append("Br-1")
        db._update_tbl_info(data_type=DATA_TYPE_GTZ,
                            updating_tbl_name=tbl_name3,
                            updating_info_cols=cols_name,
                            )
        tbl_name4 = self.current_func_name + "_samples_2"
        cols_name = []
        cols_name.append("Co-6")
        cols_name.append("Co-7")
        cols_name.append("Br-8")
        cols_name.append("Br-9")
        db._update_tbl_info(data_type=DATA_TYPE_GTZ,
                            updating_tbl_name=tbl_name4,
                            updating_info_cols=cols_name,
                            )
        samples_id = db.get_samples_id()
        self.assertEqual(len(samples_id),
                         7,
                         "SQLiteDB cannot correctly update system table info")
        samples_id = db.get_samples_id(gtz_tbl_name=tbl_name4)
        self.assertEqual(len(samples_id),
                         4,
                         "SQLiteDB cannot correctly update system table info")
        self.assertEqual(samples_id[3],
                         "Br-9",
                         "SQLiteDB cannot correctly update system table info")
        cols_name.append("Br-20")
        db._update_tbl_info(data_type=DATA_TYPE_GTZ,
                            updating_tbl_name=tbl_name4,
                            updating_info_cols=cols_name,
                            )
        samples_id = db.get_samples_id()
        self.assertEqual(len(samples_id),
                         8,
                         "SQLiteDB cannot correctly update system table info")
        samples_id = db.get_samples_id(gtz_tbl_name=tbl_name4)
        self.assertEqual(len(samples_id),
                         5,
                         "SQLiteDB cannot correctly update system table info")
        self.assertEqual(samples_id[4],
                         "Br-20",
                         "SQLiteDB cannot correctly update system table info")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_join_all_gtz_anno_1(self):
        """ test joining all gtz tables and all annotation tables """

        self.init_test(self.current_func_name)
        raw_db_file = join_path(self.data_dir,
                                "input.db")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        self.copy_file(raw_db_file,
                       db_file)
        db = SQLiteDB(db_file, verbose=False)
        row_count = db.join_all_gtz_anno()
        self.assertEqual(row_count,
                         18,
                         "SQLiteDB cannot correctly join all genotype and annotations")
        self.assertEqual(len(db.get_col_names(TBL_NAME_ALL_GTZ_ANNOS)),
                         417,
                         "SQLiteDB cannot correctly join all genotype and annotations")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_anno_col_to_db_col_1(self):
        """
        test converting ui col name to db col name
        - ui col: gnomAD_genome_ALL
        - db col: _gnomAD_genome_ALL
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            "input.db")
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(db.anno_col_to_db_col(GNOMAD_GENOME_ALL_COL_NAME),
                         "_"+GNOMAD_GENOME_ALL_COL_NAME,
                         "SQLiteDB cannot correctly map UI column name to DB column name")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_anno_col_to_db_col_2(self):
        """
        test converting ui col name to db col name
        - ui col: gnomAD_genome_ALL
        - db col: gnomAD_genome_ALL
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            "input.db")
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(db.anno_col_to_db_col(GNOMAD_GENOME_ALL_COL_NAME),
                         GNOMAD_GENOME_ALL_COL_NAME,
                         "SQLiteDB cannot correctly map UI column name to DB column name")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_anno_col_to_db_col_3(self):
        """
        test converting ui column that doesn't exist
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            "input.db")
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(db.anno_col_to_db_col("NOT_EXIST"),
                         None,
                         "SQLiteDB cannot correctly map UI column name to DB column name")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_sample_id_to_db_col_1(self):
        """ test converting sample id to db column name """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            "input.db")
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(db.sample_id_to_db_col("Co-35"),
                         "_Co_35",
                         "SQLiteDB cannot correctly map sample id to DB column name")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_sample_id_to_tbl_1(self):
        """ test converting sample id to db column name """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            "input.db")
        db = SQLiteDB(db_file, verbose=False)
        self.assertEqual(db.sample_id_to_tbl("Co-35"),
                         "gtz_WES294",
                         "SQLiteDB cannot correctly map sample id to table name")
        self.assertEqual(db.sample_id_to_tbl("Co-309"),
                         "fam119_new_sample",
                         "SQLiteDB cannot correctly map sample id to table name")
        self.assertEqual(db.sample_id_to_tbl("2016-18116"),
                         "gtz_THYRCA",
                         "SQLiteDB cannot correctly map sample id to table name")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_cal_max_ref_maf_1(self):
        """
        test finding maximum allele frequency
        from all reference databases
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        raw_db_file = join_path(self.data_dir,
                                "input.db")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        self.copy_file(raw_db_file,
                       db_file)
        db = SQLiteDB(db_file, verbose=False)
        db.cal_max_ref_maf()
        sql = "SELECT " + MAX_REF_MAF_COL_NAME
        sql += " FROM " + TBL_NAME_ALL_GTZ_ANNOS
        rows = db.read_rows(sql=sql)
        row = rows.next()
        self.assertEqual(row[0],
                         0.116,
                         "SQLiteDB cannot correctly calculate max reference allele frequency")
        row = rows.next()
        row = rows.next()
        row = rows.next()
        self.assertEqual(row[0],
                         0.0382,
                         "SQLiteDB cannot correctly calculate max reference allele frequency")
        row = rows.next()
        row = rows.next()
        row = rows.next()
        self.assertEqual(row[0],
                         0,
                         "SQLiteDB cannot correctly calculate max reference allele frequency")
        row = rows.next()
        row = rows.next()
        row = rows.next()
        self.assertEqual(row[0],
                         0.0565,
                         "SQLiteDB cannot correctly calculate max reference allele frequency")
        row = rows.next()
        row = rows.next()
        row = rows.next()
        self.assertEqual(row[0],
                         0.3504,
                         "SQLiteDB cannot correctly calculate max reference allele frequency")
        row = rows.next()
        row = rows.next()
        row = rows.next()
        self.assertEqual(row[0],
                         0,
                         "SQLiteDB cannot correctly calculate max reference allele frequency")
        row = rows.next()
        row = rows.next()
        row = rows.next()
        self.assertEqual(round(row[0], 4),
                         0.0099,
                         "SQLiteDB cannot correctly calculate max reference allele frequency")

    @unittest.skipUnless(FULL_SYSTEM_TEST or CONNECTOR_TEST or DB_TEST, "taking too long time to test")
    def test_set_ref_mutated_1(self):
        """
        test finding maximum allele frequency
        from all reference databases
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        raw_db_file = join_path(self.data_dir,
                                "input.db")
        db_file = join_path(self.working_dir,
                            self.current_func_name+".db")
        self.copy_file(raw_db_file,
                       db_file)
        db = SQLiteDB(db_file, verbose=False)
        db.set_ref_mutated()
        sql = "SELECT " + REF_MUTATED_COL_NAME
        sql += " FROM " + TBL_NAME_ALL_GTZ_ANNOS
        sql += " WHERE POS > 89684278"
        sql += " OR POS = 57427955"
        rows = db.read_rows(sql=sql)
        row = rows.next()
        self.assertEqual(row[0],
                         0,
                         "SQLiteDB cannot correctly identify mutation in reference")
        row = rows.next()
        row = rows.next()
        row = rows.next()
        self.assertEqual(row[0],
                         1,
                         "SQLiteDB cannot correctly identify mutation in reference")
        row = rows.next()
        self.assertEqual(row[0],
                         0,
                         "SQLiteDB cannot correctly identify mutation in reference")
        row = rows.next()
        row = rows.next()
        self.assertEqual(row[0],
                         0,
                         "SQLiteDB cannot correctly identify mutation in reference")
        row = rows.next()
        row = rows.next()
        self.assertEqual(round(row[0], 4),
                         1,
                         "SQLiteDB cannot correctly identify mutation in reference")
