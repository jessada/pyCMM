import unittest
from os.path import join as join_path
#from os.path import dirname
from pycmm.template import SafeTester
from pycmm.settings import MAX_REF_MAF_COL_NAME
from pycmm.proc.mutrep.dbreader import SQLiteDBReader
from pycmm.proc.mutrep.dbreader import GT_WT
from pycmm.proc.mutrep.dbreader import GT_HOM
from pycmm.proc.mutrep.dbreader import GT_HET
from pycmm.proc.mutrep.mutrep import MutRepController
from pycmm.proc.db.dbms import create_dbms_jobs_setup_file as create_jobs_setup_file
from pycmm.proc.db.dbms import SQLiteDBController
from pycmm.proc.db.connector import TBL_NAME_ALL_GTZ_ANNOS
#from pycmm.settings import KG2014OCT_ALL_COL_NAME
#from pycmm.settings import EST_KVOT_EARLYONSET_VS_BRC_COL_NAME
#from pycmm.settings import EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME
#from pycmm.settings import EST_KVOT_EARLYONSET_VS_KG_EUR_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_SYN_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_SYN_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_SYN_Z_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_MIS_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_MIS_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_MIS_Z_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_EXP_LOF_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_N_LOF_COL_NAME
from pycmm.settings import EXAC03_CONSTRAINT_PLI_COL_NAME
from pycmm.settings import INTERVAR_CLASS_COL_NAME
from pycmm.settings import INTERVAR_EVIDENCE_COL_NAME
from pycmm.settings import PATHOGENIC_COUNT_COL_NAME
#from pycmm.settings import MAX_REF_MAF_COL_NAME
#from pycmm.settings import WES294_OAF_EARLYONSET_AF_COL_NAME
#from pycmm.settings import WES294_OAF_BRCS_AF_COL_NAME
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_BENIGN
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_LIKELY_BENIGN
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_LIKELY_PATHOGENIC
from pycmm.cmmlib.intervarlib import INTERVAR_CLASS_PATHOGENIC


class TestQryCall(SafeTester):

    def __init__(self, methodName):
        super(TestQryCall, self).__init__(methodName=methodName,
                                            test_module_name=__name__,
                                            )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self, *args, **kwargs):
        if 'project_name' not in kwargs:
            kwargs['project_name'] = self.test_function
        if 'db_file' not in kwargs:
            kwargs['db_file'] = join_path(self.data_dir, "input.db")
        kwargs['project_out_dir'] = self.working_dir
        kwargs['jobs_setup_file'] = join_path(self.working_dir,
                                              kwargs['project_name']+'_jobs_setup.txt')
        create_jobs_setup_file(*args, **kwargs)
        return kwargs['jobs_setup_file']

    def test_raw_zygo_1(self):
        """ test if raw_zygosity can be correctly parsed """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records(chrom="6")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].raw_zygo,
                         GT_HOM,
                         "raw zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].raw_zygo,
                         GT_HET,
                         "raw zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].raw_zygo,
                         GT_WT,
                         "raw zygosity cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].raw_zygo,
                         GT_HOM,
                         "raw zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].raw_zygo,
                         GT_HET,
                         "raw zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].raw_zygo,
                         GT_WT,
                         "raw zygosity cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].raw_zygo,
                         GT_WT,
                         "raw zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].raw_zygo,
                         GT_WT,
                         "raw zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].raw_zygo,
                         GT_WT,
                         "raw zygosity cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].raw_zygo,
                         GT_HOM,
                         "raw zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].raw_zygo,
                         GT_HET,
                         "raw zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].raw_zygo,
                         GT_WT,
                         "raw zygosity cannot be correctly determined")

    def test_actual_zygo_1(self):
        """
        test if true zygosity can be correctly determined
        - very random cases
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records(chrom="6")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].actual_zygo,
                         GT_WT,
                         "actual zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].actual_zygo,
                         GT_HET,
                         "actual zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].actual_zygo,
                         GT_HOM,
                         "actual zygosity cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].actual_zygo,
                         GT_HOM,
                         "actual zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].actual_zygo,
                         GT_HET,
                         "actual zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].actual_zygo,
                         GT_WT,
                         "actual zygosity cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].actual_zygo,
                         GT_WT,
                         "actual zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].actual_zygo,
                         GT_WT,
                         "actual zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].actual_zygo,
                         GT_WT,
                         "actual zygosity cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].actual_zygo,
                         GT_HOM,
                         "actual zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].actual_zygo,
                         GT_HET,
                         "actual zygosity cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].actual_zygo,
                         GT_WT,
                         "actual zygosity cannot be correctly determined")

    def test_mutated_1(self):
        """
        test if mutation can be identified
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records(chrom="6")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].mutated,
                         False,
                         "mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].mutated,
                         True,
                         "mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].mutated,
                         True,
                         "mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].mutated,
                         True,
                         "mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].mutated,
                         True,
                         "mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].mutated,
                         False,
                         "mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].mutated,
                         False,
                         "mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].mutated,
                         False,
                         "mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].mutated,
                         False,
                         "mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["AB-14"].mutated,
                         True,
                         "mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1141"].mutated,
                         True,
                         "mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["5053-10D"].mutated,
                         False,
                         "mutation cannot be correctly determined")

    def test_shared_mutation_1(self):
        """
        test if shared mutation in a family can be identified
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        dbr = SQLiteDBReader(db_file, family_infos=dbc.families_info, verbose=False)
        self.assertEqual(dbr.count_rows(TBL_NAME_ALL_GTZ_ANNOS),
                         18,
                         "SQLiteDBController cannot correctly execute db jobs")
        qry_records = dbr.get_qry_records()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["813-06o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["91-04o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-37"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1467"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1681"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1116"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-358"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-364"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["813-06o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["91-04o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-37"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1467"].shared_mutation,
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1681"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1116"].shared_mutation,
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-358"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-364"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["813-06o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["91-04o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-37"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1467"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1681"].shared_mutation,
                         True,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1116"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-358"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-364"].shared_mutation,
                         True,
                         "shared mutation cannot be correctly determined")
        # Test w/o families_info
        dbr = SQLiteDBReader(db_file, verbose=False)
        qry_records = dbr.get_qry_records()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["813-06o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["91-04o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-37"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1467"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1681"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1116"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-358"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-364"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["813-06o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["91-04o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-37"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1467"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1681"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1116"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-358"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-364"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.gt["813-06o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["91-04o"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-37"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1467"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1681"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-1116"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-358"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")
        self.assertEqual(qry_record.gt["Co-364"].shared_mutation,
                         False,
                         "shared mutation cannot be correctly determined")

class TestSQLiteDBReader(SafeTester):

    def __init__(self, methodName):
        super(TestSQLiteDBReader, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self, *args, **kwargs):
        if 'project_name' not in kwargs:
            kwargs['project_name'] = self.test_function
        if 'db_file' not in kwargs:
            kwargs['db_file'] = join_path(self.data_dir, "input.db")
        kwargs['project_out_dir'] = self.working_dir
        kwargs['jobs_setup_file'] = join_path(self.working_dir,
                                              self.test_function+'_jobs_setup.txt')
        create_jobs_setup_file(*args, **kwargs)
        return kwargs['jobs_setup_file']

    def test_filter_non_intergenic_1(self):
        """ test counting non-intergenic variants """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        db_reader.filter_non_intergenic = True
        qry_records = db_reader.get_qry_records()
        self.assertEqual(len(list(qry_records)),
                         31,
                         "non-intergenic variants cannot be correctly determined")

    def test_filter_non_intronic_1(self):
        """ test counting non-intronic variants """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        db_reader.filter_non_intronic = True
        qry_records = db_reader.get_qry_records()
        self.assertEqual(len(list(qry_records)),
                         11,
                         "non-intronic variants cannot be correctly determined")

    def test_filter_non_upstream_1(self):
        """ test counting non-upstream variants """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        db_reader.filter_non_upstream = True
        qry_records = db_reader.get_qry_records()
        self.assertEqual(len(list(qry_records)),
                         27,
                         "non-upstream variants cannot be correctly determined")

    def test_filter_non_downstream_1(self):
        """ test counting non-downstream variants """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        db_reader.filter_non_downstream = True
        qry_records = db_reader.get_qry_records()
        self.assertEqual(len(list(qry_records)),
                         15,
                         "non-downstream variants cannot be correctly determined")

    def test_filter_non_utr_1(self):
        """ test counting non-utr variants """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        db_reader.filter_non_utr = True
        qry_records = db_reader.get_qry_records()
        self.assertEqual(len(list(qry_records)),
                         11,
                         "non-utr variants cannot be correctly determined")

    def test_filter_non_synonymous_1(self):
        """ test counting non-synonymous variants """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        db_reader.filter_non_synonymous = True
        qry_records = db_reader.get_qry_records()
        self.assertEqual(len(list(qry_records)),
                         11,
                         "non-synonymous variants cannot be correctly determined")

class TestQryRecord(SafeTester):

    def __init__(self, methodName):
        super(TestQryRecord, self).__init__(methodName=methodName,
                                            test_module_name=__name__,
                                            )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self, *args, **kwargs):
        if 'project_name' not in kwargs:
            kwargs['project_name'] = self.test_function
        if 'db_file' not in kwargs:
            kwargs['db_file'] = join_path(self.data_dir, "input.db")
        kwargs['project_out_dir'] = self.working_dir
        kwargs['jobs_setup_file'] = join_path(self.working_dir,
                                              self.test_function+'_jobs_setup.txt')
        create_jobs_setup_file(*args, **kwargs)
        return kwargs['jobs_setup_file']

    def test_has_mutation_1(self):
        """ test a mutation can be correctly identified in one sample """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        samples = ["Co-1141"]
        qry_records = db_reader.get_qry_records()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         True,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         True,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         True,
                         "mutations in a sample cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in a sample cannot be correctly determined")

    def test_has_mutation_2(self):
        """ test a mutation can be correctly identified in two samples """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        samples = ["Co-454", "Co-669"]
        qry_records = db_reader.get_qry_records()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         True,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         True,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         True,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_mutation(samples),
                         False,
                         "mutations in samples cannot be correctly determined")

    def test_has_shared_1(self):
        """ test if a family with shared mutation can be detected """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        dbc= SQLiteDBController(jobs_setup_file, verbose=False)
        db_file = join_path(self.data_dir,
                            'input.db')
        dbr = SQLiteDBReader(db_file, family_infos=dbc.families_info, verbose=False)
        qry_records = dbr.get_qry_records()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         False,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.has_shared(),
                         True,
                         "shared mutation cannot be correctly determined")

    def test_qry_eval_1(self):
        """ test general evaluation of qry_eval """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        expr1 = '("dpsi_zscore"!=\'\')and(float("dpsi_zscore")>-2)and("ExonicFunc_refGene"==\'synonymous_SNV\')and(float("dpsi_zscore")<2)'
        expr2 = '(float("1000g2014oct_all")>0.01)'
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.qry_eval(expr1),
                         True,
                         "cannot perform qry expression evaluation correctly")
        self.assertEqual(qry_record.qry_eval(expr2),
                         True,
                         "cannot perform qry expression evaluation correctly")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.qry_eval(expr1),
                         False,
                         "cannot perform qry expression evaluation correctly")
        self.assertEqual(qry_record.qry_eval(expr2),
                         True,
                         "cannot perform qry expression evaluation correctly")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.qry_eval(expr1),
                         True,
                         "cannot perform qry expression evaluation correctly")
        self.assertEqual(qry_record.qry_eval(expr2),
                         False,
                         "cannot perform qry expression evaluation correctly")

    def test_max_ref_maf_1(self):
        """
        test finding pre-calculated maximum reference minor allele frequency
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records(chrom='6')
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(MAX_REF_MAF_COL_NAME),
                         0.235,
                         "values of MAX_REF_MAF cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(MAX_REF_MAF_COL_NAME),
                         0.4963,
                         "values of MAX_REF_MAF cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(MAX_REF_MAF_COL_NAME),
                         0.0155,
                         "values of MAX_REF_MAF cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(MAX_REF_MAF_COL_NAME),
                         0.4969,
                         "values of MAX_REF_MAF cannot be correctly determined")

    def test_max_ref_maf_2(self):
        """
        test finding on-demand maximum reference minor allele frequency
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records(chrom='6')
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(MAX_REF_MAF_COL_NAME),
                         0.235,
                         "values of MAX_REF_MAF cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(MAX_REF_MAF_COL_NAME),
                         0.4963,
                         "values of MAX_REF_MAF cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(MAX_REF_MAF_COL_NAME),
                         0.0155,
                         "values of MAX_REF_MAF cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(MAX_REF_MAF_COL_NAME),
                         0.4969,
                         "values of MAX_REF_MAF cannot be correctly determined")

    def test_cal_est_kvot_1(self):
        """
        test finding maximum reference minor allele frequency
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records()
#    def test_qry_eval_2(self):
#        """
#        test qry_eval with multi alleleic INFO
#        """
#
#        self.init_test(self.current_func_name)
#        db_file = join_path(self.data_dir,
#                            'input.vcf.gz')
#        expr1 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" == \'NA\')'
#        expr2 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" == \'INF\')'
#        expr3 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" == \'\')'
#        expr4 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'NA\')'
#        expr4 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'INF\')'
#        expr4 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'\')'
#        expr4 += ' and (float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '") < 0.5)'
#        expr5 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'NA\')'
#        expr5 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'INF\')'
#        expr5 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'\')'
#        expr5 += ' and (float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '") < 0.8)'
#        expr6 = '("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'NA\')'
#        expr6 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'INF\')'
#        expr6 += ' and ("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '" != \'\')'
#        expr6 += ' and (float("' + EST_KVOT_EARLYONSET_VS_BRC_COL_NAME + '") < 2)'
#        expr7 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" == \'NA\')'
#        expr8 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" == \'INF\')'
#        expr9 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" == \'\')'
#        expr10 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'NA\')'
#        expr10 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'INF\')'
#        expr10 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'\')'
#        expr10 += ' and (float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '") < 0.5)'
#        expr11 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'NA\')'
#        expr11 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'INF\')'
#        expr11 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'\')'
#        expr11 += ' and (float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '") < 0.8)'
#        expr12 = '("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'NA\')'
#        expr12 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'INF\')'
#        expr12 += ' and ("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '" != \'\')'
#        expr12 += ' and (float("' + EST_KVOT_EARLYONSET_VS_EXAC_NFE_COL_NAME + '") < 2)'
#        db_reader = SQLiteDBReader(db_file)
#        qry_record = db_reader.next()
#        qry_record = db_reader.next()
#        # allele idx equal to 1 mean the first alternate allele
#        # allele idx equal to 2 mean the second alternate allele
#        # so on
#        self.assertTrue(qry_record.qry_eval(expr1, 1),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr2, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr3, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr4, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr5, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr6, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        qry_record = db_reader.next()
#        qry_record = db_reader.next()
#        qry_record = db_reader.next()
#        self.assertFalse(qry_record.qry_eval(expr7, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr8, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr9, 1),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr10, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr11, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr12, 1),
#                         "cannot perform vcf expression evaluation correctly")
#        qry_record = db_reader.next()
#        self.assertFalse(qry_record.qry_eval(expr1, 2),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr2, 2),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr3, 2),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr4, 2),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr5, 2),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr6, 2),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr1, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr2, 3),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr3, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr4, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr5, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr6, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr7, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr8, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr9, 3),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr10, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr11, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr12, 3),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr1, 4),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr2, 4),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr3, 4),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr4, 4),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr5, 4),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr6, 4),
#                        "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr1, 5),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr2, 5),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr3, 5),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr4, 5),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertFalse(qry_record.qry_eval(expr5, 5),
#                         "cannot perform vcf expression evaluation correctly")
#        self.assertTrue(qry_record.qry_eval(expr6, 5),
#                        "cannot perform vcf expression evaluation correctly")
#
    def test_pathogenic_count_1(self):
        """
        test if harmful pathogenic variants can be correctly count
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         1,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         9,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         1,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         1,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         0,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         9,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         0,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         5,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         9,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         3,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         5,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         4,
                         "number of harmful pathogenic prediction cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(PATHOGENIC_COUNT_COL_NAME),
                         9,
                         "number of harmful pathogenic prediction cannot be correctly determined")

    def test_parse_exac_constraint_1(self):
        """
        test parsing exac constraint from the annotated vcf
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records()
        qry_record = qry_records.next()
        qry_record = qry_records.next()
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_EXP_SYN_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_N_SYN_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_SYN_Z_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_EXP_MIS_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_N_MIS_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_MIS_Z_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_EXP_LOF_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_N_LOF_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_PLI_COL_NAME) == '',
                        "values of ExAC constraint cannot be correctly determined")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_EXP_SYN_COL_NAME),
                         '45.7717977506',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_N_SYN_COL_NAME),
                         '40',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_SYN_Z_COL_NAME),
                         '0.528885612026385',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_EXP_MIS_COL_NAME),
                         '115.696060891',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_N_MIS_COL_NAME),
                         '34',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_MIS_Z_COL_NAME),
                         '3.71499650338098',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_EXP_LOF_COL_NAME),
                         '15.533928734',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_N_LOF_COL_NAME),
                         '1',
                         "values of ExAC constraint cannot be correctly determined")
        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_PLI_COL_NAME),
                         '0.975506865848027',
                         "values of ExAC constraint cannot be correctly determined")

#    def test_parse_exac_constraint_2(self):
#        """
#        test parsing exac constraint for vcf without the annotation
#        """
#
#        self.init_test(self.current_func_name)
#        db_file = join_path(self.data_dir,
#                            'input.vcf.gz')
#        db_reader = SQLiteDBReader(db_file)
#        qry_record = db_reader.next()
#        self.assertTrue(qry_record.get_anno(EXAC03_CONSTRAINT_EXP_SYN_COL_NAME) == '',
#                        "values of ExAC constraint cannot be correctly determined")
#
#    def test_parse_exac_constraint_4(self):
#        """
#        test parsing exac constraint for vcf without the annotation
#        """
#
#        self.init_test(self.current_func_name)
#        db_file = join_path(self.data_dir,
#                            'input.vcf.gz')
#        db_reader = SQLiteDBReader(db_file)
#        qry_record = db_reader.next()
#        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_N_LOF_COL_NAME, 1),
#                         '',
#                         "values of ExAC constraint cannot be correctly determined")
#        self.assertEqual(qry_record.get_anno(EXAC03_CONSTRAINT_PLI_COL_NAME, 2),
#                         '',
#                         "values of ExAC constraint cannot be correctly determined")
#
    def test_parse_intervar_1(self):
        """
        test parsing intervar
        """

        self.init_test(self.current_func_name)
        db_file = join_path(self.data_dir,
                            'input.db')
        db_reader = SQLiteDBReader(db_file, verbose=False)
        qry_records = db_reader.get_qry_records()
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(INTERVAR_CLASS_COL_NAME),
                         INTERVAR_CLASS_BENIGN,
                         "Incorect intervar value")
        self.assertEqual(qry_record.get_anno(INTERVAR_EVIDENCE_COL_NAME),
                         "BA1, BS1, BP4, BP7",
                         "Incorect intervar value")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(INTERVAR_CLASS_COL_NAME),
                         INTERVAR_CLASS_LIKELY_BENIGN,
                         "Incorect intervar value")
        self.assertEqual(qry_record.get_anno(INTERVAR_EVIDENCE_COL_NAME),
                         "PM1, BS2, BP4",
                         "Incorect intervar value")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(INTERVAR_CLASS_COL_NAME),
                         INTERVAR_CLASS_UNCERTAIN_SIGNIFICANCE,
                         "Incorect intervar value")
        self.assertEqual(qry_record.get_anno(INTERVAR_EVIDENCE_COL_NAME),
                         "",
                         "Incorect intervar value")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(INTERVAR_CLASS_COL_NAME),
                         INTERVAR_CLASS_BENIGN,
                         "Incorect intervar value")
        self.assertEqual(qry_record.get_anno(INTERVAR_EVIDENCE_COL_NAME),
                         "BS1, BS2",
                         "Incorect intervar value")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(INTERVAR_CLASS_COL_NAME),
                         INTERVAR_CLASS_PATHOGENIC,
                         "Incorect intervar value")
        self.assertEqual(qry_record.get_anno(INTERVAR_EVIDENCE_COL_NAME),
                         "PVS1, PM2, PP5",
                         "Incorect intervar value")
        qry_record = qry_records.next()
        self.assertEqual(qry_record.get_anno(INTERVAR_CLASS_COL_NAME),
                         INTERVAR_CLASS_LIKELY_PATHOGENIC,
                         "Incorect intervar value")
        self.assertEqual(qry_record.get_anno(INTERVAR_EVIDENCE_COL_NAME),
                         "PVS1, PM2",
                         "Incorect intervar value")
