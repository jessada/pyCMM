import unittest
from os.path import join as join_path
from pycmm.settings import PHENOTYPE_MISSING
from pycmm.settings import PHENOTYPE_UNAFFECTED
from pycmm.settings import PHENOTYPE_AFFECTED
from pycmm.template import SafeTester
from pycmm.flow import CMMPipeline
from pycmm.flow import create_jobs_setup_file
from pycmm.cmmlib.samplelib import NO_FAMILY
from pycmm.cmmlib.samplelib import SAMPLE_GROUP_MISSING


class TestSample(SafeTester):

    def __init__(self, methodName):
        super(TestSample, self).__init__(methodName=methodName,
                                         test_module_name=__name__,
                                         )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 sample_info=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        project_name = self.test_function
        project_out_dir = self.working_dir
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=project_out_dir,
                               sample_info=sample_info,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    def test_sample_group_1(self):
        """ test if sample group can be determined if there is no information """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="18:Co-345:Co-37,12:Co-890:Co-290,13:Co-95,266:Co-131:Co-1355,314:1793-11o,987:Co-218:Co-2588,911:Co-1454:Co-4700,prostate:Pro001:Pro002:Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.families_info),
                         8,
                         "familiylib can't handle families information")
        self.assertEqual('prostate' in pl.families_info,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('breast' in pl.families_info,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["18"].members[0].sample_id,
                         "Co-345",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["18"].members[0].sample_group,
                         SAMPLE_GROUP_MISSING,
                         "sample group cannot be determined")
        self.assertEqual(pl.families_info["911"].fam_id,
                         "911",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["prostate"].members[0].fam_id,
                         "prostate",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["prostate"].members[0].sample_group,
                         SAMPLE_GROUP_MISSING,
                         "sample group cannot be determined")

    def test_sample_group_2(self):
        """ test backward compatible if loading from file (no sample group) """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual('breast' in pl.families_info,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["87"].members[0].sample_id,
                         "Co-218",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["87"].members[0].sample_group,
                         SAMPLE_GROUP_MISSING,
                         "sample group cannot be determined")
        self.assertEqual(pl.families_info["8"].members[1].sample_group,
                         SAMPLE_GROUP_MISSING,
                         "sample group cannot be determined")

    def test_sample_group_3(self):
        """ test group can be loaded if they are identified """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual('breast' in pl.families_info,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["13"].members[1].sample_id,
                         "Co-94",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["13"].members[1].sample_group,
                         SAMPLE_GROUP_MISSING,
                         "sample group cannot be determined")
        self.assertEqual(pl.families_info["87"].members[0].sample_id,
                         "Co-218",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["87"].members[0].sample_group,
                         3,
                         "sample group cannot be determined")
        self.assertEqual(pl.families_info["8"].members[1].sample_group,
                         3,
                         "sample group cannot be determined")
        self.assertEqual(pl.families_info[NO_FAMILY].members[3].sample_id,
                         "1199-05o",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info[NO_FAMILY].members[3].sample_group,
                         1,
                         "sample group cannot be determined")
        self.assertEqual(pl.families_info[NO_FAMILY].members[4].sample_group,
                         1,
                         "sample group cannot be determined")

    def test_sample_phenotype_1(self):
        """ test if sample phenotype can be determined if there is no information """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="18:Co-345:Co-37,12:Co-890:Co-290,13:Co-95,266:Co-131:Co-1355,314:1793-11o,987:Co-218:Co-2588,911:Co-1454:Co-4700,prostate:Pro001:Pro002:Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.families_info),
                         8,
                         "familiylib can't handle families information")
        self.assertEqual('prostate' in pl.families_info,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('breast' in pl.families_info,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["18"].members[0].sample_id,
                         "Co-345",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["18"].members[0].phenotype,
                         PHENOTYPE_MISSING,
                         "sample phenotype cannot be determined")
        self.assertEqual(pl.families_info["911"].fam_id,
                         "911",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["prostate"].members[0].fam_id,
                         "prostate",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["prostate"].members[0].phenotype,
                         SAMPLE_GROUP_MISSING,
                         "sample phenotype cannot be determined")

    def test_sample_phenotype_2(self):
        """ test backward compatible if loading from file (no phenotype) """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual('breast' in pl.families_info,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["87"].members[0].sample_id,
                         "Co-218",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["87"].members[0].phenotype,
                         SAMPLE_GROUP_MISSING,
                         "sample phenotype cannot be determined")
        self.assertEqual(pl.families_info["8"].members[1].phenotype,
                         SAMPLE_GROUP_MISSING,
                         "sample phenotype cannot be determined")

    def test_sample_phenotype_3(self):
        """ test phenotype can be loaded if they are identified """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual('breast' in pl.families_info,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["13"].members[1].sample_id,
                         "Co-94",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["13"].members[1].phenotype,
                         PHENOTYPE_MISSING,
                         "sample phenotype cannot be determined")
        self.assertEqual(pl.families_info["87"].members[0].sample_id,
                         "Co-218",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["87"].members[0].phenotype,
                         PHENOTYPE_MISSING,
                         "sample phenotype cannot be determined")
        self.assertEqual(pl.families_info["87"].members[1].phenotype,
                         PHENOTYPE_UNAFFECTED,
                         "sample phenotype cannot be determined")
        self.assertEqual(pl.families_info[NO_FAMILY].members[3].sample_id,
                         "1199-05o",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info[NO_FAMILY].members[3].phenotype,
                         PHENOTYPE_AFFECTED,
                         "sample phenotype cannot be determined")
        self.assertEqual(pl.families_info[NO_FAMILY].members[4].phenotype,
                         PHENOTYPE_AFFECTED,
                         "sample phenotype cannot be determined")

    def tearDown(self):
        self.remove_working_dir()

class TestFamily(SafeTester):

    def __init__(self, methodName):
        super(TestFamily, self).__init__(methodName=methodName,
                                         test_module_name=__name__,
                                         )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 sample_info=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        project_name = self.test_function
        project_out_dir = self.working_dir
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=project_out_dir,
                               sample_info=sample_info,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    def test_extract_families_info_1(self):
        """ test if Family objects can be extract from jobs setup data """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="18:Co-345:Co-37,12:Co-890:Co-290,13:Co-95,266:Co-131:Co-1355,314:1793-11o,987:Co-218:Co-2588,911:Co-1454:Co-4700,prostate:Pro001:Pro002:Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.families_info),
                         8,
                         "familiylib can't handle families information")
        self.assertEqual('prostate' in pl.families_info,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('breast' in pl.families_info,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["18"].members[0].sample_id,
                         "Co-345",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["911"].fam_id,
                         "911",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["prostate"].members[0].fam_id,
                         "prostate",
                         "familiylib can't handle families information")

    def test_extract_families_info_2(self):
        """ test if Family objects can be extract from no data """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(pl.families_info,
                         None,
                         "familiylib can't handle families information")

    def tearDown(self):
        self.remove_working_dir()

class TestSamplesInfo(SafeTester):

    def __init__(self, methodName):
        super(TestSamplesInfo, self).__init__(methodName=methodName,
                                              test_module_name=__name__,
                                              )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 sample_info=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        project_name = self.test_function
        project_out_dir = self.working_dir
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=project_out_dir,
                               sample_info=sample_info,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

    def test_has_info_1(self):
        """ test if SamplesInfo has no information"""

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertFalse(pl.has_samples_info,
                         "SamplesInfo cannot tell correctly if it has information")

    def test_has_info_2(self):
        """ test if SamplesInfo has information"""

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="18:Co-345:Co-37,12:Co-890:Co-290,13:Co-95,266:Co-131:Co-1355,314:1793-11o,987:Co-218:Co-2588,911:Co-1454:Co-4700,prostate:Pro001:Pro002:Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertTrue(pl.has_samples_info,
                        "SamplesInfo cannot tell correctly if it has information")

    def test_extract_samples_id_1(self):
        """ test if samples list can be extract from jobs setup data """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="18:Co-345:Co-37,12:Co-890:Co-290,13:Co-95,266:Co-131:Co-1355,314:1793-11o,987:Co-218:Co-2588,911:Co-1454:Co-4700,prostate:Pro001:Pro002:Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.samples_id),
                         15,
                         "familiylib can't handle families information")
        self.assertEqual('Co-2588' in pl.samples_id,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('Co-163' in pl.samples_id,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual('Pro003' in pl.samples_id,
                         True,
                         "familiylib can't handle families information")

    def test_extract_samples_id_2(self):
        """ test if samples list can be extract from no data """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(pl.samples_id,
                         None,
                         "familiylib can't handle families information")

    def test_extract_samples_id_3(self):
        """ test if samples list can be extract from jobs setup data without family id"""

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="Co-345,Co-37,Co-890,Co-290,Co-95,Co-131,Co-1355,1793-11o,Co-218,Co-2588,Co-1454,Co-4700,Pro001,Pro002,Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.samples_id),
                         15,
                         "familiylib can't handle families information")
        self.assertEqual('Co-2588' in pl.samples_id,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('Co-163' in pl.samples_id,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual('Pro003' in pl.samples_id,
                         True,
                         "familiylib can't handle families information")

    def test_extract_samples_id_w_fam_pref_1(self):
        """ 
        test if samples id can be extract from jobs setup data """

        self.init_test(self.current_func_name)
        sample_info_list = []
        sample_info_list.append("18:Co-345:Co-37")
        sample_info_list.append("12:Co-890:Co-290")
        sample_info_list.append("13:Co-95")
        sample_info_list.append("266:Co-131:Co-1355")
        sample_info_list.append("314:1793-11o")
        sample_info_list.append("987:Co-218:Co-2588")
        sample_info_list.append("prostate:Pro001:Pro002:Pro003")
        sample_info_list.append(NO_FAMILY+":AB-27:AB-28:1003-06o:1068-05o")
        sample_info_list.append(NO_FAMILY+"_232:789-04o:AB-33:Co-1258")
        sample_info = ",".join(sample_info_list)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.samples_id_w_fam_pref),
                         20,
                         "familiylib can't handle families information")
        self.assertEqual('987-Co-2588' in pl.samples_id_w_fam_pref,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('Co-2588' in pl.samples_id_w_fam_pref,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual('prostate-Pro003' in pl.samples_id_w_fam_pref,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('AB-27' in pl.samples_id_w_fam_pref,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('1068-05o' in pl.samples_id_w_fam_pref,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('789-04o' in pl.samples_id_w_fam_pref,
                         True,
                         "familiylib can't handle families information")

#    def test_samples_phenotype_1(self):
#        """ test if sample phenotype can be determined if there is no information """
#
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(sample_info="18:Co-345:Co-37,12:Co-890:Co-290,13:Co-95,266:Co-131:Co-1355,314:1793-11o,987:Co-218:Co-2588,911:Co-1454:Co-4700,prostate:Pro001:Pro002:Pro003")
#        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
#        self.assertEqual(len(pl.affected_samples),
#                         0,
#                         "incorrect number of affected samples")
#
#    def test_samples_phenotype_2(self):
#        """ test backward compatible if loading from file (no phenotype) """
#
#        self.init_test(self.current_func_name)
#        sample_info = join_path(self.data_dir,
#                                "sample.info")
#        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
#        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
#        self.assertEqual(len(pl.affected_samples),
#                         0,
#                         "incorrect number of affected samples")
#
#    def test_samples_phenotype_3(self):
#        """ test sample phenotype can be loaded if they are identified """
#
#        self.init_test(self.current_func_name)
#        sample_info = join_path(self.data_dir,
#                                "sample.info")
#        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
#        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
#        self.assertEqual(len(pl.affected_samples),
#                         8,
#                         "incorrect number of affected samples")
#        self.assertEqual(len(pl.unaffected_samples),
#                         7,
#                         "incorrect number of unaffected samples")
#        self.assertEqual(len(pl.samples_groups[1].affected),
#                         8,
#                         "incorrect number of affected samples")
#        self.assertEqual(len(pl.samples_groups[1].unaffected),
#                         0,
#                         "incorrect number of unaffected samples")
#        self.assertEqual(len(pl.samples_groups[2].affected),
#                         0,
#                         "incorrect number of affected samples")
#        self.assertEqual(len(pl.samples_groups[2].unaffected),
#                         2,
#                         "incorrect number of unaffected samples")
#        self.assertEqual(len(pl.samples_groups[3].affected),
#                         0,
#                         "incorrect number of affected samples")
#        self.assertEqual(len(pl.samples_groups[3].unaffected),
#                         2,
#                         "incorrect number of unaffected samples")
#        self.assertEqual(len(pl.samples_groups[4].affected),
#                         0,
#                         "incorrect number of affected samples")
#        self.assertEqual(len(pl.samples_groups[4].unaffected),
#                         2,
#                         "incorrect number of unaffected samples")
#
    def test_samples_groups_1(self):
        """ test if sample group; can be determined if there is no information """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="18:Co-345:Co-37,12:Co-890:Co-290,13:Co-95,266:Co-131:Co-1355,314:1793-11o,987:Co-218:Co-2588,911:Co-1454:Co-4700,prostate:Pro001:Pro002:Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.samples_groups[3]),
                         0,
                         "incorrect number of affected samples")

    def test_samples_groups_2(self):
        """ test backward compatible if loading from file (no group) """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        pl = CMMPipeline(jobs_setup_file)
        self.assertEqual(len(pl.samples_groups[3]),
                         0,
                         "incorrect number of affected samples")

    def test_samples_groups_3(self):
        """ test sample group can be loaded if they are identified """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        jobs_setup_file = self.__create_jobs_setup_file(sample_info=sample_info)
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.samples_groups[1]),
                         8,
                         "incorrect number of affected samples")
        self.assertEqual(len(pl.samples_groups[2]),
                         4,
                         "incorrect number of affected samples")
        self.assertEqual(len(pl.samples_groups[3]),
                         6,
                         "incorrect number of affected samples")
        self.assertEqual(len(pl.samples_groups[4]),
                         4,
                         "incorrect number of affected samples")

    def tearDown(self):
        self.remove_working_dir()
