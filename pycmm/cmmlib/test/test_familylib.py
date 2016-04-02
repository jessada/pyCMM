import unittest
from os.path import join as join_path
from pycmm.template import SafeTester
from pycmm.flow import CMMPipeline
from pycmm.flow import create_jobs_setup_file


class TestFamilyLib(SafeTester):

    def __init__(self, methodName):
        super(TestFamilyLib, self).__init__(methodName=methodName,
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
        self.assertEqual(pl.families_info["18"].members_info[0].sample_id,
                         "Co-345",
                         "familiylib can't handle families information")
        self.assertEqual(pl.families_info["911"].fam_id,
                         "911",
                         "familiylib can't handle families information")

    def test_extract_families_info_2(self):
        """ test if Family objects can be extract from no data """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(pl.families_info,
                         None,
                         "familiylib can't handle families information")

    def test_extract_samples_list_1(self):
        """ test if samples list can be extract from jobs setup data """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="18:Co-345:Co-37,12:Co-890:Co-290,13:Co-95,266:Co-131:Co-1355,314:1793-11o,987:Co-218:Co-2588,911:Co-1454:Co-4700,prostate:Pro001:Pro002:Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.samples_list),
                         15,
                         "familiylib can't handle families information")
        self.assertEqual('Co-2588' in pl.samples_list,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('Co-163' in pl.samples_list,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual('Pro003' in pl.samples_list,
                         True,
                         "familiylib can't handle families information")

    def test_extract_samples_list_2(self):
        """ test if samples list can be extract from no data """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(pl.samples_list,
                         None,
                         "familiylib can't handle families information")

    def test_extract_samples_list_3(self):
        """ test if samples list can be extract from jobs setup data without family id"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(sample_info="Co-345,Co-37,Co-890,Co-290,Co-95,Co-131,Co-1355,1793-11o,Co-218,Co-2588,Co-1454,Co-4700,Pro001,Pro002,Pro003")
        pl = CMMPipeline(jobs_setup_file=jobs_setup_file)
        self.assertEqual(len(pl.samples_list),
                         15,
                         "familiylib can't handle families information")
        self.assertEqual('Co-2588' in pl.samples_list,
                         True,
                         "familiylib can't handle families information")
        self.assertEqual('Co-163' in pl.samples_list,
                         False,
                         "familiylib can't handle families information")
        self.assertEqual('Pro003' in pl.samples_list,
                         True,
                         "familiylib can't handle families information")

    def tearDown(self):
        self.remove_working_dir()
