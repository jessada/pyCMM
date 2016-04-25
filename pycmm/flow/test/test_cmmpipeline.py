import unittest
from os.path import join as join_path
from os.path import dirname
from pycmm.settings import DFLT_FLOW_ALLOC_TIME
from pycmm.settings import DFLT_RPT_ALLOC_TIME
from pycmm.template import SafeTester
from pycmm.flow import CMMPipeline
from pycmm.flow import create_jobs_setup_file

class TestCMMPipeline(SafeTester):

    def __init__(self, methodName):
        super(TestCMMPipeline, self).__init__(methodName=methodName,
                                              test_module_name=__name__,
                                              )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 project_name=None,
                                 project_out_dir=None,
                                 project_code=None,
                                 flow_alloc_time=None,
                                 rpt_alloc_time=None,
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if project_name is None:
            project_name = self.test_function
        if project_out_dir is None:
            project_out_dir = self.working_dir
        create_jobs_setup_file(project_name=project_name,
                               project_out_dir=project_out_dir,
                               project_code=project_code,
                               flow_alloc_time=flow_alloc_time,
                               rpt_alloc_time=rpt_alloc_time,
                               )
        return jobs_setup_file

    def test_load_jobs_info_1(self):
        """ test if default basic CMM pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = CMMPipeline(jobs_setup_file)
        self.assertEqual(pl.project_out_dir,
                         self.working_dir,
                         "CMMPipeline cannot correctly identify 'project output dir' from jobs setup file")
        self.assertEqual(pl.project_name,
                         self.current_func_name,
                         "CMMPipeline cannot correctly identify 'project name' from jobs setup file")
        self.assertEqual(pl.project_code,
                         None,
                         "CMMPipeline cannot correctly identify 'project code' from jobs setup file")
        
    def test_load_jobs_info_2(self):
        """ test if non-default basic CMM pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code="b2011097",
                                                        )
        pl = CMMPipeline(jobs_setup_file)
        self.assertEqual(pl.project_code,
                         "b2011097",
                         "CMMPipeline cannot correctly identify 'project code' from jobs setup file")
        self.assertEqual(pl.flow_alloc_time,
                         DFLT_FLOW_ALLOC_TIME,
                         "CMMPipeline cannot correctly identify 'flow allocation time' from jobs setup file")
        self.assertEqual(pl.rpt_alloc_time,
                         DFLT_RPT_ALLOC_TIME,
                         "CMMPipeline cannot correctly identify 'report allocation time' from jobs setup file")

    def test_load_jobs_info_3(self):
        """ test if special basic CMM pipeline configurations are loaded correctly """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code="b2012247",
                                                        flow_alloc_time="12:00:00",
                                                        rpt_alloc_time="10:00:00",
                                                        )
        pl = CMMPipeline(jobs_setup_file)
        self.assertEqual(pl.project_code,
                         "b2012247",
                         "CMMPipeline cannot correctly identify 'project code' from jobs setup file")
        self.assertEqual(pl.flow_alloc_time,
                         "12:00:00",
                         "CMMPipeline cannot correctly identify 'flow allocation time' from jobs setup file")
        self.assertEqual(pl.rpt_alloc_time,
                         "10:00:00",
                         "CMMPipeline cannot correctly identify 'report allocation time' from jobs setup file")

    def tearDown(self):
        self.remove_working_dir()
