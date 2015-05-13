import unittest
import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.utils.jobman import JobManager
from pycmm.utils.jobman import JOB_STATUS_PENDING

TEST_PROJECT_CODE = 'b2011097'
TEST_PARTITION_TYPE = 'core'
TEST_NTASKS = '1'
TEST_ALLOC_TIME = '10:00:00'

class TestJobManager(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'jobman'

    def __create_db_instance(self):
        return None

    def test_submit_job(self):
        """ check if a job has been correctly submitted to SLURM """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        job_params = 'hello testt_submit_job'
        jobman = JobManager()
        jobman.submit_job(job_name,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_file,
                          job_script,
                          job_params,
                          )
        job_id = jobman.get_job_id(job_name)
        self.assertTrue(job_id.isdigit(),
                        "Invalid number of queried records")
        jobman.cancel_job(job_name)

    def test_get_job_status(self):
        """ check if job status can be correctly retrieved """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        job_params = 'I want to test_get_job_status'
        jobman = JobManager()
        jobman.submit_job(job_name,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_file,
                          job_script,
                          job_params,
                          )
        job_status = jobman.get_job_status(job_name)
        self.assertEqual(job_status,
                         JOB_STATUS_PENDING,
                         "Invalid job status"
                         )
        jobman.cancel_job(job_name)

    def test_prerequisite(self):
        """ check if 'afterok' dependcy can be applied """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name_1 = self.test_function + "_1"
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        tmp_file = join_path(self.working_dir,
                             self.test_function+'.txt')
        job_params = '3 2000000 ' + tmp_file
        self.delete_file(tmp_file)
        jobman = JobManager()
        jobman.submit_job(job_name_1,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_file,
                          job_script,
                          job_params,
                          )
        job_name_2 = self.test_function + "_2"
        job_params = '100 1200000 ' + tmp_file
        jobman.submit_job(job_name_2,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_file,
                          job_script,
                          job_params,
                          prerequisite=[job_name_1],
                          )
        job_name_3 = self.test_function + "_3"
        job_params = '410 420 ' + tmp_file
        jobman.submit_job(job_name_3,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_file,
                          job_script,
                          job_params,
                          prerequisite=[job_name_1, job_name_2],
                          )
        jobman.cancel_job(job_name_1)
        jobman.cancel_job(job_name_2)
        jobman.cancel_job(job_name_3)

    def tearDown(self):
        self.remove_working_dir()
