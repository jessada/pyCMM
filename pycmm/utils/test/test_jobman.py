import unittest
import filecmp
from os.path import join as join_path
from os.path import dirname
from pycmm.settings import SLURM_TEST
from pycmm.template import SafeTester
from pycmm.utils import exec_sh
from pycmm.utils.jobman import JobManager
from pycmm.utils.jobman import JOB_STATUS_PENDING
from pycmm.utils.jobman import JOB_STATUS_RUNNING

TEST_PROJECT_CODE = 'b2011097'
TEST_PARTITION_TYPE = 'core'
TEST_NTASKS = '1'
TEST_ALLOC_TIME = '10:00'

class JobManagerHouseKeeping(JobManager):
    """ A test class to test garbage collecting of JobManager """

    def __init__(self,
                 jobs_report_file=None,
                 ):
        JobManager.__init__(self,
                            jobs_report_file,
                            )

    def get_raw_repr(self):
        return {"NA1": "NA1",
                "NA2": "NA2",
                }

    def monitor_action(self):
        JobManager.monitor_action(self)


class TestJobManager(SafeTester):

    def __init__(self, methodName):
        super(TestJobManager, self).__init__(methodName=methodName,
                                             test_module_name=__name__,
                                             )

    def setUp(self):
        pass

    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
    def test_submit_job_1(self):
        """ check if a job has been correctly submitted to SLURM """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_prefix = join_path(self.working_dir,
                                     self.test_function)
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        tmp_file = join_path(self.working_dir,
                             self.test_function+'.txt')
        job_params = '5 100000 ' + tmp_file
        jobman = JobManager()
        jobman.submit_job(job_name,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          job_script,
                          job_params,
                          )
        job_id = jobman.get_job_id(job_name)
        self.assertTrue(job_id.isdigit(),
                        "Invalid number of queried records")
        jobman.cancel_job(job_name)

    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
    def test_get_job_status(self):
        """ check if job status can be correctly retrieved """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_prefix = join_path(self.working_dir,
                                     self.test_function)
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        tmp_file = join_path(self.working_dir,
                             self.test_function+'.txt')
        job_params = '3 200000 ' + tmp_file
        jobman = JobManager()
        jobman.submit_job(job_name,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          job_script,
                          job_params,
                          )
        job_status = jobman.get_job_status(job_name)
        self.assertTrue((job_status == JOB_STATUS_PENDING) or
                        (job_status == JOB_STATUS_RUNNING),
                        "Invalid job status"
                        )
        jobman.cancel_job(job_name)

    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
    def test_get_job_usage_mail(self):
        """
        - check if SLURM usage_mail can be applied with the flow
        NOTE -> need to manually check the result
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_prefix = join_path(self.working_dir,
                                     self.test_function)
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        tmp_file = join_path(self.working_dir,
                             self.test_function+'.txt')
        job_params = '3 20000 ' + tmp_file
        self.delete_file(tmp_file)
        jobman = JobManager()
        jobman.submit_job(job_name,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          job_script,
                          job_params,
                          email=True,
                          )

    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
    def test_monitor_jobs(self):
        """
        - check JobManager can correctly monitor all the jobs
        NOTE -> need to manually check the result
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name_1 = self.test_function + "_1"
        slurm_log_prefix = join_path(self.working_dir,
                                     self.test_function)
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        fail_script = join_path(self.data_dir,
                                'test_fail_script.sh')
        report_file = join_path(self.working_dir,
                                self.test_function+'_report.txt')
        tmp_file = join_path(self.working_dir,
                             self.test_function+'.txt')
        job_params = '3 2000 ' + tmp_file
        self.delete_file(tmp_file)
        jobman = JobManagerHouseKeeping(jobs_report_file=report_file)
        jobman.submit_job(job_name_1,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          fail_script,
                          job_params,
                          )
        job_name_2 = self.test_function + "_2"
        job_params = '100 1200 ' + tmp_file
        jobman.submit_job(job_name_2,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          job_script,
                          job_params,
                          prereq=[job_name_1],
                          )
        job_name_3 = self.test_function + "_3"
        job_params = '410 420 ' + tmp_file
        jobman.submit_job(job_name_3,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          fail_script,
                          job_params,
                          prereq=[job_name_1, job_name_2],
                          )
        jobman.monitor_jobs()

    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
    def test_prereq(self):
        """
        - check if 'afterok' dependency can be applied
        NOTE -> need to manually check the result
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name_1 = self.test_function + "_1"
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        tmp_file = join_path(self.working_dir,
                             self.test_function+'.txt')
        job_params = '3 20000 ' + tmp_file
        self.delete_file(tmp_file)
        jobman = JobManager()
        slurm_log_prefix = join_path(self.working_dir,
                                     job_name_1)
        jobman.submit_job(job_name_1,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          job_script,
                          job_params,
                          )
        job_name_2 = self.test_function + "_2"
        job_params = '100 12000 ' + tmp_file
        slurm_log_prefix = join_path(self.working_dir,
                                     job_name_2)
        jobman.submit_job(job_name_2,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          job_script,
                          job_params,
                          prereq=[job_name_1],
                          )
        job_name_3 = self.test_function + "_3"
        job_params = '410 420 ' + tmp_file
        slurm_log_prefix = join_path(self.working_dir,
                                     job_name_3)
        jobman.submit_job(job_name_3,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          job_script,
                          job_params,
                          prereq=[job_name_1, job_name_2],
                          )
        jobman.cancel_job(job_name_3)
        jobman.cancel_job(job_name_2)
        jobman.cancel_job(job_name_1)

    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
    def test_nodelist(self):
        """ test sbatch '--nodelist' feature """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        tmp_file = join_path(self.working_dir,
                             self.test_function+'.txt')
        job_params = '3 20000 ' + tmp_file
        self.delete_file(tmp_file)
        jobman = JobManager()
        slurm_log_prefix = join_path(self.working_dir,
                                     job_name)
        jobman.submit_job(job_name,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_prefix,
                          job_script,
                          job_params,
                          nodelist="m[159]",
                          )

    def tearDown(self):
        self.remove_working_dir()
