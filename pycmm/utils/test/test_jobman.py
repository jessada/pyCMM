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

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            test_module_name=__name__,
                            )

    def setUp(self):
        pass

    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
    def test_submit_job_1(self):
        """ check if a job has been correctly submitted to SLURM """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
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
                          slurm_log_file,
                          job_script,
                          job_params,
                          )
        job_id = jobman.get_job_id(job_name)
        self.assertTrue(job_id.isdigit(),
                        "Invalid number of queried records")
        jobman.cancel_job(job_name)

#    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
#    def test_submit_job_2(self):
#        """ check if a real GATK alignment job can be submitted """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#
#        gatk_alloc_time = "72:00:00"
#        sample_name = "test_sample"
#        sample_prefix = join_path(self.data_dir,
#                                  sample_name)
#        sample_group = self.current_func_name
#        ref = join_path(self.data_dir,
#                        "ref.fa")
#        raw_aligned_reads = join_path(self.working_dir,
#                                      sample_name+"_raw_aligned.sam")
#        sorted_reads = join_path(self.working_dir,
#                                 sample_name+"_sorted.bam")
#        dedup_reads = join_path(self.working_dir,
#                                sample_name+"_dedup.bam")
#        raw_gvcf = join_path(self.working_dir,
#                             sample_name+"_raw.gvcf")
#        metrics_file = join_path(self.working_dir,
#                                 sample_name+"_metrics.txt")
#        jobman = JobManager()
#
#        bwa_mem_job_name = "test_bwa_mem"
#        slurm_log_file = join_path(self.working_dir,
#                                   bwa_mem_job_name+'.log')
#        job_script = "$PYCMM/bash/GATKBP_bwa_mem.sh"
#        job_params = " -S " + sample_prefix 
#        job_params += " -g " + sample_group
#        job_params += " -r " + ref
#        job_params += " -o " + raw_aligned_reads
#        jobman.submit_job(bwa_mem_job_name,
#                          TEST_PROJECT_CODE,
#                          TEST_PARTITION_TYPE,
#                          "1",
#                          gatk_alloc_time,
#                          slurm_log_file,
#                          job_script,
#                          job_params,
#                          )
#
#        sort_reads_job_name = "test_sorted_reads"
#        slurm_log_file = join_path(self.working_dir,
#                                   sort_reads_job_name+'.log')
#        job_script = "$PYCMM/bash/GATKBP_sort_reads.sh"
#        job_params = " -S " + raw_aligned_reads 
#        job_params += " -o " + sorted_reads
#        jobman.submit_job(sort_reads_job_name,
#                          TEST_PROJECT_CODE,
#                          TEST_PARTITION_TYPE,
#                          "2",
#                          gatk_alloc_time,
#                          slurm_log_file,
#                          job_script,
#                          job_params,
#                          prereq=[bwa_mem_job_name],
#                          )
#
#        dedup_reads_job_name = "test_dedup_reads"
#        slurm_log_file = join_path(self.working_dir,
#                                   dedup_reads_job_name+'.log')
#        job_script = "$PYCMM/bash/GATKBP_markdup.sh"
#        job_params = " -B " + sorted_reads 
#        job_params += " -o " + dedup_reads
#        job_params += " -m " + metrics_file
#        jobman.submit_job(dedup_reads_job_name,
#                          TEST_PROJECT_CODE,
#                          TEST_PARTITION_TYPE,
#                          "4",
#                          gatk_alloc_time,
#                          slurm_log_file,
#                          job_script,
#                          job_params,
#                          prereq=[sort_reads_job_name],
#                          )
#
#        hap_cal_job_name = "test_hap_cal"
#        slurm_log_file = join_path(self.working_dir,
#                                   hap_cal_job_name+'.log')
#        job_script = "$PYCMM/bash/GATKBP_haplotype_caller.sh"
#        job_params = " -B " + dedup_reads 
#        job_params += " -o " + raw_gvcf
#        job_params += " -r " + ref
#        jobman.submit_job(hap_cal_job_name,
#                          TEST_PROJECT_CODE,
#                          TEST_PARTITION_TYPE,
#                          "2",
#                          gatk_alloc_time,
#                          slurm_log_file,
#                          job_script,
#                          job_params,
#                          prereq=[dedup_reads_job_name],
#                          )
#
#        job_id = jobman.get_job_id(bwa_mem_job_name)
#        jobman.cancel_job(bwa_mem_job_name)
#
    @unittest.skipUnless(SLURM_TEST, "taking too long time to test")
    def test_get_job_status(self):
        """ check if job status can be correctly retrieved """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
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
                          slurm_log_file,
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
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
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
                          slurm_log_file,
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
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
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
                          slurm_log_file,
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
                          slurm_log_file,
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
                          slurm_log_file,
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
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        job_script = join_path(self.data_dir,
                               'test_job_script.sh')
        tmp_file = join_path(self.working_dir,
                             self.test_function+'.txt')
        job_params = '3 20000 ' + tmp_file
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
        job_params = '100 12000 ' + tmp_file
        jobman.submit_job(job_name_2,
                          TEST_PROJECT_CODE,
                          TEST_PARTITION_TYPE,
                          TEST_NTASKS,
                          TEST_ALLOC_TIME,
                          slurm_log_file,
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
                          slurm_log_file,
                          job_script,
                          job_params,
                          prereq=[job_name_1, job_name_2],
                          )
        jobman.cancel_job(job_name_3)
        jobman.cancel_job(job_name_2)
        jobman.cancel_job(job_name_1)

    def tearDown(self):
        self.remove_working_dir()
