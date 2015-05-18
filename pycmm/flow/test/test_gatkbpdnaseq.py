import unittest
import filecmp
import datetime
from os.path import join as join_path
from os.path import dirname
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.utils import exec_sh
from pycmm.flow.gatkbpdnaseq import bwa_mem as gatk_bwa_mem
from pycmm.flow.gatkbpdnaseq import sort_sam as gatk_sort_sam
from pycmm.flow.gatkbpdnaseq import mark_dup as gatk_mark_dup
from pycmm.flow.gatkbpdnaseq import haplotype_caller as gatk_haplotype_caller
from pycmm.flow.gatkbpdnaseq import preprocess_sample

TEST_PROJECT_CODE = 'b2011097'

class TestFunctions(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'gatkbpdnaseq'

    def __create_db_instance(self):
        return None

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_bwa_mem(self):
        """ test bwa mem """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        sample_name = "test_sample"
        sample_prefix = join_path(self.data_dir,
                                  sample_name)
        sample_group = self.current_func_name
        ref = join_path(self.data_dir,
                        "ref.fa")
        raw_aligned_reads = join_path(self.working_dir,
                                      sample_name+"_raw_aligned.sam")
        gatk_bwa_mem(job_name,
                     TEST_PROJECT_CODE,
                     slurm_log_file,
                     sample_prefix,
                     sample_group,
                     ref,
                     raw_aligned_reads,
                     )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_sort_sam(self):
        """ test sort reads """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        sample_name = "test_sample"
        raw_aligned_reads = join_path(self.data_dir,
                                      sample_name+"_raw_aligned_reads.sam")
        sorted_reads = join_path(self.working_dir,
                                 sample_name+"_sorted_reads.bam")
        gatk_sort_sam(job_name,
                        TEST_PROJECT_CODE,
                        slurm_log_file,
                        raw_aligned_reads,
                        sorted_reads,
                        )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_mark_dup(self):
        """ test dedup reads """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        sample_name = "test_sample"
        sorted_reads = join_path(self.data_dir,
                                 sample_name+"_sorted_reads.bam")
        dedup_reads = join_path(self.working_dir,
                                sample_name+"_dedup_reads.bam")
        metrics_file = join_path(self.working_dir,
                                 sample_name+"_metrics.txt")
        gatk_mark_dup(job_name,
                      TEST_PROJECT_CODE,
                      slurm_log_file,
                      sorted_reads,
                      dedup_reads,
                      metrics_file,
                      )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_haplotype_caller(self):
        """ test haplotype caller """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        sample_name = "test_sample"
        dedup_reads = join_path(self.data_dir,
                                sample_name+"_dedup_reads.bam")
        ref = join_path(self.data_dir,
                        "ref.fa")
        out_gvcf = join_path(self.working_dir,
                             sample_name+".gvcf")
        gatk_haplotype_caller(job_name,
                              TEST_PROJECT_CODE,
                              slurm_log_file,
                              dedup_reads,
                              out_gvcf,
                              ref,
                              )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_1(self):
        """ test sample pre-processing workflow """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        time_stamp = datetime.datetime.now()
        sample_name = "test_sample"
        sample_prefix = join_path(self.data_dir,
                                  sample_name)
        sample_group = self.current_func_name
        ref = join_path(self.data_dir,
                        "ref.fa")
        out_gvcf = join_path(self.working_dir,
                             sample_name+".gvcf")
        preprocess_sample(TEST_PROJECT_CODE,
                          sample_name,
                          sample_prefix,
                          sample_group,
                          ref,
                          self.working_dir,
                          self.working_dir,
                          time_stamp,
                          out_gvcf,
                          )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_2(self):
        """ test sample pre-processing workflow """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        time_stamp = datetime.datetime.now()
        sample_name = "test_big_sample"
        sample_prefix = join_path(self.data_dir,
                                  sample_name)
        sample_group = self.current_func_name
        ref = join_path(self.data_dir,
                        "ref.fa")
        out_gvcf = join_path(self.working_dir,
                             sample_name+".gvcf")
        preprocess_sample(TEST_PROJECT_CODE,
                          sample_name,
                          sample_prefix,
                          sample_group,
                          ref,
                          self.working_dir,
                          self.working_dir,
                          time_stamp,
                          out_gvcf,
                          email=True,
                          )

    def tearDown(self):
        self.remove_working_dir()
