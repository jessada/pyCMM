import unittest
import filecmp
import datetime
from os.path import join as join_path
from os.path import dirname
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.utils import mylogger
from pycmm.flow.gatkbpdnaseq import bwa_mem as gatk_bwa_mem
from pycmm.flow.gatkbpdnaseq import sort_sam as gatk_sort_sam
from pycmm.flow.gatkbpdnaseq import mark_dup as gatk_mark_dup
from pycmm.flow.gatkbpdnaseq import haplotype_caller as gatk_haplotype_caller
from pycmm.flow.gatkbpdnaseq import combine_gvcfs as gatk_combine_gvcfs
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
        """ test mark duplicate """

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
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        gatk_haplotype_caller(job_name,
                              TEST_PROJECT_CODE,
                              slurm_log_file,
                              dedup_reads,
                              out_gvcf,
                              ref,
                              targets_interval_list=targets_interval_list,
                              email=True,
                              )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_combine_gvcfs_1(self):
        """ test combine small gvcfs """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        gvcf_file_1 = join_path(self.data_dir,
                                "test_M_sample_1.g.vcf")
        gvcf_file_2 = join_path(self.data_dir,
                                "test_M_sample_2.g.vcf")
        gvcf_file_3 = join_path(self.data_dir,
                                "test_M_sample_3.g.vcf")
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        gvcf_list = [gvcf_file_1, gvcf_file_2, gvcf_file_3]
        ref = join_path(self.data_dir,
                        "ref.fa")
        out_merged_gvcf = join_path(self.working_dir,
                                    self.test_function+"_merged.gvcf")
        gatk_combine_gvcfs(job_name,
                           TEST_PROJECT_CODE,
                           slurm_log_file,
                           gvcf_list,
                           out_merged_gvcf,
                           ref,
                           self.working_dir,
                           email=True,
                           )

    @unittest.skip("reserved for load test")
    def test_combine_gvcfs_2(self):
        """ test combine gvcfs """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        gvcf_file_prefix = join_path(self.data_dir,
                                     self.test_function+"_")
        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')

        gvcf_list = map(lambda x: gvcf_file_prefix+"{:03d}".format(x)+".gvcf",
                xrange(1,401))
        ref = join_path(self.data_dir,
                        "ref.fa")
        out_merged_gvcf = join_path(self.working_dir,
                                    self.test_function+"_merged.gvcf")
        gatk_combine_gvcfs(job_name,
                           TEST_PROJECT_CODE,
                           slurm_log_file,
                           gvcf_list,
                           out_merged_gvcf,
                           ref,
                           self.working_dir,
                           )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_genotype_gvcfs(self):
        """ test genotype gvcfs """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        gvcf_file = join_path(self.data_dir,
                              "test_M_sample.g.vcf")

        slurm_log_file = join_path(self.working_dir,
                                   self.test_function+'.log')
        gvcf_list = [gvcf_file]
        ref = join_path(self.data_dir,
                        "ref.fa")
        out_genotyped_vcf = join_path(self.working_dir,
                                    self.test_function+"_genotyped.vcf")
        gatk_combine_gvcfs(job_name,
                           TEST_PROJECT_CODE,
                           slurm_log_file,
                           gvcf_list,
                           out_genotyped_vcf,
                           ref,
                           self.working_dir,
                           email=True,
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
                          self.working_dir,
                          time_stamp,
                          out_gvcf,
                          email=True,
                          )

    def tearDown(self):
        self.remove_working_dir()
