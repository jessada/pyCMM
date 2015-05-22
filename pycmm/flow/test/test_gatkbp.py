import unittest
import filecmp
import datetime
from os.path import join as join_path
from os.path import dirname
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.utils import mylogger
#from pycmm.flow.gatkbp import bwa_mem as gatk_bwa_mem
#from pycmm.flow.gatkbp import sort_sam as gatk_sort_sam
#from pycmm.flow.gatkbp import mark_dup as gatk_mark_dup
#from pycmm.flow.gatkbp import haplotype_caller as gatk_haplotype_caller
#from pycmm.flow.gatkbp import combine_gvcfs as gatk_combine_gvcfs
#from pycmm.flow.gatkbp import preprocess_sample
from pycmm.flow.gatkbp import GATKBPPipeline
from pycmm.flow.gatkbp import JOBS_SETUP_DATASET_NAME_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_PROJECT_CODE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_REFFERENCE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_OUTPUT_DIR_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_VARIANTS_CALLING_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_JOBS_REPORT_FILE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_USAGE_MAIL_KEY

TEST_PROJECT_CODE = 'b2011097'

class TestGATKBPPipeline(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'gatkbp'

    def __create_jobs_setup_file(self,
                                 dataset_name=None,
                                 project_code=TEST_PROJECT_CODE,
                                 variants_calling="YES",
                                 usage_mail="NO",
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        reference_file = join_path(self.data_dir,
                                   'ref.fa')
        jobs_report_file = join_path(self.working_dir,
                                     self.test_function+'_rpt.txt')
        f_rpt = open(jobs_setup_file, "w")
        if dataset_name is None:
            f_rpt.write("##" + JOBS_SETUP_DATASET_NAME_KEY + "=" + self.test_function + "\n")
        else:
            f_rpt.write("##" + JOBS_SETUP_DATASET_NAME_KEY + "=" + dataset_name + "\n")
        f_rpt.write("##" + JOBS_SETUP_PROJECT_CODE_KEY + "=" + project_code + "\n")
        f_rpt.write("##" + JOBS_SETUP_REFFERENCE_KEY + "=" + reference_file + "\n")
        f_rpt.write("##" + JOBS_SETUP_OUTPUT_DIR_KEY + "=" + self.working_dir + "\n")
        f_rpt.write("##" + JOBS_SETUP_VARIANTS_CALLING_KEY + "=" + variants_calling + "\n")
        f_rpt.write("##" + JOBS_SETUP_JOBS_REPORT_FILE_KEY + "=" + jobs_report_file + "\n")
        f_rpt.write("##" + JOBS_SETUP_USAGE_MAIL_KEY + "=" + usage_mail + "\n")
        sample_header = "#SAMPLE"
        sample_header += "\tFASTQ1"
        sample_header += "\tFASTQ2"
        sample_header += "\tSAMPLE_GROUP"
        sample_header += "\tPLATFORM"
        sample_header += "\tLIBRARY"
        sample_header += "\tUNIT"
        sample_header += "\tUSAGE_MAIL"
        sample_header += "\tPREPROCESS_SAMPLE"
        sample_header += "\n"
        f_rpt.write(sample_header)
        f_rpt.close()
        return jobs_setup_file

    def __add_sample_jobs_setup_file(self,
                                     jobs_setup_file,
                                     sample_name,
                                     fastq1,
                                     fastq2,
                                     sample_group,
                                     usage_mail="NO",
                                     preprocess_sample="YES",
                                     ):
        sample_info = []
        sample_info.append(sample_name)
        sample_info.append(fastq1)
        sample_info.append(fastq2)
        sample_info.append(sample_group)
        sample_info.append("Illumina")
        sample_info.append("lib1")
        sample_info.append("unit1")
        sample_info.append(usage_mail)
        sample_info.append(preprocess_sample)
        f_rpt = open(jobs_setup_file, "a")
        f_rpt.write("\t".join(sample_info) + "\n")
        f_rpt.close()

    def test_load_jobs_info(self):
        """ test if meta data and sample info are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        jobs_setup_file = join_path(self.data_dir,
                                    self.test_function+'.txt')
        pl = GATKBPPipeline(jobs_setup_file)
        self.assertEqual(pl.dataset_name,
                         "test-dataset",
                         "GATKBPPipeline cannot correctly read meta info 'dataset name' from jobs setup file")
        self.assertEqual(pl.project_code,
                         "b1234567",
                         "GATKBPPipeline cannot correctly read meta info 'project code' from jobs setup file")
        self.assertEqual(pl.reference,
                         "/glob/jessada/Homo_sapiens.GRCh37.57.dna.concat.fa",
                         "GATKBPPipeline cannot correctly read meta info 'reference' from jobs setup file")
        self.assertEqual(pl.output_dir,
                         "/glob/jessada/private/projects/",
                         "GATKBPPipeline cannot correctly read meta info 'output directory' from jobs setup file")
        self.assertEqual(pl.variants_calling,
                         True,
                         "GATKBPPipeline cannot correctly read meta info 'variants calling' from jobs setup file")
        self.assertEqual(pl.jobs_report_file,
                         "/glob/jessada/job_rpt.txt",
                         "GATKBPPipeline cannot correctly read meta info 'jobs report file' from jobs setup file")
        self.assertEqual(pl.usage_mail,
                         False,
                         "GATKBPPipeline cannot correctly read meta info 'usage mail' from jobs setup file")
        self.assertEqual(len(pl.samples),
                         3,
                         "GATKBPPipeline cannot correctly identify number of sampels from jobs setup file")
        self.assertEqual(pl.samples['test_sample3'].sample_name,
                         'test_sample3',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample3'].library,
                         'lib1',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample3'].unit,
                         'unit1',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample3'].fastq1_file,
                         '/proj/b2012247/private/axeq_custom/1303AHS-0016/Co-1116_1.fastq.gz',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample_1'].sample_name,
                         'test_sample_1',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample_1'].platform,
                         'Illumina',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample_1'].fastq2_file,
                         '/proj/b2012247/private/axeq_custom/1303AHS-0016/Co-1747_2.fastq.gz',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample10000'].sample_name,
                         'test_sample10000',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample10000'].sample_group,
                         'test_group',
                         "GATKBPPipeline cannot correctly sample information from jobs setup file")

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_bwa_mem(self):
        """ test bwa mem """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        sample_group = self.current_func_name
        fastq1_file = join_path(self.data_dir,
                           sample_name + "_1.fastq.gz")
        fastq2_file = join_path(self.data_dir,
                           sample_name + "_2.fastq.gz")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          fastq1_file,
                                          fastq2_file,
                                          sample_group,
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.bwa_mem(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_sort_sam(self):
        """ test sort sam """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        sam_file = join_path(self.data_dir,
                             sample_name+"_raw_aligned_reads.sam")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.sort_sam(sample_name, sam_file=sam_file)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_mark_dup(self):
        """ test mark duplicate """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        bam_file = join_path(self.data_dir,
                             sample_name+"_sorted_reads.bam")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.mark_dup(sample_name, bam_file=bam_file)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_haplotype_caller(self):
        """ test haplotype caller """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        dedup_reads_file = join_path(self.data_dir,
                                     sample_name+"_dedup_reads.bam")
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.hap_cal(sample_name, 
                   bam_file=dedup_reads_file,
                   targets_interval_list=targets_interval_list)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_combine_gvcfs_1(self):
        """ test combine small gvcfs """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        gvcf_file_1 = join_path(self.data_dir,
                                "test_M_sample_1.g.vcf")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          "test_M_sample_1",
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        gvcf_file_2 = join_path(self.data_dir,
                                "test_M_sample_2.g.vcf")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          "test_M_sample_2",
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        gvcf_file_3 = join_path(self.data_dir,
                                "test_M_sample_3.g.vcf")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          "test_M_sample_3",
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        gvcf_list = [gvcf_file_1, gvcf_file_2, gvcf_file_3]
        pl = GATKBPPipeline(jobs_setup_file)
        pl.combine_gvcfs(gvcf_list=gvcf_list)

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
        jobs_setup_file = self.__create_jobs_setup_file()
        gvcf_file = join_path(self.data_dir,
                              "test_M_sample.g.vcf")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          "test_M_sample",
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        gvcf_list = [gvcf_file]
        pl = GATKBPPipeline(jobs_setup_file)
        pl.genotype_gvcfs(gvcf_list=gvcf_list)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_1(self):
        """ test sample pre-processing workflow (targeted sequencing) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test_small_sample"
        sample_group = self.current_func_name
        fastq1_file = join_path(self.data_dir,
                           sample_name + "_1.fastq.gz")
        fastq2_file = join_path(self.data_dir,
                           sample_name + "_2.fastq.gz")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          fastq1_file,
                                          fastq2_file,
                                          sample_group,
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_sample(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_2(self):
        """ test sample pre-processing workflow (whole exome) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test_big_sample"
        sample_group = self.current_func_name
        fastq1_file = join_path(self.data_dir,
                           sample_name + "_1.fastq.gz")
        fastq2_file = join_path(self.data_dir,
                           sample_name + "_2.fastq.gz")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          fastq1_file,
                                          fastq2_file,
                                          sample_group,
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_sample(sample_name)


    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_dataset_1(self):
        """ test dataset pre-processing workflow (targeted sequencing) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_group = self.current_func_name
        for sample_idx in xrange(1,4):
            sample_name = "test_sample_M_" + str(sample_idx)
            fastq1_file = join_path(self.data_dir,
                               sample_name + "_1.fastq.gz")
            fastq2_file = join_path(self.data_dir,
                               sample_name + "_2.fastq.gz")
            self.__add_sample_jobs_setup_file(jobs_setup_file,
                                              sample_name,
                                              fastq1_file,
                                              fastq2_file,
                                              sample_group,
                                              )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_dataset()

