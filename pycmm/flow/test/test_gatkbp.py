import unittest
import filecmp
import datetime
from os.path import exists as path_exists
from os.path import join as join_path
from os.path import dirname
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.utils import mylogger
from pycmm.flow.gatkbp import GATKBPPipeline
from pycmm.flow.gatkbp import JOBS_SETUP_DATASET_NAME_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_PROJECT_CODE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_KNOWN_INDELS_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_REFFERENCE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_OUTPUT_DIR_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_VARIANTS_CALLING_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_JOBS_REPORT_FILE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_USAGE_MAIL_KEY
#from pycmm.settings import KNOWN_INDELS_1000G_PHASE1
#from pycmm.settings import KNOWN_INDELS_MILLS_AND_1000G_GOLD_STANDARD

TEST_PROJECT_CODE = 'b2011158'

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
#                                 known_indels=[KNOWN_INDELS_1000G_PHASE1,
#                                               KNOWN_INDELS_MILLS_AND_1000G_GOLD_STANDARD,
#                                               ],
                                 variants_calling="YES",
                                 targets_interval_list=None,
                                 usage_mail="NO",
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        known_indels = join_path(self.data_dir,
                                 '1000g_phase1.vcf')
        known_indels += "," + join_path(self.data_dir,
                                        'mills_n_1000g_gold.vcf')
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
        f_rpt.write("##" + JOBS_SETUP_KNOWN_INDELS_KEY + "=" + known_indels + "\n")
        f_rpt.write("##" + JOBS_SETUP_REFFERENCE_KEY + "=" + reference_file + "\n")
        f_rpt.write("##" + JOBS_SETUP_OUTPUT_DIR_KEY + "=" + self.working_dir + "\n")
        f_rpt.write("##" + JOBS_SETUP_VARIANTS_CALLING_KEY + "=" + variants_calling + "\n")
        f_rpt.write("##" + JOBS_SETUP_JOBS_REPORT_FILE_KEY + "=" + jobs_report_file + "\n")
        if targets_interval_list is not None:
            f_rpt.write("##" + JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY + "=" + targets_interval_list + "\n")
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
        self.assertEqual(pl.known_indels[0],
                         "/glob/jessada/private/master_data/known_indels_SNPs/1000G_phase1.indels.b37.vcf",
                         "GATKBPPipeline cannot correctly read meta info 'known indels[0]' from jobs setup file")
        self.assertEqual(pl.known_indels[1],
                         "/glob/jessada/private/master_data/known_indels_SNPs/Mills_and_1000G_gold_standard.indels.b37.vcf",
                         "GATKBPPipeline cannot correctly read meta info 'known indels[1]' from jobs setup file")
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
        self.assertEqual(pl.targets_interval_list,
                         "/glob/jessada/target_interval",
                         "GATKBPPipeline cannot correctly read meta info 'targets.interval_list' from jobs setup file")
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
    def test_indels_target_intervals(self):
        """ test creating indels target intervals """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        bam_file = join_path(self.data_dir,
                             sample_name+"_dedup_reads.bam")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          usage_mail="YES",
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.indels_target_intervals(sample_name, bam_file=bam_file)

#    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_indels_realign(self):
        """ test creating indels target intervals """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        bam_file = join_path(self.data_dir,
                             sample_name+"_dedup_reads.bam")
        indels_target_intervals = join_path(self.data_dir,
                                            'indels_target_intervals.list')
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          usage_mail="YES",
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.indels_realign(sample_name,
                          bam_file=bam_file,
                          indels_target_intervals=indels_target_intervals,
                          )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_haplotype_caller_1(self):
        """ test haplotype caller with small data size """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(targets_interval_list=targets_interval_list)
        sample_name = "test-small-sample"
        dedup_reads_file = join_path(self.data_dir,
                                     sample_name+"_dedup_reads.bam")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.hap_cal(sample_name, 
                   bam_file=dedup_reads_file,
                   )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_haplotype_caller_2(self):
        """ test haplotype caller with big data size """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-big-sample"
        dedup_reads_file = join_path(self.data_dir,
                                     sample_name+"_dedup_reads.bam")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.hap_cal(sample_name, 
                   bam_file=dedup_reads_file,
                   )

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
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(targets_interval_list=targets_interval_list)
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
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(targets_interval_list=targets_interval_list)
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

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_garbage_collecting(self):
        """ test if the workflow can correctly collecting garbage """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          )
        pl = GATKBPPipeline(jobs_setup_file)
        sample_rec = pl.samples[sample_name]
        raw_aligned_source = join_path(self.data_dir,
                                       sample_name+"_raw_aligned_reads.sam")
        raw_aligned_dest = join_path(self.working_dir,
                                     sample_rec.raw_aligned_reads_file)
        self.copy_file(raw_aligned_source, raw_aligned_dest)
        sorted_reads_source = join_path(self.data_dir,
                                        sample_name+"_sorted_reads.bam")
        sorted_reads_dest = join_path(self.working_dir,
                                      sample_rec.sorted_reads_file)
        self.copy_file(sorted_reads_source, sorted_reads_dest)
        pl.mark_dup(sample_name)
        pl.monitor_jobs()
        self.assertFalse(path_exists(raw_aligned_dest),
                         " garbage collecting process doesn't work properly")
        self.assertFalse(path_exists(sorted_reads_dest),
                         " garbage collecting process doesn't work properly")

    def tearDown(self):
        self.remove_working_dir()
