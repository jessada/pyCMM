import unittest
import filecmp
import datetime
from os.path import exists as path_exists
from os.path import join as join_path
from os.path import dirname
from os.path import isdir
from pycmm import settings
from pycmm.template import SafeTester
from pycmm.utils import mylogger
from pycmm.settings import FAST_PROJECT_CODE
from pycmm.settings import SLOW_PROJECT_CODE
from pycmm.flow.gatkbp import GATKBPPipeline
from pycmm.flow.gatkbp import JOBS_SETUP_DATASET_NAME_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_PROJECT_CODE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_KNOWN_INDELS_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_DBSNP_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_REFFERENCE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_OUTPUT_DIR_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_VARIANTS_CALLING_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_JOBS_REPORT_FILE_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY
from pycmm.flow.gatkbp import JOBS_SETUP_DATASET_USAGE_MAIL_KEY
from pycmm.flow.gatkbp import create_jobs_setup_file

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
                                 sample_group='test_group',
                                 project_code=SLOW_PROJECT_CODE,
                                 variants_calling="YES",
                                 targets_interval_list=None,
                                 dataset_usage_mail="NO",
                                 sample_usage_mail={},
                                 ):
        jobs_setup_file = join_path(self.working_dir,
                                    self.test_function+'_jobs_setup.txt')
        if dataset_name is None:
            dataset_name = self.test_function
        known_indels = []
        known_indels.append(join_path(self.data_dir,
                                      '1000G_phase1.indels.b37.vcf'))
        known_indels.append(join_path(self.data_dir,
                                      'Mills_and_1000G_gold_standard.indels.b37.vcf'))
        dbsnp_file = join_path(self.data_dir,
                               'dbsnp_138.b37.vcf')
        reference_file = join_path(self.data_dir,
                                   'ref.fa')
        create_jobs_setup_file(dataset_name=dataset_name,
                               sample_group=sample_group,
                               project_code=project_code,
                               reference_file=reference_file,
                               project_out_dir=self.working_dir,
                               samples_root_dir=self.data_dir,
                               known_indels_file=known_indels,
                               dbsnp_file=dbsnp_file,
                               targets_interval_list=targets_interval_list,
                               dataset_usage_mail=dataset_usage_mail,
                               sample_usage_mail=sample_usage_mail,
                               out_jobs_setup_file=jobs_setup_file,
                               )
        return jobs_setup_file

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
        self.assertEqual(pl.dbsnp,
                         "/glob/jessada/private/master_data/known_indels_SNPs/dbsnp.vcf",
                         "GATKBPPipeline cannot correctly read meta info 'dbsnp' from jobs setup file")
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
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample3'].library,
                         'lib1',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample3'].unit,
                         'unit1',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample3'].fastq_pairs[0]['R1'],
                         '/proj/b2012247/private/axeq_custom/1303AHS-0016/Co-1116_1.fastq.gz',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample_1'].sample_name,
                         'test_sample_1',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample_1'].preprocess_sample,
                         True,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample_1'].platform,
                         'Illumina',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample_1'].fastq_pairs[3]['R2'],
                         '/proj/b2012247/private/axeq_custom/1303AHS-0016/Co-1747_2.fastq.gz',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample10000'].sample_name,
                         'test_sample10000',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample10000'].sample_group,
                         'test_group',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples['test_sample10000'].usage_mail,
                         False,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_bwa_mem_1(self):
        """ test basic bwa mem, a normal pair of fastq files """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
        sample_name = "test-M_sample"
        pl = GATKBPPipeline(jobs_setup_file)
        pl.bwa_mem(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_bwa_mem_4(self):
        """ test bwa mem with concatenate gz files """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_name = "test-L_sample"
        sample_usage_mail = {sample_name: "YES"}
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
                                                        sample_usage_mail=sample_usage_mail)
        pl = GATKBPPipeline(jobs_setup_file)
        pl.bwa_mem(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_bwa_mem_5(self):
        """ test bwa mem with problematics Ill-Br """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_name = "test-Ill-Prob"
        sample_usage_mail = {sample_name: "YES"}
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
#                                                        sample_usage_mail=sample_usage_mail,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.bwa_mem(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_sort_sam_1(self):
        """ test basic sort sam, a sam file from one pair of fastq files """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        pl = GATKBPPipeline(jobs_setup_file)
        sam_file = join_path(self.data_dir,
                             sample_name+"_concatted_reads.sam")
        pl.sort_sam(sample_name, sam_file=sam_file)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_sort_sam_2(self):
        """ test a little advanced sort sam, 6 pairs of fastq files """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-6-pairs"
        pl = GATKBPPipeline(jobs_setup_file)
        sam_file = join_path(self.data_dir,
                             sample_name+"_concatted_reads.sam")
        pl.sort_sam(sample_name, sam_file=sam_file)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_sort_sam_3(self):
        """ test a little more advanced sort sam, 6 abnormal pairs of fastq files """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-6-unpairs"
        pl = GATKBPPipeline(jobs_setup_file)
        sam_file = join_path(self.data_dir,
                             sample_name+"_concatted_reads.sam")
        pl.sort_sam(sample_name, sam_file=sam_file)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_sort_sam_4(self):
        """ test sort sam with concatenate gz files """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_name = "test-sample"
        pl = GATKBPPipeline(jobs_setup_file)
        sam_file = join_path(self.data_dir,
                             sample_name+"_concatted_reads.sam")
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
    def test_create_intervals_1(self):
        """ test standard creating indels target intervals """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
        sample_name = "test-sample"
        pl = GATKBPPipeline(jobs_setup_file)
        bam_file = join_path(self.data_dir,
                             sample_name+"_dedup_reads.bam")
        pl.create_intervals(sample_name, bam_file=bam_file)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_create_intervals_2(self):
        """ test creating indels target intervals with bam file that has extremely high score"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
        sample_name = "hi_score_sample"
        pl = GATKBPPipeline(jobs_setup_file)
        bam_file = join_path(self.data_dir,
                             sample_name+"_dedup_reads.bam")
        pl.create_intervals(sample_name, bam_file=bam_file)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_indels_realign_1(self):
        """ test standard indel realignment """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
        sample_name = "test-sample"
        pl = GATKBPPipeline(jobs_setup_file)
        bam_file = join_path(self.data_dir,
                             sample_name+"_dedup_reads.bam")
        indels_target_intervals_file = join_path(self.data_dir,
                                                 'indels_target_intervals.list')
        pl.indels_realign(sample_name,
                          bam_file=bam_file,
                          indels_target_intervals_file=indels_target_intervals_file,
                          )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_indels_realign_2(self):
        """ test indel realignment with high quality score"""

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
        sample_name = "hi_score_sample"
        pl = GATKBPPipeline(jobs_setup_file)
        bam_file = join_path(self.data_dir,
                             sample_name+"_dedup_reads.bam")
        indels_target_intervals_file = join_path(self.data_dir,
                                                 'indels_target_intervals.list')
        pl.indels_realign(sample_name,
                          bam_file=bam_file,
                          indels_target_intervals_file=indels_target_intervals_file,
                          )

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_base_recal_1(self):
        """ test standard base recalibration """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
        sample_name = "test-sample"
        pl = GATKBPPipeline(jobs_setup_file)
        bam_file = join_path(self.data_dir,
                             sample_name+"_realigned_reads.bam")
        pl.base_recal(sample_name, bam_file=bam_file)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_base_recal_2(self):
        """ test base recalibration with high quality score """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
        sample_name = "hi_score_sample"
        pl = GATKBPPipeline(jobs_setup_file)
        bam_file = join_path(self.data_dir,
                             sample_name+"_realigned_reads.bam")
        pl.base_recal(sample_name, bam_file=bam_file)

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
                                     sample_name+"_recal_reads.bam")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          usage_mail="YES",
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
                                     sample_name+"_recal_reads.bam")
        self.__add_sample_jobs_setup_file(jobs_setup_file,
                                          sample_name,
                                          "NA",
                                          "NA",
                                          "NA",
                                          usage_mail="YES",
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
                           FAST_PROJECT_CODE,
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
        sample_name = "test-M_sample"
        sample_usage_mail = {sample_name: "YES"}
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
                                                   #     sample_usage_mail=sample_usage_mail,
                                                        targets_interval_list=targets_interval_list,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_sample(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_2(self):
        """ test sample pre-processing workflow (whole exome) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code='b2011097')
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
    def test_preprocess_sample_4(self):
        """ test sample pre-processing workflow (Illumina WES) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_name = "test-MR_sample"
        sample_usage_mail = {sample_name: "YES"}
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(project_code=SLOW_PROJECT_CODE,
                                                        sample_usage_mail=sample_usage_mail,
#                                                        targets_interval_list=targets_interval_list,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_sample(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_5(self):
        """ test sample pre-processing workflow (one of problematic sample from Illimina WES) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_name = "test-Ill_prob"
        sample_usage_mail = {sample_name: "YES"}
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_sample(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_Br_652(self):
        """ test sample pre-processing workflow (one of problematic sample from Illimina WES) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_name = "Br-652"
        sample_usage_mail = {sample_name: "YES"}
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_sample(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_Br_724(self):
        """ test sample pre-processing workflow (one of problematic sample from Illimina WES) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_name = "Br-724"
        sample_usage_mail = {sample_name: "YES"}
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_sample(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_Br_732(self):
        """ test sample pre-processing workflow (one of problematic sample from Illimina WES) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_name = "Br-732"
        sample_usage_mail = {sample_name: "YES"}
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_sample(sample_name)

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_dataset_1(self):
        """ test preprocess target sequencing dataset """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
                                                        dataset_usage_mail="YES",
                                                        )
        pl = GATKBPPipeline(jobs_setup_file)
        pl.preprocess_dataset()

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_dataset_2(self):
        """
        test dataset pre-processing in the case that only one out of three
        need to be processed
        """

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
            if sample_idx > 2:
                self.__add_sample_jobs_setup_file(jobs_setup_file,
                                                  sample_name,
                                                  fastq1_file,
                                                  fastq2_file,
                                                  sample_group,
                                                  )
            else:
                self.__add_sample_jobs_setup_file(jobs_setup_file,
                                                  sample_name,
                                                  fastq1_file,
                                                  fastq2_file,
                                                  sample_group,
                                                  preprocess_sample="NO",
                                                  )
        pl = GATKBPPipeline(jobs_setup_file)
        samples_working_dir = join_path(self.working_dir,
                                        "tmp")
        for sample_idx in xrange(1,3):
            sample_name = "test_sample_M_" + str(sample_idx)
            sample_working_dir = join_path(samples_working_dir,
                                           sample_name)
            raw_vcf_source = join_path(self.data_dir,
                                       sample_name+"_raw.g.vcf")
            raw_vcf_dest = join_path(sample_working_dir,
                                     sample_name+"_raw.g.vcf")
            self.copy_file(raw_vcf_source, raw_vcf_dest)
            raw_idx_source = join_path(self.data_dir,
                                       sample_name+"_raw.g.vcf.idx")
            raw_idx_dest = join_path(sample_working_dir,
                                     sample_name+"_raw.g.vcf.idx")
            self.copy_file(raw_idx_source, raw_idx_dest)
        pl.preprocess_dataset()

    @unittest.skipUnless(settings.SLURM_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_dataset_3(self):
        """ test preprocess WES dataset """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(project_code=SLOW_PROJECT_CODE,
                                                        dataset_usage_mail="YES",
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

class TestFunctions(SafeTester):

    def __init__(self, test_name):
        SafeTester.__init__(self,
                            test_name,
                            dirname(__file__),
                            )

    def setUp(self):
        self.module_name = 'gatkbp'

    @unittest.skipUnless(isdir("/proj/b2011117"), "This can only run in UPPMAX")
    def test_create_jobs_setup_file(self):
        """ test if function create_jobs_setup_file can be really used in production """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        samples_root_dir = self.data_dir
        out_jobs_setup_file = join_path(self.working_dir,
                                        'jobs_setup_file.txt')
        exp_jobs_setup_file = join_path(self.data_dir,
                                        'expected_jobs_setup_file.txt')
        create_jobs_setup_file(dataset_name=self.test_function,
                              sample_group='NA',
                              project_code=SLOW_PROJECT_CODE,
                              reference_file='/links/to/reference_file',
                              project_out_dir=self.working_dir,
                              jobs_report_file='/path/to/job_report_file',
                              samples_root_dir=samples_root_dir,
                              out_jobs_setup_file=out_jobs_setup_file,
                              )
        self.assertTrue(filecmp.cmp(out_jobs_setup_file,
                                    exp_jobs_setup_file),
                        "create_jobs_setup_file doesn't funciton correctly")

    def tearDown(self):
        self.remove_working_dir()
