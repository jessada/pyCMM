import unittest
import filecmp
from os.path import exists as path_exists
from os.path import join as join_path
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.template import SafeTester
from pycmm.settings import SLOW_PROJECT_CODE
from pycmm.flow.gatkbp import GATKBPPipeline
from pycmm.flow.gatkbp import create_jobs_setup_file
from pycmm.utils import file_size

DFLT_KNOWN_INDELS = []
DFLT_KNOWN_INDELS.append('1000G_phase1.indels.b37.vcf')
DFLT_KNOWN_INDELS.append('Mills_and_1000G_gold_standard.indels.b37.vcf')
DFLT_DBSNP_FILE = 'dbsnp_138.b37.vcf'
DFLT_REF_FILE = 'ref.fa'

GATKBP_TEST = False

def create_gatk_jobs_setup_file(test_function,
                                working_dir,
                                data_dir,
                                project_code=None,
                                flow_alloc_time=None,
                                indel_recal=True,
                                variants_calling=False,
                                known_indels=DFLT_KNOWN_INDELS,
                                targets_interval_list=None,
                                dataset_usage_mail=False,
                                sample_usage_mail={},
                                ):
    jobs_setup_file = join_path(working_dir, test_function+'_jobs_setup.txt')
    known_indels_file = []
    for item in known_indels:
        known_indels_file.append(join_path(data_dir, item))
    dbsnp_file = join_path(data_dir, DFLT_DBSNP_FILE)
    reference_file = join_path(data_dir, DFLT_REF_FILE)
    create_jobs_setup_file(project_name=test_function,
                           project_code=project_code,
                           flow_alloc_time=flow_alloc_time,
                           reference_file=reference_file,
                           project_out_dir=working_dir,
                           samples_root_dir=data_dir,
                           indel_recal=indel_recal,
                           variants_calling=variants_calling,
                           known_indels_file=known_indels_file,
                           dbsnp_file=dbsnp_file,
                           targets_interval_list=targets_interval_list,
                           dataset_usage_mail=dataset_usage_mail,
                           sample_usage_mail=sample_usage_mail,
                           out_jobs_setup_file=jobs_setup_file,
                           )
    return jobs_setup_file

class TestGATKBPPipeline(SafeTester):

    def __init__(self, methodName):
        super(TestGATKBPPipeline, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self,
                                 project_code=None,
                                 flow_alloc_time=None,
                                 indel_recal=True,
                                 variants_calling=False,
                                 known_indels=DFLT_KNOWN_INDELS,
                                 targets_interval_list=None,
                                 dataset_usage_mail=False,
                                 sample_usage_mail={},
                                 ):
        return create_gatk_jobs_setup_file(test_function=self.test_function,
                                           working_dir=self.working_dir,
                                           data_dir=self.data_dir,
                                           project_code=project_code,
                                           flow_alloc_time=flow_alloc_time,
                                           indel_recal=indel_recal,
                                           variants_calling=variants_calling,
                                           known_indels=known_indels,
                                           targets_interval_list=targets_interval_list,
                                           dataset_usage_mail=dataset_usage_mail,
                                           sample_usage_mail=sample_usage_mail,
                                           )

    def test_load_jobs_info_1(self):
        """ test if default job configurations are loaded correctly """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        job_name = self.test_function
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        exp_project_name = self.test_function
        self.assertEqual(pl.project_name,
                         exp_project_name,
                         "GATKBPPipeline cannot correctly identify 'dataset name' from jobs setup file")
        self.assertEqual(pl.project_code,
                         None,
                         "GATKBPPipeline cannot correctly identify 'project code' from jobs setup file")
        self.assertEqual(pl.flow_alloc_time,
                         None,
                         "GATKBPPipeline cannot correctly identify 'gatk allocation time' from jobs setup file")
        self.assertEqual(pl.gatk_params.known_indels[0],
                         join_path(self.data_dir, DFLT_KNOWN_INDELS[0]),
                         "GATKBPPipeline cannot correctly identify 'known indels[0]' from jobs setup file")
        self.assertEqual(pl.gatk_params.known_indels[1],
                         join_path(self.data_dir, DFLT_KNOWN_INDELS[1]),
                         "GATKBPPipeline cannot correctly identify 'known indels[1]' from jobs setup file")
        self.assertEqual(pl.gatk_params.dbsnp,
                         join_path(self.data_dir, DFLT_DBSNP_FILE),
                         "GATKBPPipeline cannot correctly identify 'dbsnp' from jobs setup file")
        self.assertEqual(pl.gatk_params.reference,
                         join_path(self.data_dir, DFLT_REF_FILE),
                         "GATKBPPipeline cannot correctly identify 'reference' from jobs setup file")
        self.assertEqual(pl.project_out_dir,
                         self.working_dir,
                         "GATKBPPipeline cannot correctly identify 'output directory' from jobs setup file")
        self.assertEqual(pl.gatk_params.indel_recal,
                         True,
                         "GATKBPPipeline cannot correctly identify 'indel recalibration' from jobs setup file")
        self.assertEqual(pl.gatk_params.variants_calling,
                         False,
                         "GATKBPPipeline cannot correctly identify 'variants calling' from jobs setup file")
        self.assertEqual(pl.jobs_report_file,
                         join_path(pl.project_out_dir, exp_project_name+"_rpt.txt"),
                         "GATKBPPipeline cannot correctly identify 'jobs report file' from jobs setup file")
        self.assertEqual(pl.gatk_params.targets_interval_list,
                         None,
                         "GATKBPPipeline cannot correctly identify 'targets.interval_list' from jobs setup file")
        self.assertEqual(pl.gatk_params.dataset_usage_mail,
                         False,
                         "GATKBPPipeline cannot correctly identify 'usage mail' from jobs setup file")
        self.assertEqual(len(pl.samples),
                         3,
                         "GATKBPPipeline cannot correctly identify number of samples from jobs setup file")
        test_sample_id_3 = 'test_sample_3'
        self.assertEqual(pl.samples[test_sample_id_3].sample_id,
                         test_sample_id_3,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples[test_sample_id_3].library,
                         'lib1',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples[test_sample_id_3].unit,
                         'unit1',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        exp_R1_file = join_path(join_path(join_path(self.data_dir,
                                                    'test_group_x'),
                                          test_sample_id_3),
                                test_sample_id_3 + '_R1.fastq.gz')
        self.assertEqual(pl.samples[test_sample_id_3].fastq_pairs[0]['R1'],
                         exp_R1_file,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        test_sample_id_1 = 'test_sample_1'
        self.assertEqual(pl.samples[test_sample_id_1].sample_id,
                         test_sample_id_1,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples[test_sample_id_1].preprocess_sample,
                         True,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples[test_sample_id_1].platform,
                         'Illumina',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples[test_sample_id_1].sample_group,
                         'test_group_1',
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples[test_sample_id_1].usage_mail,
                         False,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        exp_R2_file = join_path(join_path(join_path(self.data_dir,
                                                    'test_group_1'),
                                          test_sample_id_1),
                                test_sample_id_1 + '_R2.fastq.gz')
        self.assertEqual(pl.samples[test_sample_id_1].fastq_pairs[0]['R2'],
                         exp_R2_file,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if modified job configurations are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(project_code=SLOW_PROJECT_CODE,
                                                        indel_recal=False,
                                                        variants_calling=True,
                                                        flow_alloc_time="20:00:00",
                                                        targets_interval_list=targets_interval_list,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        exp_project_name = self.test_function
        self.assertEqual(pl.flow_alloc_time,
                         "20:00:00",
                         "GATKBPPipeline cannot correctly identify 'gatk allocation time' from jobs setup file")
        self.assertEqual(pl.gatk_params.indel_recal,
                         False,
                         "GATKBPPipeline cannot correctly identify 'indel recalibration' from jobs setup file")
        self.assertEqual(pl.gatk_params.variants_calling,
                         True,
                         "GATKBPPipeline cannot correctly identify 'variants calling' from jobs setup file")
        self.assertEqual(pl.gatk_params.targets_interval_list,
                         targets_interval_list,
                         "GATKBPPipeline cannot correctly identify 'targets.interval_list' from jobs setup file")

# ************************************** test very small sequence w/o indel **************************************
    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_bwa_mem_1(self):
        """ test basic bwa mem (w/o indel), a normal pair of small fastq files """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_id = "test-S_sample"
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        result = pl.bwa_mem(sample_id)
        exp_result = join_path(self.data_dir,
                               "exp_result")
        self.assertTrue(filecmp.cmp(result,
                                    exp_result),
                        "bwa-mem doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_sort_sam_1(self):
        """ test basic sort sam (w/o indel), a sam file from one pair of small fastq files """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_id = "test-S_sample"
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        sam_file = join_path(self.data_dir,
                             sample_id+".sam")
        result = pl.sort_sam(sample_id, sam_file=sam_file)
        exp_result = join_path(self.data_dir,
                               "exp_result")
        self.assertTrue(filecmp.cmp(result,
                                    exp_result),
                        "SortSam doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_mark_dup_1(self):
        """ test mark duplicate (w/o indel) (small bam) """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_id = "test-S_sample"
        bam_file = join_path(self.data_dir,
                             sample_id+"_sorted_reads.bam")
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        result = pl.mark_dup(sample_id, bam_file=bam_file)
        exp_result = join_path(self.data_dir,
                               "exp_result")
        self.assertTrue(filecmp.cmp(result,
                                    exp_result),
                        "MarkDuplicates doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_base_recal_1(self):
        """ test standard base recalibration (w/o indel) (small) """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_id = "test-S_sample"
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        bam_file = join_path(self.data_dir,
                             sample_id+"_realigned_reads.bam")
        result = pl.base_recal(sample_id, bam_file=bam_file)
        exp_result = join_path(self.data_dir,
                               "exp_result")
        self.assertTrue(filecmp.cmp(result,
                                    exp_result),
                        "BaseRecalibrator doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_haplotype_caller_1(self):
        """ test haplotype caller with small data size (w/o indel) """

        self.init_test(self.current_func_name)
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(targets_interval_list=targets_interval_list)
        sample_id = "test-S_sample"
        dedup_reads_file = join_path(self.data_dir,
                                     sample_id+"_recal_reads.bam")
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        result = pl.hap_cal(sample_id, 
                            bam_file=dedup_reads_file,
                            )
        exp_result = join_path(self.data_dir,
                               "exp_result")
        self.assertTrue(file_size(result) > 5000,
                        "HaplotypeCaller doesn't function correctly")
# ************************************** test very small sequence w/o indel **************************************

# ************************************** test very small sequence with indel **************************************
    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_bwa_mem_2(self):
        """ test basic bwa mem with indel, a normal pair of small fastq files """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_id = "test-S_sample_w_indel"
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        result = pl.bwa_mem(sample_id)
        exp_result = join_path(self.data_dir,
                               "exp_result")
        self.assertTrue(filecmp.cmp(result,
                                    exp_result),
                        "bwa-mem doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_sort_sam_2(self):
        """ test basic sort sam with indel, a sam file from one pair of small fastq files """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_id = "test-S_sample_w_indel"
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        sam_file = join_path(self.data_dir,
                             sample_id+".sam")
        result = pl.sort_sam(sample_id, sam_file=sam_file)
        exp_result = join_path(self.data_dir,
                               "exp_result")
        self.assertTrue(filecmp.cmp(result,
                                    exp_result),
                        "SortSam doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_mark_dup_2(self):
        """ test mark duplicate with indel (small bam) """

        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        sample_id = "test-S_sample_w_indel"
        bam_file = join_path(self.data_dir,
                             sample_id+"_sorted_reads.bam")
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        result = pl.mark_dup(sample_id, bam_file=bam_file)
        exp_result = join_path(self.data_dir,
                               "exp_result")
        self.assertTrue(filecmp.cmp(result,
                                    exp_result),
                        "MarkDuplicates doesn't function correctly")

#    # too elaborate to setup a test
#    # with known indels, it results in millions of re-alignment sites
#    # without know indels, it results in no re-alignment sites
#    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
#    def test_create_intervals_2(self):
#        """ test standard creating indels target intervals (small) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(known_indels=[])
#        sample_id = "test-S_sample_w_indel"
#        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
#        bam_file = join_path(self.data_dir,
#                             sample_id+"_dedup_reads.bam")
#        result = pl.create_intervals(sample_id, bam_file=bam_file)
#        exp_result = join_path(self.data_dir,
#                               "exp_result")
#        self.assertTrue(filecmp.cmp(result,
#                                    exp_result),
#                        "RealignerTargetCreator doesn't function correctly")
#
#    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
#    def test_indels_realign_2(self):
#        """ test standard indel realignment (small) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(known_indels=[])
#        sample_id = "test-S_sample"
#        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
#        bam_file = join_path(self.data_dir,
#                             sample_id+"_dedup_reads.bam")
#        indels_target_intervals_file = join_path(self.data_dir,
#                                                 'indels_target_intervals.list')
#        pl.indels_realign(sample_id,
#                          bam_file=bam_file,
#                          indels_target_intervals_file=indels_target_intervals_file,
#                          )
# ************************************** test very small sequence with indel **************************************



#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_bwa_mem_4(self):
#        """ test bwa mem with concatenate gz files """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        sample_id = "test-L_sample"
#        sample_usage_mail = {sample_id: "YES"}
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
#                                                        sample_usage_mail=sample_usage_mail)
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.bwa_mem(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_bwa_mem_5(self):
#        """ test bwa mem with problematics Ill-Br """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        sample_id = "test-Ill-Prob"
#        sample_usage_mail = {sample_id: "YES"}
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
##                                                        sample_usage_mail=sample_usage_mail,
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.bwa_mem(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_sort_sam_2(self):
#        """ test a little advanced sort sam, 6 pairs of fastq files """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        sample_id = "test-6-pairs"
#        pl = GATKBPPipeline(jobs_setup_file)
#        sam_file = join_path(self.data_dir,
#                             sample_id+"_concatted_reads.sam")
#        pl.sort_sam(sample_id, sam_file=sam_file)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_sort_sam_3(self):
#        """ test a little more advanced sort sam, 6 abnormal pairs of fastq files """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        sample_id = "test-6-unpairs"
#        pl = GATKBPPipeline(jobs_setup_file)
#        sam_file = join_path(self.data_dir,
#                             sample_id+"_concatted_reads.sam")
#        pl.sort_sam(sample_id, sam_file=sam_file)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_sort_sam_4(self):
#        """ test sort sam with concatenate gz files """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        sample_id = "test-sample"
#        pl = GATKBPPipeline(jobs_setup_file)
#        sam_file = join_path(self.data_dir,
#                             sample_id+"_concatted_reads.sam")
#        pl.sort_sam(sample_id, sam_file=sam_file)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_create_intervals_2(self):
#        """ test creating indels target intervals with bam file that has extremely high score"""
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
#        sample_id = "hi_score_sample"
#        pl = GATKBPPipeline(jobs_setup_file)
#        bam_file = join_path(self.data_dir,
#                             sample_id+"_dedup_reads.bam")
#        pl.create_intervals(sample_id, bam_file=bam_file)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_indels_realign_2(self):
#        """ test indel realignment with high quality score"""
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
#        sample_id = "hi_score_sample"
#        pl = GATKBPPipeline(jobs_setup_file)
#        bam_file = join_path(self.data_dir,
#                             sample_id+"_dedup_reads.bam")
#        indels_target_intervals_file = join_path(self.data_dir,
#                                                 'indels_target_intervals.list')
#        pl.indels_realign(sample_id,
#                          bam_file=bam_file,
#                          indels_target_intervals_file=indels_target_intervals_file,
#                          )
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_base_recal_2(self):
#        """ test base recalibration with high quality score """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE)
#        sample_id = "hi_score_sample"
#        pl = GATKBPPipeline(jobs_setup_file)
#        bam_file = join_path(self.data_dir,
#                             sample_id+"_realigned_reads.bam")
#        pl.base_recal(sample_id, bam_file=bam_file)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_haplotype_caller_2(self):
#        """ test haplotype caller with big data size """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        sample_id = "test-big-sample"
#        dedup_reads_file = join_path(self.data_dir,
#                                     sample_id+"_recal_reads.bam")
#        self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                          sample_id,
#                                          "NA",
#                                          "NA",
#                                          "NA",
#                                          usage_mail="YES",
#                                          )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.hap_cal(sample_id, 
#                   bam_file=dedup_reads_file,
#                   )
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_combine_gvcfs_1(self):
#        """ test combine small gvcfs """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        gvcf_file_1 = join_path(self.data_dir,
#                                "test_M_sample_1.g.vcf")
#        self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                          "test_M_sample_1",
#                                          "NA",
#                                          "NA",
#                                          "NA",
#                                          )
#        gvcf_file_2 = join_path(self.data_dir,
#                                "test_M_sample_2.g.vcf")
#        self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                          "test_M_sample_2",
#                                          "NA",
#                                          "NA",
#                                          "NA",
#                                          )
#        gvcf_file_3 = join_path(self.data_dir,
#                                "test_M_sample_3.g.vcf")
#        self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                          "test_M_sample_3",
#                                          "NA",
#                                          "NA",
#                                          "NA",
#                                          )
#        gvcf_list = [gvcf_file_1, gvcf_file_2, gvcf_file_3]
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.combine_gvcfs(gvcf_list=gvcf_list)
#
#    @unittest.skip("reserved for load test")
#    def test_combine_gvcfs_2(self):
#        """ test combine gvcfs """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        job_name = self.test_function
#        gvcf_file_prefix = join_path(self.data_dir,
#                                     self.test_function+"_")
#        slurm_log_file = join_path(self.working_dir,
#                                   self.test_function+'.log')
#        gvcf_list = map(lambda x: gvcf_file_prefix+"{:03d}".format(x)+".gvcf",
#                        xrange(1,401))
#        ref = join_path(self.data_dir,
#                        "ref.fa")
#        out_merged_gvcf = join_path(self.working_dir,
#                                    self.test_function+"_merged.gvcf")
#        gatk_combine_gvcfs(job_name,
#                           FAST_PROJECT_CODE,
#                           slurm_log_file,
#                           gvcf_list,
#                           out_merged_gvcf,
#                           ref,
#                           self.working_dir,
#                           )
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_genotype_gvcfs(self):
#        """ test genotype gvcfs """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        gvcf_file = join_path(self.data_dir,
#                              "test_M_sample.g.vcf")
#        self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                          "test_M_sample",
#                                          "NA",
#                                          "NA",
#                                          "NA",
#                                          )
#        gvcf_list = [gvcf_file]
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.genotype_gvcfs(gvcf_list=gvcf_list)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_sample_1(self):
#        """ test sample pre-processing workflow (targeted sequencing) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        sample_id = "test-M_sample"
#        sample_usage_mail = {sample_id: "YES"}
#        targets_interval_list = join_path(self.data_dir,
#                                          "targets.interval_list") 
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
#                                                   #     sample_usage_mail=sample_usage_mail,
#                                                        targets_interval_list=targets_interval_list,
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_sample(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_sample_2(self):
#        """ test sample pre-processing workflow (whole exome) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(project_code='b2011097')
#        sample_id = "test_big_sample"
#        sample_group = self.current_func_name
#        fastq1_file = join_path(self.data_dir,
#                           sample_id + "_1.fastq.gz")
#        fastq2_file = join_path(self.data_dir,
#                           sample_id + "_2.fastq.gz")
#        self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                          sample_id,
#                                          fastq1_file,
#                                          fastq2_file,
#                                          sample_group,
#                                          )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_sample(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_sample_4(self):
#        """ test sample pre-processing workflow (Illumina WES) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        sample_id = "test-MR_sample"
#        sample_usage_mail = {sample_id: "YES"}
#        targets_interval_list = join_path(self.data_dir,
#                                          "targets.interval_list") 
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=SLOW_PROJECT_CODE,
#                                                        sample_usage_mail=sample_usage_mail,
##                                                        targets_interval_list=targets_interval_list,
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_sample(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_sample_5(self):
#        """ test sample pre-processing workflow (one of problematic sample from Illimina WES) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        sample_id = "test-Ill_prob"
#        sample_usage_mail = {sample_id: "YES"}
#        targets_interval_list = join_path(self.data_dir,
#                                          "targets.interval_list") 
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_sample(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_sample_Br_652(self):
#        """ test sample pre-processing workflow (one of problematic sample from Illimina WES) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        sample_id = "Br-652"
#        sample_usage_mail = {sample_id: "YES"}
#        targets_interval_list = join_path(self.data_dir,
#                                          "targets.interval_list") 
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_sample(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_sample_Br_724(self):
#        """ test sample pre-processing workflow (one of problematic sample from Illimina WES) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        sample_id = "Br-724"
#        sample_usage_mail = {sample_id: "YES"}
#        targets_interval_list = join_path(self.data_dir,
#                                          "targets.interval_list") 
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_sample(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_sample_Br_732(self):
#        """ test sample pre-processing workflow (one of problematic sample from Illimina WES) """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        sample_id = "Br-732"
#        sample_usage_mail = {sample_id: "YES"}
#        targets_interval_list = join_path(self.data_dir,
#                                          "targets.interval_list") 
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_sample(sample_id)
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_dataset_1(self):
#        """ test preprocess target sequencing dataset """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=FAST_PROJECT_CODE,
#                                                        dataset_usage_mail="YES",
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_dataset()
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_dataset_2(self):
#        """
#        test dataset pre-processing in the case that only one out of three
#        need to be processed
#        """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        targets_interval_list = join_path(self.data_dir,
#                                          "targets.interval_list") 
#        jobs_setup_file = self.__create_jobs_setup_file(targets_interval_list=targets_interval_list)
#        sample_group = self.current_func_name
#        for sample_idx in xrange(1,4):
#            sample_id = "test_sample_M_" + str(sample_idx)
#            fastq1_file = join_path(self.data_dir,
#                               sample_id + "_1.fastq.gz")
#            fastq2_file = join_path(self.data_dir,
#                               sample_id + "_2.fastq.gz")
#            if sample_idx > 2:
#                self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                                  sample_id,
#                                                  fastq1_file,
#                                                  fastq2_file,
#                                                  sample_group,
#                                                  )
#            else:
#                self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                                  sample_id,
#                                                  fastq1_file,
#                                                  fastq2_file,
#                                                  sample_group,
#                                                  preprocess_sample="NO",
#                                                  )
#        pl = GATKBPPipeline(jobs_setup_file)
#        samples_working_dir = join_path(self.working_dir,
#                                        "tmp")
#        for sample_idx in xrange(1,3):
#            sample_id = "test_sample_M_" + str(sample_idx)
#            sample_working_dir = join_path(samples_working_dir,
#                                           sample_id)
#            raw_vcf_source = join_path(self.data_dir,
#                                       sample_id+"_raw.g.vcf")
#            raw_vcf_dest = join_path(sample_working_dir,
#                                     sample_id+"_raw.g.vcf")
#            self.copy_file(raw_vcf_source, raw_vcf_dest)
#            raw_idx_source = join_path(self.data_dir,
#                                       sample_id+"_raw.g.vcf.idx")
#            raw_idx_dest = join_path(sample_working_dir,
#                                     sample_id+"_raw.g.vcf.idx")
#            self.copy_file(raw_idx_source, raw_idx_dest)
#        pl.preprocess_dataset()
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_preprocess_dataset_3(self):
#        """ test preprocess WES dataset """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file(project_code=SLOW_PROJECT_CODE,
#                                                        dataset_usage_mail="YES",
#                                                        )
#        pl = GATKBPPipeline(jobs_setup_file)
#        pl.preprocess_dataset()
#
#    @unittest.skipUnless(settings.SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
#    def test_garbage_collecting(self):
#        """ test if the workflow can correctly collecting garbage """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        sample_id = "test-sample"
#        self.__add_sample_jobs_setup_file(jobs_setup_file,
#                                          sample_id,
#                                          "NA",
#                                          "NA",
#                                          "NA",
#                                          )
#        pl = GATKBPPipeline(jobs_setup_file)
#        sample_rec = pl.samples[sample_id]
#        raw_aligned_source = join_path(self.data_dir,
#                                       sample_id+"_raw_aligned_reads.sam")
#        raw_aligned_dest = join_path(self.working_dir,
#                                     sample_rec.raw_aligned_reads_file)
#        self.copy_file(raw_aligned_source, raw_aligned_dest)
#        sorted_reads_source = join_path(self.data_dir,
#                                        sample_id+"_sorted_reads.bam")
#        sorted_reads_dest = join_path(self.working_dir,
#                                      sample_rec.sorted_reads_file)
#        self.copy_file(sorted_reads_source, sorted_reads_dest)
#        pl.mark_dup(sample_id)
#        pl.monitor_jobs()
#        self.assertFalse(path_exists(raw_aligned_dest),
#                         " garbage collecting process doesn't work properly")
#        self.assertFalse(path_exists(sorted_reads_dest),
#                         " garbage collecting process doesn't work properly")

    def tearDown(self):
        self.remove_working_dir()
