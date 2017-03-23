import unittest
import filecmp
from os.path import exists as path_exists
from os.path import join as join_path
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.template import SafeTester
from pycmm.settings import TEST_PROJECT_CODE
from pycmm.flow.gatkbp import GATKBPPipeline
from pycmm.flow.gatkbp import create_jobs_setup_file
from pycmm.utils import file_size
from pycmm.utils import find_files
from pycmm.cmmlib.fastqlib import ENCODING_ILLUMINA_1_5
from pycmm.cmmlib.fastqlib import ENCODING_ILLUMINA_1_8

DFLT_KNOWN_INDELS = []
DFLT_KNOWN_INDELS.append('1000G_phase1.indels.b37.vcf')
DFLT_KNOWN_INDELS.append('Mills_and_1000G_gold_standard.indels.b37.vcf')
DFLT_DBSNP_FILE = 'dbsnp_138.b37.vcf'
DFLT_REF_FILE = 'ref.fa'

GATKBP_TEST = False


def create_gatk_jobs_setup_file(*args, **kwargs):
    known_indels_file = []
    for item in DFLT_KNOWN_INDELS:
        known_indels_file.append(join_path(kwargs['samples_root_dir'], item))
    kwargs['known_indels_file'] = known_indels_file
    kwargs['dbsnp_file'] = join_path(kwargs['samples_root_dir'], DFLT_DBSNP_FILE)
    kwargs['reference_file'] = join_path(kwargs['samples_root_dir'], DFLT_REF_FILE)
    kwargs['out_jobs_setup_file'] = join_path(kwargs['project_out_dir'],
                                              kwargs['project_name']+'_jobs_setup.txt')
    create_jobs_setup_file(*args, **kwargs)
    return kwargs['out_jobs_setup_file']

class TestGATKBPPipeline(SafeTester):

    def __init__(self, methodName):
        super(TestGATKBPPipeline, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self, *args, **kwargs):
        kwargs['project_name'] = self.test_function
        kwargs['project_out_dir'] = self.working_dir
        kwargs['samples_root_dir'] = self.data_dir
        return create_gatk_jobs_setup_file(*args, **kwargs)

    def test_load_jobs_info_1(self):
        """ test if default job configurations are loaded correctly """

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
        self.assertEqual(pl.job_alloc_time,
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
        self.assertEqual(pl.gatk_params.variants_calling,
                         False,
                         "GATKBPPipeline cannot correctly identify 'variants calling' from jobs setup file")
        self.assertEqual(pl.jobs_report_file,
                         join_path(pl.project_out_dir, exp_project_name+"_rpt.txt"),
                         "GATKBPPipeline cannot correctly identify 'jobs report file' from jobs setup file")
        self.assertEqual(pl.gatk_params.targets_interval_list,
                         None,
                         "GATKBPPipeline cannot correctly identify 'targets.interval_list' from jobs setup file")
        self.assertEqual(pl.gatk_params.split_regions_file,
                         None,
                         "GATKBPPipeline cannot correctly identify 'split regions list' from jobs setup file")
        self.assertEqual(pl.gatk_params.dataset_usage_mail,
                         False,
                         "GATKBPPipeline cannot correctly identify 'usage mail' from jobs setup file")
        self.assertEqual(len(pl.samples),
                         6,
                         "GATKBPPipeline cannot correctly identify number of samples from jobs setup file")
        test_sample_id_3 = 'test_sample_3'
        self.assertEqual(pl.samples[test_sample_id_3].sample_id,
                         test_sample_id_3,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples[test_sample_id_3].encoding,
                         ENCODING_ILLUMINA_1_8,
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
        self.assertEqual(pl.samples[test_sample_id_1].encoding,
                         ENCODING_ILLUMINA_1_8,
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
        test_sample_id_2 = 'test_sample_2'
        self.assertEqual(pl.samples[test_sample_id_2].encoding,
                         ENCODING_ILLUMINA_1_5,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(pl.samples[test_sample_id_2].preprocess_sample,
                         True,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        test_sample_id_4 = 'test_sample_4'
        self.assertEqual(pl.samples[test_sample_id_4].encoding,
                         ENCODING_ILLUMINA_1_8,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        test_sample_id_5 = 'test_sample_5'
        self.assertEqual(pl.samples[test_sample_id_5].encoding,
                         ENCODING_ILLUMINA_1_5,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        test_sample_id_6 = 'test_sample_6'
        self.assertEqual(pl.samples[test_sample_id_6].encoding,
                         ENCODING_ILLUMINA_1_8,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")

    def test_load_jobs_info_2(self):
        """ test if modified job configurations are loaded correctly """

        self.init_test(self.current_func_name)
        job_name = self.test_function
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(project_code=TEST_PROJECT_CODE,
                                                        variants_calling=True,
                                                        job_alloc_time="20:00:00",
                                                        targets_interval_list=targets_interval_list,
                                                        split_regions_file=split_regions_file,
                                                        preprocess_sample=False,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        exp_project_name = self.test_function
        self.assertEqual(pl.job_alloc_time,
                         "20:00:00",
                         "GATKBPPipeline cannot correctly identify 'gatk allocation time' from jobs setup file")
        self.assertEqual(pl.gatk_params.variants_calling,
                         True,
                         "GATKBPPipeline cannot correctly identify 'variants calling' from jobs setup file")
        self.assertEqual(pl.gatk_params.targets_interval_list,
                         targets_interval_list,
                         "GATKBPPipeline cannot correctly identify 'targets.interval_list' from jobs setup file")
        self.assertEqual(pl.gatk_params.split_regions_file,
                         split_regions_file,
                         "GATKBPPipeline cannot correctly identify 'split regions list' from jobs setup file")
        test_sample_id_2 = 'test_sample_2'
        self.assertEqual(pl.samples[test_sample_id_2].preprocess_sample,
                         False,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(len(pl.gatk_params.split_regions_txt_list),
                         3,
                         "GATKBPPipeline cannot correctly identify 'split regions list' from jobs setup file")
        self.assertEqual(pl.gatk_params.split_regions_txt_list[1],
                         '19_57742321_57742390',
                         "GATKBPPipeline cannot correctly identify 'split regions list' from jobs setup file")
        self.assertEqual(len(pl.samples[test_sample_id_2].fastq_pairs[0]),
                         2,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        test_ng_id_1 = 'test_negative_control_sample_1'
        self.assertEqual(len(pl.samples[test_ng_id_1].fastq_pairs[0]),
                         1,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(len(pl.samples[test_ng_id_1].fastq_pairs),
                         2,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        test_piper_id_1 = 'test_piper_sample_1'
        R1_file = pl.samples[test_piper_id_1].fastq_pairs[0]['R1']
        exp_R2_file = R1_file.replace('_1.fastq.gz', '_2.fastq.gz')
        self.assertEqual(pl.samples[test_piper_id_1].fastq_pairs[0]['R2'],
                         exp_R2_file,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")
        self.assertEqual(len(pl.samples[test_piper_id_1].fastq_pairs[0]),
                         2,
                         "GATKBPPipeline cannot correctly read sample information from jobs setup file")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_split_gvcfs_01(self):
        """
        test basic split gvcf produced by sample preprocessing.
        input are result from normal pairs of small fastq files
        """

        self.init_test(self.current_func_name)
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(split_regions_file=split_regions_file)
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        for sample_id in pl.samples_id:
            gvcf_file = join_path(join_path(self.data_dir,
                                            "gvcf"),
                                  sample_id+".g.vcf.gz")
            pl.copy_file(gvcf_file,
                         pl.gvcf_out_dir)
            tbi_file = join_path(join_path(self.data_dir,
                                           "gvcf"),
                                  sample_id+".g.vcf.gz.tbi")
            pl.copy_file(tbi_file,
                         pl.gvcf_out_dir)
        pl.split_gvcfs()
        region_txt = pl.gatk_params.split_regions_txt_list[0]
        split_gz = join_path(join_path(pl.split_gvcfs_out_dir,
                                       region_txt),
                             pl.samples_id[0]+"."+region_txt+".g.vcf.gz")
        self.assertEqual(file_size(split_gz),
                         4523,
                         "split_gvcfs doesn't function correctly")
        region_txt = pl.gatk_params.split_regions_txt_list[1]
        split_tbi = join_path(join_path(pl.split_gvcfs_out_dir,
                                        region_txt),
                              pl.samples_id[0]+"."+region_txt+".g.vcf.gz.tbi")
        self.assertTrue(file_size(split_tbi) > 10,
                        "split_gvcfs doesn't function correctly")
                                          
    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_split_gvcfs_02(self):
        """
        test basic split gvcf produced by piper.
        """

        self.init_test(self.current_func_name)
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(split_regions_file=split_regions_file)
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        for sample_id in pl.samples_id:
            gvcf_file = join_path(join_path(self.data_dir,
                                            "gvcf"),
                                  sample_id+".g.vcf.gz")
            pl.copy_file(gvcf_file,
                         pl.gvcf_out_dir)
            tbi_file = join_path(join_path(self.data_dir,
                                           "gvcf"),
                                  sample_id+".g.vcf.gz.tbi")
            pl.copy_file(tbi_file,
                         pl.gvcf_out_dir)
        pl.split_gvcfs()
        file_count = 0
        for file_name in find_files(pl.working_dir, pattern="*.g.vcf.gz"):
            file_count += 1
        self.assertEqual(file_count,
                         24,
                         "split_gvcfs doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_combine_gvcfs_01(self):
        """
        test combine split gvcf files. Input are result from normal pairs of
        small fastq files
        """

        self.init_test(self.current_func_name)
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(split_regions_file=split_regions_file)
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        src_dir = join_path(self.data_dir,
                            'split_gvcfs')
        cmd = "cp -r"
        cmd += " " + src_dir
        cmd += " " + pl.working_dir
        self.exec_sh(cmd)
        region_txt = pl.gatk_params.split_regions_txt_list[2]
        out_file = pl.combine_gvcfs(region_txt=region_txt)
        self.assertTrue(file_size(out_file) > 5000,
                        "combine gvcfs doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_combine_gvcfs_02(self):
        """
        test combine split gvcf files. Input are sampling of piper gvcfs
        """

        self.init_test(self.current_func_name)
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(split_regions_file=split_regions_file)
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        src_dir = join_path(self.data_dir,
                            'split_gvcfs')
        cmd = "cp -r"
        cmd += " " + src_dir
        cmd += " " + pl.working_dir
        self.exec_sh(cmd)
        for region_txt in pl.gatk_params.split_regions_txt_list:
            out_file = pl.combine_gvcfs(region_txt=region_txt)
            self.assertTrue(file_size(out_file) > 4000,
                            "combine gvcfs doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_genotype_gvcfs_01(self):
        """
        test GATK genotype gvcfs function. Input are result from normal pairs of
        small fastq files
        """

        self.init_test(self.current_func_name)
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(split_regions_file=split_regions_file)
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        src_dir = join_path(self.data_dir,
                            'tmp')
        cmd = "cp -r"
        cmd += " " + src_dir
        cmd += " " + pl.project_out_dir
        self.exec_sh(cmd)
        region_txt = pl.gatk_params.split_regions_txt_list[2]
        out_file = pl.genotype_gvcfs(region_txt=region_txt)
        self.assertTrue(file_size(out_file) > 5000,
                        "genotype gvcfs doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_genotype_gvcfs_02(self):
        """
        test GATK genotype gvcfs function. Input are sampling of piper gvcfs
        """

        self.init_test(self.current_func_name)
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(split_regions_file=split_regions_file)
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        src_dir = join_path(self.data_dir,
                            'tmp')
        cmd = "cp -r"
        cmd += " " + src_dir
        cmd += " " + pl.project_out_dir
        self.exec_sh(cmd)
        for region_txt in pl.gatk_params.split_regions_txt_list:
            out_file = pl.genotype_gvcfs(region_txt=region_txt)
            self.assertTrue(file_size(out_file) > 4000,
                            "genotype gvcfs doesn't function correctly")

    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_concat_vcfs_02(self):
        """
        test if split vcf cannbe concatenated. Input are sampling of piper gvcfs
        """

        self.init_test(self.current_func_name)
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(split_regions_file=split_regions_file)
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        src_dir = join_path(self.data_dir,
                            'tmp')
        cmd = "cp -r"
        cmd += " " + src_dir
        cmd += " " + pl.project_out_dir
        self.exec_sh(cmd)
        out_file = pl.concat_vcfs()
        self.assertTrue(file_size(out_file) > 4000,
                        "concat vcfs doesn't function correctly")

#    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
#    def test_vqsr_02(self):
#        """
#        test VQSR process. Input are sampling of piper gvcfs
#        """
#
#        self.individual_debug = True
#        self.init_test(self.current_func_name)
#        jobs_setup_file = self.__create_jobs_setup_file()
#        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
#        src_dir = join_path(self.data_dir,
#                            'tmp')
#        cmd = "cp -r"
#        cmd += " " + src_dir
#        cmd += " " + pl.project_out_dir
#        self.exec_sh(cmd)
#        out_file = pl.vqsr()
#        self.assertTrue(file_size(out_file) > 4000,
#                        "concat vcfs doesn't function correctly")
#
    @unittest.skipUnless(FULL_SYSTEM_TEST or GATKBP_TEST, "taking too long time to test")
    def test_preprocess_dataset_01(self):
        """
        test dataset processing without sample preprocessing steps.
        Input are result from normal pairs of small fastq files
        """

        self.init_test(self.current_func_name)
        split_regions_file = join_path(self.data_dir,
                                       "split.regions_file") 
        jobs_setup_file = self.__create_jobs_setup_file(split_regions_file=split_regions_file,
                                                        preprocess_sample=False,
                                                        variants_calling=True,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        for sample_id in pl.samples_id:
            gvcf_file = join_path(join_path(self.data_dir,
                                            "gvcf"),
                                  sample_id+".g.vcf.gz")
            pl.copy_file(gvcf_file,
                         pl.gvcf_out_dir)
            tbi_file = join_path(join_path(self.data_dir,
                                           "gvcf"),
                                  sample_id+".g.vcf.gz.tbi")
            pl.copy_file(tbi_file,
                         pl.gvcf_out_dir)
        pl.process_dataset()

    def tearDown(self):
        self.remove_working_dir()
