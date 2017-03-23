import unittest
from os.path import exists as path_exists
from os.path import join as join_path
from pycmm.settings import SLURM_GATKBP_TEST
from pycmm.settings import FULL_SYSTEM_TEST
from pycmm.template import SafeTester
from pycmm.settings import TEST_PROJECT_CODE
from pycmm.flow.gatkbp import GATKBPPipeline
from pycmm.flow.test.test_gatkbp import create_gatk_jobs_setup_file


class TestGATKBPPipelineSlurm(SafeTester):

    def __init__(self, methodName):
        super(TestGATKBPPipelineSlurm, self).__init__(methodName=methodName,
                                                 test_module_name=__name__,
                                                 )

    def setUp(self):
        pass

    def __create_jobs_setup_file(self, *args, **kwargs):
        kwargs['project_code'] = TEST_PROJECT_CODE
        kwargs['project_name'] = self.test_function
        kwargs['project_out_dir'] = self.working_dir
        kwargs['samples_root_dir'] = self.data_dir
        return create_gatk_jobs_setup_file(*args, **kwargs)

# ************************************** test very small sequence w/o indel **************************************
    @unittest.skipUnless(SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_1(self):
        """ test sample pre-processing workflow (small) (w/o indel) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_id = "test-S_sample"
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(targets_interval_list=targets_interval_list,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        pl.preprocess_sample(sample_id)
# ************************************** test very small sequence w/o indel **************************************

# ************************************** test very small sequence with indel **************************************
    @unittest.skipUnless(SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_2(self):
        """ test sample pre-processing workflow (small with indel) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_id = "test-S_sample_w_indel"
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(targets_interval_list=targets_interval_list,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        pl.preprocess_sample(sample_id)

    @unittest.skipUnless(SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
    def test_process_dataset_2(self):
        """ test a complete workflow to preprocess a dataset (small with indel) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(dataset_usage_mail="YES",
                                                        targets_interval_list=targets_interval_list,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        pl.process_dataset()
# ************************************** test very small sequence with indel **************************************

# ************************************** test medium-size data (targeted sequenceing, AXEQ) **************************************
    @unittest.skipUnless(SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_3(self):
        """ test sample pre-processing workflow with targeted sequencing sample """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_id = "test-M_sample"
        targets_interval_list = join_path(self.data_dir,
                                          "targets.interval_list") 
        jobs_setup_file = self.__create_jobs_setup_file(targets_interval_list=targets_interval_list,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        pl.preprocess_sample(sample_id)
# ************************************** test medium-size data (targeted sequenceing, AXEQ) **************************************

# ************************************** test exome-size data **************************************
    @unittest.skipUnless(SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_4(self):
        """ test sample pre-processing workflow with exome-size sample """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_id = "test-L_sample"
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        pl.preprocess_sample(sample_id)
# ************************************** test exome-size data **************************************

# ************************************** test preprocess dataset without variant calling **************************************
    @unittest.skipUnless(SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
    def test_process_dataset_5(self):
        """ test a complete workflow to preprocess a dataset (small with indel) """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file(dataset_usage_mail="YES",
                                                        variants_calling=False,
                                                        )
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        pl.process_dataset()
# ************************************** test preprocess dataset without variant calling **************************************

# ************************************** test orginally wrong encode Phred score samples **************************************
    @unittest.skipUnless(SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
    def test_preprocess_sample_Br_652(self):
        """
        test sample pre-processing workflow (one of problematic samples from
        Illimina WES probably due to incompatible encoding in the original
        fastq files)
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        sample_id = "test_Br-652"
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        pl.preprocess_sample(sample_id)
# ************************************** test orginally wrong encode Phred score samples **************************************

    @unittest.skipUnless(SLURM_GATKBP_TEST, "taking too much UPPMAX cpu-core hours")
    def test_vqsr_02(self):
        """
        test VQSR process. Input are sampling of piper gvcfs
        """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        jobs_setup_file = self.__create_jobs_setup_file()
        pl = GATKBPPipeline(jobs_setup_file=jobs_setup_file)
        src_dir = join_path(self.data_dir,
                            'tmp')
        cmd = "cp -r"
        cmd += " " + src_dir
        cmd += " " + pl.project_out_dir
        self.exec_sh(cmd)
        pl.vqsr()

    def tearDown(self):
        self.remove_working_dir()
