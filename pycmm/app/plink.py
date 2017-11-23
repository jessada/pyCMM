import sys
from collections import OrderedDict
from pycmm.settings import PLINK_CREATE_JOB_SETUP_FILE_DESCRIPTION
from pycmm.settings import PLINK_MERGE_HAP_ASSOCS_DESCRIPTION
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import set_log_file
from pycmm.app import display_configs
from pycmm.flow.plink import PlinkPipeline
from pycmm.flow.plink import create_jobs_setup_file

def app_pycmm_plink_merge_hap_assocs(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = PlinkPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    PLINK_MERGE_HAP_ASSOCS_DESCRIPTION,
                    kwargs,
                    pl,
                    )
    pl.merge_hap_assocs(hap_assoc_files=kwargs['hap_assoc_files'],
                        out_file=kwargs['out_file']
                        )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_plink_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['project name (-d)'] = kwargs['project_name']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    input_binary = kwargs['input_binary']
    input_file_prefix = kwargs['input_file_prefix']
    if kwargs['input_binary']:
        required_params['input file prefix (--bfile)'] = input_file_prefix
    else:
        required_params['input file prefix (--file)'] = input_file_prefix
    optional_params = OrderedDict()
    if kwargs['input_dna_regions'] is not None:
        optional_params['input dna regions (-R)'] = kwargs['input_dna_regions'].split(",")
    if kwargs['filter_criteria'] is not None:
        optional_params['filter criteria (--filter_criteria)'] = kwargs['filter_criteria'].split(",")
    required_params['phenotype file (--pheno)'] = kwargs['phenotype_file']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    if kwargs['sample_info'] is not None:
        optional_params['sample information (-p)'] = kwargs['sample_info'].split(",")
    optional_params['job allocation time (--job_alloc_time)'] = kwargs['job_alloc_time']
    optional_params['cut-off p-value (--cutoff_pvalue)'] = kwargs['cutoff_pvalue']
    optional_params['cut-off odds ratio (--cutoff_ors)'] = kwargs['cutoff_ors']
    optional_params['families haplotyep file prefix (--fam_hap)'] = kwargs['fam_hap_prefix']
    optional_params['haplotype windows (--hap_window)'] = kwargs['hap_window_sizes'].split(",")
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=PLINK_CREATE_JOB_SETUP_FILE_DESCRIPTION,
                     third_party_software_version=None,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    create_jobs_setup_file(project_name=kwargs['project_name'],
                           project_out_dir=kwargs['project_out_dir'],
                           input_file_prefix=kwargs['input_file_prefix'],
                           phenotype_file=kwargs['phenotype_file'],
                           input_binary=kwargs['input_binary'],
                           input_dna_regions=kwargs['input_dna_regions'],
                           cutoff_pvalue=kwargs['cutoff_pvalue'],
                           cutoff_ors=kwargs['cutoff_ors'],
                           hap_window_sizes=kwargs['hap_window_sizes'],
                           project_code=kwargs['project_code'],
                           sample_info=kwargs['sample_info'],
                           job_alloc_time=kwargs['job_alloc_time'],
                           filter_criteria=kwargs['filter_criteria'],
                           fam_hap_prefix=kwargs['fam_hap_prefix'],
                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
