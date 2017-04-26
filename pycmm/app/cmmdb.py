import sys
from collections import OrderedDict
from pycmm.settings import CMMDB_MUTSTAT_DESCRIPTION
from pycmm.settings import CMMDB_VCFAF2ANNOVAR_DESCRIPTION
from pycmm.settings import CMMDB_TABLEANNOVAR_DESCRIPTION
from pycmm.settings import CMMDB_CREATE_JOB_SETUP_FILE_DESCRIPTION
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import set_log_file
from pycmm.flow.cmmdb import CMMDBPipeline
from pycmm.flow.cmmdb import create_jobs_setup_file
from pycmm.app import display_configs


def app_pycmm_cmmdb_cal_mut_stat(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = CMMDBPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    CMMDB_MUTSTAT_DESCRIPTION,
                    kwargs,
                    pl,
                    )
    pl.cal_mut_stat()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_vcfaf_to_annovar(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = CMMDBPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    CMMDB_VCFAF2ANNOVAR_DESCRIPTION,
                    kwargs,
                    pl,
                    )
    pl.vcfaf_to_annovar()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_cmmdb_table_annovar(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = CMMDBPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    CMMDB_TABLEANNOVAR_DESCRIPTION,
                    kwargs,
                    pl,
                    )
    pl.table_annovar()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_cmmdb_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['dataset name (-d)'] = kwargs['dataset_name']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    required_params['vcf tabix file (-i)'] = kwargs['vcf_tabix_file']
    optional_params = OrderedDict()
    if kwargs['vcf2avdb_key_table'] is not None:
        optional_params['vcf2avdb key (-t)'] = kwargs['vcf2avdb_key_table']
    if kwargs['db_region'] is not None:
        optional_params['vcf region (-r)'] = kwargs['db_region']
    if kwargs['sample_info'] is not None:
        optional_params['sample information (-s)'] = kwargs['sample_info']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    if kwargs['job_alloc_time'] is not None:
        optional_params['job allocation time (--job_alloc_time)'] = kwargs['job_alloc_time']
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=CMMDB_CREATE_JOB_SETUP_FILE_DESCRIPTION,
                     third_party_software_version=None,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    kwargs['project_name'] = kwargs['dataset_name']
    create_jobs_setup_file(*args, **kwargs)
#    create_jobs_setup_file(project_name=kwargs['dataset_name'],
#                           project_out_dir=kwargs['project_out_dir'],
#                           vcf_tabix_file=kwargs['vcf_tabix_file'],
#                           db_region=kwargs['db_region'],
#                           sample_info=kwargs['sample_info'],
#                           project_code=kwargs['project_code'],
#                           job_alloc_time=kwargs['job_alloc_time'],
#                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
#                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
