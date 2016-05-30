import sys
from collections import OrderedDict
from pycmm.settings import DNASEQ_PIPELINE_DESCRIPTION
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.flow.gatkbp import create_jobs_setup_file

def app_pycmm_gatkbp_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['dataset name (-d)'] = kwargs['dataset_name']
    required_params['project code (-p)'] = kwargs['project_code']
    required_params['reference file (-R)'] = kwargs['reference_file']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    required_params['samples root directory (-I)'] = kwargs['samples_root_dir']
    optional_params = OrderedDict()
    if kwargs['gatkbp_alloc_time'] is not None:
        optional_params['GATK Best Practice allocation time (--gatkbp_alloc_time)'] = kwargs['gatkbp_alloc_time']
    if kwargs['known_indels_file'] is not None:
        known_indels_files = kwargs['known_indels_file']
        for idx in xrange(len(known_indels_files)):
            optional_params['known indels #'+str(idx+1)+' (--known_indels)'] = known_indels_files[idx]
    if kwargs['dbsnp_file'] is not None:
        optional_params['dbSNP (--dbsnp)'] = kwargs['dbsnp_file']
    optional_params['variants calling (--variants_calling)'] = kwargs['variants_calling']
    if kwargs['targets_interval_list'] is not None:
        optional_params['targets interval list (--targets_interval_list)'] = kwargs['targets_interval_list']
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=DNASEQ_PIPELINE_DESCRIPTION,
                     third_party_software_version=None,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    create_jobs_setup_file(project_name=kwargs['dataset_name'],
                           project_out_dir=kwargs['project_out_dir'],
                           reference_file=kwargs['reference_file'],
                           samples_root_dir=kwargs['samples_root_dir'],
                           project_code=kwargs['project_code'],
                           flow_alloc_time=kwargs['gatkbp_alloc_time'],
                           known_indels_file=kwargs['known_indels_file'],
                           dbsnp_file=kwargs['dbsnp_file'],
                           variants_calling=kwargs['variants_calling'],
                           targets_interval_list=kwargs['targets_interval_list'],
                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
