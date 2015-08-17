import sys
import datetime
from collections import OrderedDict
from pycmm import settings
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import log_file_with_time_stamp
from pycmm.flow.gatkbp import GATKBPPipeline
from pycmm.flow.gatkbp import create_jobs_setup_file


def app_pycmm_dnaseq_pipeline(*args, **kwargs):
    mylogger.getLogger(__name__)
    time_stamp = datetime.datetime.now()
    log_file = log_file_with_time_stamp(kwargs["log_file"],
                                        time_stamp,
                                        )
    mylogger.set_log_file(log_file)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params=OrderedDict()
    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
    optional_params=OrderedDict()
    optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=settings.DNASEQ_PIPELINE_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    pl = GATKBPPipeline(kwargs['jobs_setup_file'])
    pl.preprocess_dataset()
    pl.monitor_jobs()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    time_stamp = datetime.datetime.now()
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params=OrderedDict()
    required_params['dataset name (-d)'] = kwargs['dataset_name']
    required_params['sample group (-g)'] = kwargs['sample_group']
    required_params['project code (-p)'] = kwargs['project_code']
    required_params['reference file (-R)'] = kwargs['reference_file']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    required_params['samples root directory (-I)'] = kwargs['samples_root_dir']
    optional_params=OrderedDict()
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
    disp.show_config(app_description=settings.DNASEQ_PIPELINE_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    create_jobs_setup_file(dataset_name=kwargs['dataset_name'],
                           sample_group=kwargs['sample_group'],
                           project_code=kwargs['project_code'],
                           reference_file=kwargs['reference_file'],
                           project_out_dir=kwargs['project_out_dir'],
                           samples_root_dir=kwargs['samples_root_dir'],
                           known_indels_file=kwargs['known_indels_file'],
                           dbsnp_file=kwargs['dbsnp_file'],
                           variants_calling=kwargs['variants_calling'],
                           targets_interval_list=kwargs['targets_interval_list'],
                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
