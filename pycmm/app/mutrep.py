import sys
import datetime
from collections import OrderedDict
from pycmm import settings
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import log_file_with_time_stamp
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.mutrep import create_jobs_setup_file


def app_pycmm_mutrep_cal_mut_stat(*args, **kwargs):
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
    pl = MutRepPipeline(kwargs['jobs_setup_file'])
    pl.cal_mut_stat()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_mutrep_table_annovar(*args, **kwargs):
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
    disp.show_config(app_description=settings.MUTREP_PIPELINE_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    pl = MutRepPipeline(kwargs['jobs_setup_file'])
    pl.table_annovar()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_mutrep_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    time_stamp = datetime.datetime.now()
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params=OrderedDict()
    required_params['dataset name (-d)'] = kwargs['dataset_name']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    required_params['vcf tabix file (-i)'] = kwargs['vcf_tabix_file']
    optional_params=OrderedDict()
    if kwargs['vcf_region'] is not None:
        optional_params['vcf region (-r)'] = kwargs['vcf_region']
    if kwargs['patients_list'] is not None:
        optional_params['patients list (-c)'] = kwargs['patients_list']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=settings.DNASEQ_PIPELINE_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    create_jobs_setup_file(dataset_name=kwargs['dataset_name'],
                           project_out_dir=kwargs['project_out_dir'],
                           vcf_tabix_file=kwargs['vcf_tabix_file'],
                           vcf_region=kwargs['vcf_region'],
                           patients_list=kwargs['patients_list'],
                           project_code=kwargs['project_code'],
                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
