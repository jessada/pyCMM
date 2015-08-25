import sys
import datetime
from collections import OrderedDict
from pycmm import settings
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import log_file_with_time_stamp
from pycmm.flow.cmmdb import CMMDBPipeline
from pycmm.flow.cmmdb import create_jobs_setup_file


def app_pycmm_cmmdb_cal_mut_stat(*args, **kwargs):
    mylogger.getLogger(__name__)
    time_stamp = datetime.datetime.now()
    log_file = log_file_with_time_stamp(kwargs["log_file"],
                                        time_stamp,
                                        )
    mylogger.set_log_file(log_file)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
    optional_params = OrderedDict()
    optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=settings.CMMDB_MUTSTAT_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    pl = CMMDBPipeline(kwargs['jobs_setup_file'])
    pl.cal_mut_stat()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_cmmdb_table_annovar(*args, **kwargs):
    mylogger.getLogger(__name__)
    time_stamp = datetime.datetime.now()
    log_file = log_file_with_time_stamp(kwargs["log_file"],
                                        time_stamp,
                                        )
    mylogger.set_log_file(log_file)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
    optional_params = OrderedDict()
    optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=settings.CMMDB_TABLEANNOVAR_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    pl = CMMDBPipeline(kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    ta_params = OrderedDict()
    ta_params['vcf tabix input file'] = pl.annovar_config.input_file
    ta_params['databases folder'] = pl.annovar_config.db_folder
    ta_params['build version'] = pl.annovar_config.buildver
    ta_params['output prefix'] = pl.annovar_config.out_prefix
    ta_params['protocols'] = pl.annovar_config.protocols
    ta_params['operations'] = pl.annovar_config.operations
    ta_params['NA string'] = pl.annovar_config.nastring
    ta_params['annotated vcf file'] = pl.annovar_config.annotated_vcf
    disp.disp_params_set("Table Annovar parameters", ta_params)
    pl.table_annovar()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_cmmdb_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    time_stamp = datetime.datetime.now()
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['dataset name (-d)'] = kwargs['dataset_name']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    required_params['vcf tabix file (-i)'] = kwargs['vcf_tabix_file']
    optional_params = OrderedDict()
    if kwargs['vcf_region'] is not None:
        optional_params['vcf region (-r)'] = kwargs['vcf_region']
    if kwargs['sample_infos'] is not None:
        optional_params['patients list (-c)'] = kwargs['sample_infos']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=settings.CMMDB_PIPELINE_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    layout_params = OrderedDict()
    if kwargs['anno_cols'] is not None:
        layout_params['annotation columns (-a)'] = kwargs['anno_cols']
    else:
        layout_params['annotation columns (-a)'] = "all"
    if kwargs['annotated_vcf_tabix'] is not None:
        layout_params['annotated vcf tablx file (-A)'] = kwargs['annotated_vcf_tabix']
    if kwargs['report_regions'] is not None:
        layout_params['report regions (-R)'] = kwargs['report_regions']
    else:
        layout_params['report regions (-R)'] = "all"
    layout_params['call information (-C)'] = kwargs['call_info']
    if kwargs['frequency_ratios'] is not None:
        layout_params['rare frequency ratios (-f)'] = kwargs['frequency_ratios']
    layout_params['split chromosome (-s)'] = kwargs['split_chrom']
    disp.disp_params_set("Report layout parameters", layout_params)
    create_jobs_setup_file(dataset_name=kwargs['dataset_name'],
                           project_out_dir=kwargs['project_out_dir'],
                           vcf_tabix_file=kwargs['vcf_tabix_file'],
                           vcf_region=kwargs['vcf_region'],
                           sample_infos=kwargs['sample_infos'],
                           project_code=kwargs['project_code'],
                           anno_cols=kwargs['anno_cols'],
                           annotated_vcf_tabix=kwargs['annotated_vcf_tabix'],
                           report_regions=kwargs['report_regions'],
                           call_info=kwargs['call_info'],
                           frequency_ratios=kwargs['frequency_ratios'],
                           split_chrom=kwargs['split_chrom'],
                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
