import sys
from collections import OrderedDict
from pycmm import settings
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import set_log_file
from pycmm.flow.cmmdb import CMMDBPipeline
from pycmm.flow.cmmdb import create_jobs_setup_file


def app_pycmm_cmmdb_cal_mut_stat(*args, **kwargs):
    mylogger.getLogger(__name__)
    log_file = set_log_file(kwargs['log_file'])
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
    optional_params = None
    if log_file is not None:
        optional_params = OrderedDict()
        optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=settings.CMMDB_MUTSTAT_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    pl = CMMDBPipeline(kwargs['jobs_setup_file'])
    pl.cal_mut_stat()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_cmmdb_table_annovar(*args, **kwargs):
    mylogger.getLogger(__name__)
    log_file = set_log_file(kwargs['log_file'])
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
    optional_params = None
    if log_file is not None:
        optional_params = OrderedDict()
        optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=settings.CMMDB_TABLEANNOVAR_DESCRIPTION,
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
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['dataset name (-d)'] = kwargs['dataset_name']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    required_params['vcf tabix file (-i)'] = kwargs['vcf_tabix_file']
    optional_params = OrderedDict()
    if kwargs['db_region'] is not None:
        optional_params['vcf region (-r)'] = kwargs['db_region']
    if kwargs['sample_infos'] is not None:
        optional_params['sample information (-s)'] = kwargs['sample_infos']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    if kwargs['db_alloc_time'] is not None:
        optional_params['db allocation time (--db_alloc_time)'] = kwargs['db_alloc_time']
    if kwargs['rpt_alloc_time'] is not None:
        optional_params['report allocation time (--rpt_alloc_time)'] = kwargs['rpt_alloc_time']
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=settings.CMMDB_PIPELINE_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    layout_params = OrderedDict()
    if kwargs['anno_cols'] is not None:
        layout_params['annotation columns (-a)'] = kwargs['anno_cols']
    else:
        layout_params['annotation columns (-a)'] = "all"
    if kwargs['anno_excl_tags'] is not None:
        layout_params['annotation excluded tags (-a)'] = kwargs['anno_excl_tags']
    if kwargs['annotated_vcf_tabix'] is not None:
        layout_params['annotated vcf tablx file (-A)'] = kwargs['annotated_vcf_tabix']
    if kwargs['report_regions'] is not None:
        layout_params['report regions (-R)'] = kwargs['report_regions']
    else:
        layout_params['report regions (-R)'] = "all"
    if kwargs['frequency_ratios'] is not None:
        layout_params['rare frequency ratios (-f)'] = kwargs['frequency_ratios']
    layout_params['split chromosome (--split_chrom)'] = kwargs['split_chrom']
    layout_params['summary_families sheet (--summary_families)'] = kwargs['summary_families_sheet']
    extra_anno_cols = {}
    extra_anno_cols['call detail (--call_detail)'] = kwargs['call_detail']
    layout_params['extra annotation columns'] = extra_anno_cols
    exclusion_out = {}
    exclusion_out['common mutations (--exclude_common)'] = kwargs['exclude_common']
    exclusion_out['intergenic mutations (--exclude_intergenic)'] = kwargs['exclude_intergenic']
    exclusion_out['intronic mutations (--exclude_intronic)'] = kwargs['exclude_intronic']
    layout_params['exclusion'] = exclusion_out
    layout_params['only summary report (--only_summary)'] = kwargs['only_summary']
    layout_params['only families report (--only_families)'] = kwargs['only_families']
    disp.disp_params_set("Report layout parameters", layout_params)
    create_jobs_setup_file(dataset_name=kwargs['dataset_name'],
                           project_out_dir=kwargs['project_out_dir'],
                           vcf_tabix_file=kwargs['vcf_tabix_file'],
                           db_region=kwargs['db_region'],
                           sample_infos=kwargs['sample_infos'],
                           project_code=kwargs['project_code'],
                           db_alloc_time=kwargs['db_alloc_time'],
                           rpt_alloc_time=kwargs['rpt_alloc_time'],
                           anno_cols=kwargs['anno_cols'],
                           anno_excl_tags=kwargs['anno_excl_tags'],
                           annotated_vcf_tabix=kwargs['annotated_vcf_tabix'],
                           report_regions=kwargs['report_regions'],
                           frequency_ratios=kwargs['frequency_ratios'],
                           split_chrom=kwargs['split_chrom'],
                           summary_families_sheet=kwargs['summary_families_sheet'],
                           call_detail=kwargs['call_detail'],
                           exclude_common=kwargs['exclude_common'],
                           exclude_intergenic=kwargs['exclude_intergenic'],
                           exclude_intronic=kwargs['exclude_intronic'],
                           only_summary=kwargs['only_summary'],
                           only_families=kwargs['only_families'],
                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
