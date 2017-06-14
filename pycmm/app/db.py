import sys
from collections import OrderedDict
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.proc.db.dbms import create_dbms_jobs_setup_file
from pycmm.proc.db.dbms import SQLiteDBController
from pycmm.settings import DBMS_CREATE_JOB_SETUP_FILE_DESCRIPTION
from pycmm.settings import DBMS_EXECUTE_DB_JOBS_DESCRIPTION
from pycmm.app import display_configs

def app_pycmm_dbms_execute_db_jobs(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = SQLiteDBController(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    DBMS_EXECUTE_DB_JOBS_DESCRIPTION,
                    kwargs,
                    pl,
                    )
    pl.run_offline_pipeline()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_dbms_run_controller(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = SQLiteDBController(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    DBMS_EXECUTE_DB_JOBS_DESCRIPTION,
                    kwargs,
                    pl,
                    )
    pl.run_controller()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_dbms_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['project name (-d)'] = kwargs['project_name']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    optional_params = OrderedDict()
    if kwargs['sample_info'] is not None:
        optional_params['sample information (-s)'] = kwargs['sample_info']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    if kwargs['job_alloc_time'] is not None:
        optional_params['job allocation time (--job_alloc_time)'] = kwargs['job_alloc_time']
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=DBMS_CREATE_JOB_SETUP_FILE_DESCRIPTION,
                     third_party_software_version=None,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    db_params = OrderedDict()
#    if kwargs['anno_cols'] is not None:
#        layout_params['annotation columns (-a)'] = kwargs['anno_cols'].split(",")
#    else:
#        layout_params['annotation columns (-a)'] = "all"
#    if kwargs['anno_excl_tags'] is not None:
#        layout_params['annotation excluded tags (-E)'] = kwargs['anno_excl_tags'].split(",")
#    if kwargs['rows_filter_actions'] is not None:
#        layout_params['rows filtering criteria (--filter_actions)'] = kwargs['rows_filter_actions'].split(",")
#    if kwargs['expression_patterns'] is not None:
#        layout_params['expression_patterns (--expression_patterns)'] = kwargs['expression_patterns'].split(",")
#    if kwargs['expression_usages'] is not None:
#        layout_params['expression_usages (--expression_usages)'] = kwargs['expression_usages'].split(",")
#    if kwargs['annotated_vcf_tabix'] is not None:
#        layout_params['annotated vcf tablx file (-A)'] = kwargs['annotated_vcf_tabix']
#    if kwargs['report_regions'] is not None:
#        layout_params['report regions (-R)'] = kwargs['report_regions'].split(",")
#    else:
#        layout_params['report regions (-R)'] = "all"
#    if kwargs['frequency_ratios'] is not None:
#        layout_params['rare frequency ratios (-f)'] = OrderedDict(item.split(":") for item in kwargs['frequency_ratios'].split(","))
#    if kwargs['filter_genes'] is not None:
#        layout_params['filter genes (--filter_genes)'] = kwargs['filter_genes'].split(',')
#    if kwargs['color_genes'] is not None:
#        layout_params['color genes (--color_genes)'] = kwargs['color_genes'].split(',')
    db_params['db file (--db_fuke)'] = kwargs['db_file']
    db_params['db jobs (--db_jobs)'] = kwargs['db_jobs']
#    layout_params['summary_families sheet (--summary_families)'] = kwargs['summary_families_sheet']
#    extra_anno_cols = {}
#    extra_anno_cols['call detail (--call_detail)'] = kwargs['call_detail']
#    extra_anno_cols['genotype quality (--call_gq)'] = kwargs['call_gq']
#    layout_params['extra annotation columns'] = extra_anno_cols
#    layout_params['coloring shared variants (--coloring_shared)'] = kwargs['coloring_shared']
#    layout_params['coloring variant zygosities (--coloring_zygosity)'] = kwargs['coloring_zygosity']
#    layout_params['show shared variants (--show_shared_variants)'] = kwargs['show_shared_variants']
    disp.disp_params_set("database configuration parameters", db_params)
#    kwargs['project_name'] = kwargs['project_name']
    create_dbms_jobs_setup_file(**kwargs)
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
