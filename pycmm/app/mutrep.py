import sys
from collections import OrderedDict
from pycmm.settings import MUTREP_CREATE_JOB_SETUP_FILE_DESCRIPTION
from pycmm.settings import MUTREP_FAMILY_REPORT_DESCRIPTION
from pycmm.settings import MUTREP_SUMMARY_REPORT_DESCRIPTION
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.mutrep import DNARegion
from pycmm.flow.mutrep import create_jobs_setup_file
from pycmm.app import display_configs

def app_pycmm_family_report(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = MutRepPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    fam_id = kwargs['fam_id']
    raw_report_regions = kwargs['report_regions']
    if raw_report_regions is None:
        report_regions = None
    else:
        report_regions = map(lambda x: DNARegion(x),
                             raw_report_regions.split(","))
    out_file = kwargs['out_file']

    custom_params = OrderedDict()
    custom_params['family id'] = fam_id
    if report_regions is None:
        custom_params['report region(s)'] = report_regions
    else:
        custom_params['report region(s)'] = map(lambda x: x.get_raw_obj_str(),
                                                report_regions)
    custom_params['output file'] = out_file
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    MUTREP_FAMILY_REPORT_DESCRIPTION,
                    kwargs,
                    pl,
                    custom_params=custom_params,
                    )
    pl.gen_family_report(fam_id, report_regions, out_file=out_file)
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_summary_report(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = MutRepPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    raw_report_regions = kwargs['report_regions']
    if raw_report_regions is None:
        report_regions = None
    else:
        report_regions = map(lambda x: DNARegion(x),
                             raw_report_regions.split(","))
    out_file = kwargs['out_file']

    custom_params = OrderedDict()
    if report_regions is None:
        custom_params['report region(s)'] = report_regions
    else:
        custom_params['report region(s)'] = map(lambda x: x.get_raw_obj_str(),
                                                report_regions)
    custom_params['output file'] = out_file
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    MUTREP_FAMILY_REPORT_DESCRIPTION,
                    kwargs,
                    pl,
                    custom_params=custom_params,
                    )
    pl.gen_summary_report(report_regions, out_file=out_file)
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_mutation_reports(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = MutRepPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    MUTREP_FAMILY_REPORT_DESCRIPTION,
                    kwargs,
                    pl,
                    )
    pl.monitor_jobs()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_mutrep_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['dataset name (-d)'] = kwargs['dataset_name']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    optional_params = OrderedDict()
    if kwargs['sample_info'] is not None:
        optional_params['sample information (-s)'] = kwargs['sample_info']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    if kwargs['job_alloc_time'] is not None:
        optional_params['job allocation time (--job_alloc_time)'] = kwargs['job_alloc_time']
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=MUTREP_CREATE_JOB_SETUP_FILE_DESCRIPTION,
                     third_party_software_version=None,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    layout_params = OrderedDict()
    if kwargs['anno_cols'] is not None:
        layout_params['annotation columns (-a)'] = kwargs['anno_cols']
    else:
        layout_params['annotation columns (-a)'] = "all"
    if kwargs['anno_excl_tags'] is not None:
        layout_params['annotation excluded tags (-E)'] = kwargs['anno_excl_tags'].split(",")
    if kwargs['rows_filter_actions'] is not None:
        layout_params['rows filtering criteria (--filter_actions)'] = kwargs['rows_filter_actions'].split(",")
    if kwargs['expression_patterns'] is not None:
        layout_params['expression_patterns (--expression_patterns)'] = kwargs['expression_patterns'].split(",")
    if kwargs['expression_usages'] is not None:
        layout_params['expression_usages (--expression_usages)'] = kwargs['expression_usages'].split(",")
    if kwargs['annotated_vcf_tabix'] is not None:
        layout_params['annotated vcf tablx file (-A)'] = kwargs['annotated_vcf_tabix']
    if kwargs['report_regions'] is not None:
        layout_params['report regions (-R)'] = kwargs['report_regions'].split(",")
    else:
        layout_params['report regions (-R)'] = "all"
    if kwargs['frequency_ratios'] is not None:
        layout_params['rare frequency ratios (-f)'] = OrderedDict(item.split(":") for item in kwargs['frequency_ratios'].split(","))
    if kwargs['filter_genes'] is not None:
        layout_params['filter genes (--filter_genes)'] = kwargs['filter_genes'].split(',')
    layout_params['split chromosome (--split_chrom)'] = kwargs['split_chrom']
    layout_params['summary_families sheet (--summary_families)'] = kwargs['summary_families_sheet']
    extra_anno_cols = {}
    extra_anno_cols['call detail (--call_detail)'] = kwargs['call_detail']
    extra_anno_cols['genotype quality (--call_gq)'] = kwargs['call_gq']
    layout_params['extra annotation columns'] = extra_anno_cols
    layout_params['only summary report (--only_summary)'] = kwargs['only_summary']
    layout_params['only families report (--only_families)'] = kwargs['only_families']
    layout_params['coloring shared variants (--coloring_shared)'] = kwargs['coloring_shared']
    layout_params['coloring variant zygosities (--coloring_zygosity)'] = kwargs['coloring_zygosity']
    disp.disp_params_set("Report layout parameters", layout_params)
    kwargs['project_name'] = kwargs['dataset_name']
    create_jobs_setup_file(**kwargs)
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
