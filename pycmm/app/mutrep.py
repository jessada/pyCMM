import sys
from collections import OrderedDict
from pycmm import settings
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import set_log_file
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.mutrep import ReportRegion

def __display_report_config(func_name,
                            kwargs,
                            pl,
                            ):
    log_file = set_log_file(kwargs['log_file'])
    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
    if 'fam_id' in kwargs:
        required_params['family id (-f)'] = kwargs['fam_id']
    if 'report_regions' in kwargs:
        raw_report_regions = kwargs['report_regions']
        if raw_report_regions is None:
            required_params['report regions (-r)'] = "All" 
        else:
            required_params['report regions (-r)'] = raw_report_regions 
    if 'out_file' in kwargs:
        required_params['output file (-o)'] = kwargs['out_file']
    optional_params = None
    if log_file is not None:
        optional_params = OrderedDict()
        optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=settings.MUTREP_FAMILY_REPORT_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    layout_params = OrderedDict()
    layout_params['annotated columns'] = pl.report_layout.anno_cols
    layout_params['annoated vcf tabix file'] = pl.annotated_vcf_tabix
    report_regions = pl.report_layout.report_regions
    if report_regions is not None:
        report_region_outs = []
        for report_region in report_regions:
            report_region_out = report_region.chrom
            if report_region.start_pos is not None:
                report_region_out += ":" + report_region.start_pos
                report_region_out += "-" + report_region.end_pos
            report_region_outs.append(report_region_out)
        layout_params['report regions'] = report_region_outs
    else:
        layout_params['report regions'] = "all"
    layout_params['frequency ratios'] = pl.report_layout.freq_ratios
    layout_params['reports output folder'] = pl.rpts_out_dir
    layout_params['split chromosome'] = pl.report_layout.split_chrom
    layout_params['summary_families sheet'] = pl.report_layout.summary_families_sheet
    layout_params['call detail'] = pl.report_layout.call_detail
    if pl.report_layout.anno_excl_tags:
        layout_params['columns exclusion tags'] = pl.report_layout.anno_excl_tags
    filter_actions_out = {}
    filter_actions_out['Rare'] = pl.report_layout.filter_rare
    filter_actions_out['Non-Intergenic'] = pl.report_layout.filter_non_intergenic
    filter_actions_out['Non-Intronic'] = pl.report_layout.filter_non_intronic
    filter_actions_out['Has-mutation'] = pl.report_layout.filter_has_mutation
    filter_actions_out['Has-shared'] = pl.report_layout.filter_has_shared
    layout_params['filter actions'] = filter_actions_out
    report_exclusion_out = {}
    layout_params['report exclusion'] = report_exclusion_out
    report_exclusion_out['only summary'] = pl.report_layout.only_summary
    report_exclusion_out['only families'] = pl.report_layout.only_families
    disp.disp_params_set("Report layout parameters", layout_params)
    disp.new_section_txt(" . . . G E N E R A T I N G   R E P O R T S . . . ")

def app_pycmm_family_report(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = MutRepPipeline(kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    __display_report_config(func_name, kwargs, pl)
    fam_id = kwargs['fam_id']
    raw_report_regions = kwargs['report_regions']
    out_file = kwargs['out_file']
    if raw_report_regions is None:
        report_regions = None
    else:
        report_regions = map(lambda x: ReportRegion(x), raw_report_regions.split(","))
    pl.gen_family_report(fam_id, report_regions, out_file=out_file)
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_summary_report(*args, **kwargs):
    mylogger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name
    pl = MutRepPipeline(kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    __display_report_config(func_name, kwargs, pl)
    raw_report_regions = kwargs['report_regions']
    out_file = kwargs['out_file']
    if raw_report_regions is None:
        report_regions = None
    else:
        report_regions = map(lambda x: ReportRegion(x), raw_report_regions.split(","))
    pl.gen_summary_report(report_regions, out_file=out_file)
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_mutation_reports(*args, **kwargs):
    mylogger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name
    pl = MutRepPipeline(kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    __display_report_config(func_name, kwargs, pl)
    if pl.report_layout.only_summary:
        pl.gen_summary_reports()
    elif pl.report_layout.only_families:
        pl.gen_families_reports()
    else:
        pl.gen_families_reports()
        pl.gen_summary_reports()
    pl.monitor_jobs()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
