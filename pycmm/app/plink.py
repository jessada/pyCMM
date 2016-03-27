import sys
from collections import OrderedDict
from pycmm.settings import PLINK_CREATE_JOB_SETUP_FILE_DESCRIPTION
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.flow.plink import PlinkPipeline
from pycmm.flow.plink import create_jobs_setup_file

#def __display_report_config(func_name,
#                            kwargs,
#                            pl,
#                            ):
#    log_file = set_log_file(kwargs['log_file'])
#    disp.new_section_txt("S T A R T <" + func_name + ">")
#    required_params = OrderedDict()
#    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
#    if 'fam_id' in kwargs:
#        required_params['family id (-f)'] = kwargs['fam_id']
#    if 'report_regions' in kwargs:
#        raw_report_regions = kwargs['report_regions']
#        if raw_report_regions is None:
#            required_params['report regions (-r)'] = "All" 
#        else:
#            required_params['report regions (-r)'] = raw_report_regions 
#    if 'out_file' in kwargs:
#        required_params['output file (-o)'] = kwargs['out_file']
#    optional_params = None
#    if log_file is not None:
#        optional_params = OrderedDict()
#        optional_params['log file (-l)'] = log_file
#    disp.show_config(app_description=settings.MUTREP_FAMILY_REPORT_DESCRIPTION,
#                     required_params=required_params,
#                     optional_params=optional_params,
#                     )
#    layout_params = OrderedDict()
#    layout_params['annotated columns'] = pl.report_layout.anno_cols
#    layout_params['annoated vcf tabix file'] = pl.annotated_vcf_tabix
#    report_regions = pl.report_layout.report_regions
#    if report_regions is not None:
#        report_region_outs = []
#        for report_region in report_regions:
#            report_region_out = report_region.chrom
#            if report_region.start_pos is not None:
#                report_region_out += ":" + report_region.start_pos
#                report_region_out += "-" + report_region.end_pos
#            report_region_outs.append(report_region_out)
#        layout_params['report regions'] = report_region_outs
#    else:
#        layout_params['report regions'] = "all"
#    layout_params['frequency ratios'] = pl.report_layout.freq_ratios
#    layout_params['reports output folder'] = pl.rpts_out_dir
#    layout_params['split chromosome'] = pl.report_layout.split_chrom
#    layout_params['summary_families sheet'] = pl.report_layout.summary_families_sheet
#    layout_params['call detail'] = pl.report_layout.call_detail
#    if pl.report_layout.anno_excl_tags:
#        layout_params['columns exclusion tags'] = pl.report_layout.anno_excl_tags
#    filter_actions_out = OrderedDict()
#    filter_actions_out['Rare'] = pl.report_layout.filter_rare
#    filter_actions_out['Non-Intergenic'] = pl.report_layout.filter_non_intergenic
#    filter_actions_out['Non-Intronic'] = pl.report_layout.filter_non_intronic
#    filter_actions_out['Non-Upstream'] = pl.report_layout.filter_non_upstream
#    filter_actions_out['Non-Downstream'] = pl.report_layout.filter_non_downtream
#    filter_actions_out['Non-UTR'] = pl.report_layout.filter_non_utr
#    filter_actions_out['Non-Synonymous'] = pl.report_layout.filter_non_synonymous
#    filter_actions_out['Has-mutation'] = pl.report_layout.filter_has_mutation
#    filter_actions_out['Has-shared'] = pl.report_layout.filter_has_shared
#    layout_params['filter actions'] = filter_actions_out
#    report_exclusion_out = {}
#    layout_params['report exclusion'] = report_exclusion_out
#    report_exclusion_out['only summary'] = pl.report_layout.only_summary
#    report_exclusion_out['only families'] = pl.report_layout.only_families
#    disp.disp_params_set("Report layout parameters", layout_params)
#    disp.new_section_txt(" . . . G E N E R A T I N G   R E P O R T S . . . ")

def app_pycmm_plink_hap_assocs_slurm(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = PlinkPipeline(kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
#    __display_report_config(func_name, kwargs, pl)
#    raw_report_regions = kwargs['report_regions']
#    out_file = kwargs['out_file']
#    if raw_report_regions is None:
#        report_regions = None
#    else:
#        report_regions = map(lambda x: ReportRegion(x), raw_report_regions.split(","))
#    pl.gen_family_report(fam_id, report_regions, out_file=out_file)
    if len(pl.job_nodelist) != 0:
        pl.monitor_jobs()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_plink_hap_assocs(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = PlinkPipeline(kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
#    __display_report_config(func_name, kwargs, pl)
#    raw_report_regions = kwargs['report_regions']
#    out_file = kwargs['out_file']
#    if raw_report_regions is None:
#        report_regions = None
#    else:
#        report_regions = map(lambda x: ReportRegion(x), raw_report_regions.split(","))
#    pl.gen_family_report(fam_id, report_regions, out_file=out_file)
    if pl.project_code is not None:
        pl.run_hap_assocs_slurm(jobs_setup_file=kwargs['jobs_setup_file'],
                                log_file=kwargs['log_file'])
    else:
        pl.run_hap_assocs_offline()
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
    required_params['phenotype file (--pheno)'] = kwargs['phenotype_file']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    optional_params['flow allocation time (--flow_alloc_time)'] = kwargs['flow_alloc_time']
    optional_params['cut-off p-value (-c)'] = kwargs['cutoff_pvalue']
    optional_params['haplotype windows (--hap_window)'] = kwargs['hap_window_sizes'].split(",")
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=PLINK_CREATE_JOB_SETUP_FILE_DESCRIPTION,
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
                           hap_window_sizes=kwargs['hap_window_sizes'],
                           project_code=kwargs['project_code'],
                           flow_alloc_time=kwargs['flow_alloc_time'],
                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
