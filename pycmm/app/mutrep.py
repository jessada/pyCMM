import sys
import datetime
from collections import OrderedDict
from pycmm import settings
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import log_file_with_time_stamp
from pycmm.flow.mutrep import MutRepPipeline
from pycmm.flow.mutrep import ReportRegion


def app_pycmm_family_report(*args, **kwargs):
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
    fam_id = kwargs['fam_id']
    raw_report_regions = kwargs['report_regions']
    required_params['family id (-f)'] = fam_id
    required_params['report regions (-r)'] = raw_report_regions 
    optional_params = OrderedDict()
    optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=settings.CMMDB_MUTATIONREPORTS_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    pl = MutRepPipeline(kwargs['jobs_setup_file'])
    layout_params = OrderedDict()
    layout_params['annoated vcf tabix file'] = pl.annotated_vcf_tabix
    layout_params['annotated columns'] = pl.report_layout.anno_cols
    layout_params['extra call information'] = pl.report_layout.call_info
    layout_params['frequency ratios'] = pl.report_layout.freq_ratios
    layout_params['reports output folder'] = pl.rpts_out_dir
    mylogger.getLogger(__name__)
    disp.disp_params_set("Report layout parameters", layout_params)
    disp.new_section_txt(" . . . G E N E R A T I N G   R E P O R T S . . . ")
    report_regions = map(lambda x: ReportRegion(x), raw_report_regions.split(","))
    pl.gen_family_report(fam_id, report_regions)
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_mutation_reports(*args, **kwargs):
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
    disp.show_config(app_description=settings.CMMDB_MUTATIONREPORTS_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    pl = MutRepPipeline(kwargs['jobs_setup_file'])
    layout_params = OrderedDict()
    layout_params['annoated vcf tabix file'] = pl.annotated_vcf_tabix
    layout_params['annotated columns'] = pl.report_layout.anno_cols
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
    layout_params['extra call information'] = pl.report_layout.call_info
    layout_params['frequency ratios'] = pl.report_layout.freq_ratios
    layout_params['reports output folder'] = pl.rpts_out_dir
    mylogger.getLogger(__name__)
    disp.disp_params_set("Report layout parameters", layout_params)
    disp.new_section_txt(" . . . G E N E R A T I N G   R E P O R T S . . . ")
    pl.gen_family_reports()
    pl.gen_summary_report()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
