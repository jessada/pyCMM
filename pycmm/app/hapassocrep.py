import sys
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.flow.hapassocrep import HapAssocRepPipeline

def app_pycmm_plink_hap_assocs_report(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = HapAssocRepPipeline(hap_assoc_file=kwargs['hap_assoc_file'],
                             snp_info_file=kwargs['snp_info_file'],
                             jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
#    __display_report_config(func_name, kwargs, pl)
#    raw_report_regions = kwargs['report_regions']
#    out_file = kwargs['out_file']
#    if raw_report_regions is None:
#        report_regions = None
#    else:
#        report_regions = map(lambda x: ReportRegion(x), raw_report_regions.split(","))
#    pl.gen_family_report(fam_id, report_regions, out_file=out_file)
    pl.gen_report()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
