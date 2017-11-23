import sys
from pycmm.settings import PLINK_HAP_ASSOCS_REPORT_DESCRIPTION
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.flow.hapassocrep import HapAssocRepPipeline
from pycmm.app import display_configs

def app_pycmm_plink_hap_assocs_report(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = HapAssocRepPipeline(hap_assoc_file=kwargs['hap_assoc_file'],
                             snp_info_file=kwargs['snp_info_file'],
                             jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    PLINK_HAP_ASSOCS_REPORT_DESCRIPTION,
                    kwargs,
                    pl,
                    )
    pl.gen_report(out_file=kwargs['out_file'])
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
