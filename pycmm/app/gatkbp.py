import sys
import datetime
from collections import OrderedDict
from pycmm import settings
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import log_file_with_time_stamp
from pycmm.flow.gatkbp import GATKBPPipeline


def app_pycmm_dnaseq_pipeline(*args, **kwargs):
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
    pl = GATKBPPipeline(kwargs['jobs_setup_file'])
    pl.preprocess_dataset()
    pl.monitor_jobs()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
