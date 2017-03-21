import sys
from collections import OrderedDict
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import set_log_file

PIPELINE_ALLOC_TIME_VAR = "pipeline_alloc_time"
PIPELINE_ALLOC_TIME_DFLT = "10-00:00:00"

def display_configs(func_name,
                    app_description,
                    kwargs,
                    pl,
                    custom_params=None,
                    ):
    log_file = set_log_file(kwargs['log_file'])
    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
    optional_params = None
    if log_file is not None:
        optional_params = OrderedDict()
        optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=app_description,
                     third_party_software_version=pl.get_third_party_software_version(),
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.disp_params_set("Pipeline parameters", pl.get_raw_obj_str())

    if hasattr(pl, "gatk_params"):
        disp.disp_params_set("GATK DNA-Seq Best Practice parameters", pl.gatk_params.get_raw_obj_str())
    if hasattr(pl, "plink_params"):
        disp.disp_params_set("Plink parameters", pl.plink_params.get_raw_obj_str())
    if hasattr(pl, "mutstat_params"):
        disp.disp_params_set("Mutation statistics database parameters", pl.mutstat_params.get_raw_obj_str())
    if hasattr(pl, "annovar_params"):
        disp.disp_params_set("Annovar parameters", pl.annovar_params.get_raw_obj_str())
    if hasattr(pl, "report_layout"):
        disp.disp_params_set("report layout parameters", pl.report_layout.get_raw_obj_str())
    if hasattr(pl, "rpt_params"):
        disp.disp_params_set("report parameters", pl.rpt_params.get_raw_obj_str())
    if custom_params is not None:
        disp.disp_params_set("custom parameters", custom_params)

    disp.new_section_txt(" . . . E X E C U T I N G . . . ")

def app_pycmm_slurm_monitor_pipeline(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pipeline = kwargs['pipeline']
    pipeline_description = kwargs['pipeline_description']
    pl = pipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    pipeline_description,
                    kwargs,
                    pl,
                    )
    # check if itself is in job mode 
    if len(pl.job_nodelist) != 0:
        pl.monitor_jobs()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_pipeline(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pipeline = kwargs['pipeline']
    pipeline_description = kwargs['pipeline_description']
    pl = pipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_configs(func_name,
                    pipeline_description,
                    kwargs,
                    pl,
                    )
    pipeline_bin = kwargs['pipeline_bin']
    if PIPELINE_ALLOC_TIME_VAR in kwargs:
        pipeline_alloc_time = kwargs[PIPELINE_ALLOC_TIME_VAR]
    else:
        pipeline_alloc_time = PIPELINE_ALLOC_TIME_DFLT
    if pl.project_code is not None:
        pl.run_slurm_monitor_pipeline(class_slurm_bin=pipeline_bin,
                                      alloc_time=pipeline_alloc_time,
                                      log_file=kwargs['log_file'],
                                      )
    else:
        pl.run_offline_pipeline()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
