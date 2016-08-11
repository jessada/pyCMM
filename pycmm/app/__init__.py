import sys
from collections import OrderedDict
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import set_log_file

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
    disp.disp_params_set("Pipeline parameters", pl.get_raw_repr())

    if hasattr(pl, "gatk_params"):
## *********************************************************************************************** Need refactoring ***********************************************************************************************
        gatk_params = OrderedDict()
        gatk_params['reference file'] = pl.gatk_params.reference
        gatk_params['known indels'] = pl.gatk_params.known_indels
        gatk_params['dbsnp file'] = pl.gatk_params.dbsnp
        gatk_params['variants calling'] = pl.gatk_params.variants_calling
        gatk_params['targets interval list'] = pl.gatk_params.targets_interval_list
        gatk_params['split chromosome regions'] = pl.gatk_params.split_regions_file
        gatk_params['dataset usage mail'] = pl.gatk_params.dataset_usage_mail
        disp.disp_params_set("GATK DNA-Seq Best Practice parameters", gatk_params)
## *********************************************************************************************** Need refactoring ***********************************************************************************************
    if hasattr(pl, "plink_params"):
        disp.disp_params_set("Plink parameters", pl.plink_params.get_raw_repr())
    if hasattr(pl, "mutstat_params"):
        disp.disp_params_set("Mutation statistics database parameters", pl.mutstat_params.get_raw_repr())
    if hasattr(pl, "annovar_params"):
        disp.disp_params_set("Annovar parameters", pl.annovar_params.get_raw_repr())
    if hasattr(pl, "report_layout"):
        disp.disp_params_set("report layout parameters", pl.report_layout.get_raw_repr())
    if hasattr(pl, "rpt_params"):
        disp.disp_params_set("report parameters", pl.rpt_params.get_raw_repr())
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
    pl.monitor_jobs()
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
    if pl.project_code is not None:
        pl.run_slurm_monitor_pipeline(class_slurm_bin=pipeline_bin,
                                      log_file=kwargs['log_file']
                                      )
    else:
        pl.run_offline_pipeline()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
