from pycmm.utils import mylogger
import pkg_resources
import sys

param_display_fmt = "  {name:<45}{value}"

def new_section_txt(txt):
    adj_txt = " " + txt + " "
    mylogger.info("")
    mylogger.info(adj_txt.center(140,"*"))

def disp_header(header_txt):
    mylogger.info("")
    mylogger.info(header_txt)

def disp_subheader(subheader_txt):
    mylogger.info("  " + subheader_txt)

def disp_param(param_name, param_value):
    mylogger.info(param_display_fmt.format(name=param_name+":",
                                           value=param_value))

def debug_param(param_name, param_value):
    mylogger.debug(param_display_fmt.format(name=param_name+":",
                                            value=param_value))

def disp_subparam(subparam_name, subparam_value):
    disp_param("  "+subparam_name, subparam_value)

def show_config(app_description,
                     time_stamp,
                     required_params,
                     optional_params,
                     ):
    disp_header("Description")
    mylogger.info("  " + app_description)
    disp_header("Version and environment configuration")
    disp_param("pyCMM version", pkg_resources.get_distribution("pycmm").version)
    disp_param("parameters", " ".join(sys.argv[1:]))
    debug_param("debug mode", "ON")
    disp_param("time stamp", time_stamp)
    disp_params_set("Required parameters", required_params) 
    disp_params_set("Optional parameters", optional_params) 

def disp_params_set(params_name,
                    params,
                    ):
    disp_header(params_name)
    for key in params:
        disp_param(key, params[key])


