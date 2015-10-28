from pycmm.utils import mylogger
import pkg_resources
import sys

param_display_fmt = "  {name:<50}{value}"

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
                required_params,
                optional_params,
                ):
    disp_header("Description")
    mylogger.info("  " + app_description)
    disp_header("Version and environment configuration")
    disp_param("pyCMM version", pkg_resources.get_distribution("pycmm").version)
    disp_param("parameters", " ".join(sys.argv[1:]))
    debug_param("debug mode", "ON")
    disp_params_set("Required parameters", required_params) 
    if optional_params is not None:
        disp_params_set("Optional parameters", optional_params) 

def disp_params_set(params_name,
                    params,
                    ):
    disp_header(params_name)
    for key in params:
        val = params[key]
        if type(val) is list:
            disp_subheader(key)
            for entry_idx in xrange(len(val)):
                disp_subparam("#"+str(entry_idx+1), val[entry_idx])
        elif type(val) is dict:
            disp_subheader(key)
            for sub_key in val:
                disp_subparam(sub_key, val[sub_key])
        else:
            disp_param(key, params[key])

