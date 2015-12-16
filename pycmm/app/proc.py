import sys
from collections import OrderedDict
from pycmm.settings import CMMDB_TABLEANNOVAR_DESCRIPTION
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import set_log_file
from pycmm.proc.annovarlib import Annovar


def app_pycmm_dummy_table_annovar(*args, **kwargs):
    mylogger.getLogger(__name__)
    if 'log_file' in kwargs:
        log_file = set_log_file(kwargs['log_file'])
    else:
        log_file = None
    func_name = sys._getframe().f_code.co_name

    dataset_name = kwargs['dataset_name']
    input_file = kwargs['input_file']
    db_folder = kwargs['db_folder']
    buildver = kwargs['buildver']
    protocols = kwargs['protocols']
    operations = kwargs['operations']
    nastring = kwargs['nastring']
    data_out_folder = kwargs['data_out_folder']
    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = None
    required_params = OrderedDict()
    required_params['dataset name (--dataset_name)'] = dataset_name
    required_params['input file (--input_file)'] = input_file
    required_params['protocols (--protocols)'] = protocols
    required_params['operations (--operations)'] = operations
    required_params['data output folder (--data_out_folder)'] = data_out_folder
    optional_params = OrderedDict()
    optional_params['build version (--buildver)'] = buildver
    optional_params['NA string (--nastring)'] = nastring
    if log_file is not None:
        optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=CMMDB_TABLEANNOVAR_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    av = Annovar(dataset_name=dataset_name,
                 input_file=input_file,
                 db_folder=db_folder,
                 buildver=buildver,
                 protocols=protocols,
                 operations=operations,
                 nastring=nastring,
                 data_out_folder=data_out_folder,
                 )
    av.run_table_annovar()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
