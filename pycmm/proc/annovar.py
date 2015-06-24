import sys
from pycmm.template import pyCMMBase

ANNOVAR_PARAMS_INPUT_FILE_KEY = "input_file"
ANNOVAR_PARAMS_DB_FOLDER_KEY = "db_folder"
ANNOVAR_PARAMS_BUILDVER_KEY = "buildver"
ANNOVAR_PARAMS_OUT_PREFIX_KEY = "out_prefix"
ANNOVAR_PARAMS_DB_LIST_KEY = "db_list"
ANNOVAR_PARAMS_DB_NAME_KEY = "db_name"
ANNOVAR_PARAMS_DB_OP_KEY = "db_op"
ANNOVAR_PARAMS_NASTRING_KEY = "nastring"

class Annovar(pyCMMBase):
    """ A structure to parse and keep sample information """

    def __init__(self,
                 annovar_params={},
                 ):
        self.__annovar_params = annovar_params

    def table_annovar(self,
                      method,
                      params=None,
                      ):
        pass

