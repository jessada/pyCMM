import sys
from pycmm.template import pyCMMBase

ANNOVAR_PARAMS_INPUT_FILE_KEY = "input_file"
ANNOVAR_PARAMS_DB_FOLDER_KEY = "db_folder"
ANNOVAR_PARAMS_BUILDVER_KEY = "buildver"
ANNOVAR_PARAMS_OUT_PREFIX_KEY = "out_prefix"
ANNOVAR_PARAMS_DB_NAMES_KEY = "db_names"
ANNOVAR_PARAMS_DB_OPS_KEY = "db_ops"
ANNOVAR_PARAMS_NASTRING_KEY = "nastring"

class Annovar(pyCMMBase):
    """ A structure to parse and keep sample information """

    def __init__(self,
                 annovar_params,
                 ):
        self.__annovar_params = annovar_params

    @property
    def input_file(self):
        return self.__annovar_params[ANNOVAR_PARAMS_INPUT_FILE_KEY]

    @property
    def db_folder(self):
        return self.__annovar_params[ANNOVAR_PARAMS_DB_FOLDER_KEY]

    @property
    def buildver(self):
        return self.__annovar_params[ANNOVAR_PARAMS_BUILDVER_KEY]

    @property
    def out_prefix(self):
        return self.__annovar_params[ANNOVAR_PARAMS_OUT_PREFIX_KEY]

    @property
    def protocols(self):
        return self.__annovar_params[ANNOVAR_PARAMS_DB_NAMES_KEY]

    @property
    def operations(self):
        return self.__annovar_params[ANNOVAR_PARAMS_DB_OPS_KEY]

    @property
    def nastring(self):
        return self.__annovar_params[ANNOVAR_PARAMS_NASTRING_KEY]

    @property
    def annotated_vcf(self):
        return self.out_prefix + "." + self.buildver + "_multianno.vcf"

    @property
    def table_annovar_cmd(self):
        cmd = "table_annovar.pl"
        cmd += " " + self.input_file
        cmd += " " + self.db_folder
        cmd += " -buildver " + self.buildver
        cmd += " -out " + self.out_prefix
        cmd += " -remove"
        cmd += " -protocol " + self.protocols
        cmd += " -operation " + self.operations
        cmd += " -nastring " + self.nastring
        cmd += " -vcfinput"
        return cmd
