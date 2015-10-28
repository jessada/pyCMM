import sys
from pycmm.template import pyCMMBase
from pycmm.settings import LJB_SIFT_PREDICTION_COL_NAME as SIFT_PRED_COL
from pycmm.settings import LJB_POLYPHEN2_HDIV_PREDICTION_COL_NAME as POLYPHEN2_HDIV_PRED_COL
from pycmm.settings import LJB_POLYPHEN2_HVAR_PREDICTION_COL_NAME as POLYPHEN2_HVAR_PRED_COL
from pycmm.settings import LJB_LRT_PREDICTION_COL_NAME as LRT_PRED_COL
from pycmm.settings import LJB_MUTATIONTASTER_PREDICTION_COL_NAME as MUTATIONTASTER_PRED_COL
from pycmm.settings import LJB_MUTATIONASSESSOR_PREDICTION_COL_NAME as MUTATIONASSESSOR_PRED_COL
from pycmm.settings import LJB_FATHMM_PREDICTION_COL_NAME as FATHMM_PRED_COL
from pycmm.settings import LJB_RADIALSVM_PREDICTION_COL_NAME as RADIALSVM_PRED_COL
from pycmm.settings import LJB_LR_PREDICTION_COL_NAME as LR_PRED_COL

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

class PredictionTranslator(pyCMMBase):
    """
    A class to translate codes from effect predictors by using informaiton
    from http://annovar.openbioinformatics.org/en/latest/user-guide/filter/
    As 20 August 2015 the translation is as below

    SIFT (sift)                        D: Deleterious (sift<=0.05); T: tolerated (sift>0.tolerated05)
    PolyPhen 2 HDIV (pp2_hdiv)         D: Probably damaging (>=0.957),damaging P: possibly damaging (0.453<=pp2_hdiv<=0.956); B: benign (pp2_hdiv<=0.452)
    PolyPhen 2 HVar (pp2_hvar)         D: Probably damaging (>=0.909), P: possibly damaging (0.447<=pp2_hdiv<=0.909); B: benign (pp2_hdiv<=0.446)
    LRT (lrt)                          D: Deleterious; N: Neutral; U: Unknown
    MutationTaster (mt)                A" ("disease_causing_automatic"); "D" ("disease_causing"); "N" ("PolyPhenmorphism"); "P" ("polymorphism_automatic")
    MutationAssessor (ma)              H: high; M: medium; L: low; N: neutral. H/M means functional and L/N   means non-functional
    FATHMM (fathmm)                    D: Deleterious; T: Tolerated
    MetaSVM (metasvm)                  D: Deleterious; T: Tolerated
    MetaLR (metasvmalr)                D: Deleterious; T: Tolerated
    GERP++ (gerp++)                    higher scores are more deleterious
    PhyloP (phylop)                    higher scores are more deleterious
    SiPhy (siphy)                      higher scores are more deleterious
    """

    def __init__(self):
        self.__expl = self.__init_explanation()

    def __init_explanation(self):
        expl = {}
        expl[SIFT_PRED_COL] = {'D': 'Deleterious', 'T': 'Tolerated'}
        expl[POLYPHEN2_HDIV_PRED_COL] = {'D': 'Probably damaging', 'P': 'Possibly damaging', 'B': 'Benign'}
        expl[POLYPHEN2_HVAR_PRED_COL] = {'D': 'Probably damaging', 'P': 'Possibly damaging', 'B': 'Benign'}
        expl[LRT_PRED_COL] = {'D': 'Deleterious', 'N': 'Neutral', 'U': 'Unknown'}
        expl[MUTATIONTASTER_PRED_COL] = {'A': 'Disease causing automatic', 'D': 'Disease causing', 'N': 'PolyPhenmorphism', 'P': 'Polymorphism_automatic'}
        expl[MUTATIONASSESSOR_PRED_COL] = {'H': 'High', 'M': 'Medium', 'L': 'Low', 'N': 'Neutral', 'H/M': 'Functional', 'L/N': 'Non-functional'}
        expl[FATHMM_PRED_COL] = {'D': 'Deleterious', 'T': 'Tolerated'}
        expl[RADIALSVM_PRED_COL] = {'D': 'Deleterious', 'T': 'Tolerated'}
        expl[LR_PRED_COL] = {'D': 'Deleterious', 'T': 'Tolerated'}
        return expl

    def recognize(self, col_name):
        return col_name in self.__expl

    def describe(self, col_name, code):
        if code == ".":
            return "."
        if code in self.__expl[col_name]:
            return self.__expl[col_name][code]
        else:
            return "NA"
