ENCODING_SOLEXA = "Solexa"
ENCODING_ILLUMINA_1_3 = "Illumina-1.3"
ENCODING_ILLUMINA_1_5 = "Illumina-1.5"
ENCODING_ILLUMINA_1_8 = "Illumina-1.8"

RANGES = {
#    'Sanger': (33, 93),
    ENCODING_SOLEXA: (59, 104),
    ENCODING_ILLUMINA_1_3: (64, 104),
    ENCODING_ILLUMINA_1_5: (67, 105),
    ENCODING_ILLUMINA_1_8: (33, 74),
}

def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89)
    """

    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)

def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings



class AnnovarParams(pyCMMBase):
    """ A structure to parse and keep sample information """

    def __init__(self,
                 params,
                 ):
        self.__params = params

    @property
    def input_file(self):
        return self.__params[ANNOVAR_PARAMS_INPUT_FILE_KEY]

    @property
    def db_folder(self):
        return self.__params[ANNOVAR_PARAMS_DB_FOLDER_KEY]

    @property
    def buildver(self):
        return self.__params[ANNOVAR_PARAMS_BUILDVER_KEY]

    @property
    def protocols(self):
        return self.__params[ANNOVAR_PARAMS_DB_NAMES_KEY]

    @property
    def operations(self):
        return self.__params[ANNOVAR_PARAMS_DB_OPS_KEY]

    @property
    def nastring(self):
        return self.__params[ANNOVAR_PARAMS_NASTRING_KEY]

class Annovar(pyCMMBase):
    """ A structure to parse and keep sample information """

    def __init__(self,
                 dataset_name,
                 input_file,
                 db_folder,
                 buildver,
                 protocols,
                 operations,
                 nastring,
                 data_out_folder,
                 ):
        self.__dataset_name = dataset_name
        self.__input_file = input_file
        self.__db_folder = db_folder
        self.__buildver = buildver
        self.__protocols = protocols
        self.__operations = operations
        self.__nastring = nastring
        self.__data_out_folder = data_out_folder
        self.__tmp_annovar_prefix = self.local_tmp_file

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["dataset name"] = self.dataset_name
        raw_repr["input file"] = self.input_file
        raw_repr["db folder"] = self.db_folder
        raw_repr["build ver"] = self.buildver
        raw_repr["protocols"] = self.protocols
        raw_repr["operations"] = self.operations
        raw_repr["NA string"] = self.nastring
        raw_repr["data out folder"] = self.data_out_folder
        raw_repr["tmp annovar prefix"] = self.tmp_annovar_prefix
        raw_repr["annotated vcf output file"] = self.out_annotated_vcf
        return raw_repr

    @property
    def dataset_name(self):
        return self.__dataset_name

    @property
    def input_file(self):
        return self.__input_file

    @property
    def db_folder(self):
        return self.__db_folder

    @property
    def buildver(self):
        return self.__buildver

    @property
    def protocols(self):
        return self.__protocols

    @property
    def operations(self):
        return self.__operations

    @property
    def nastring(self):
        return self.__nastring

    @property
    def data_out_folder(self):
        return self.__data_out_folder

    @property
    def tmp_annovar_prefix(self):
        return self.__tmp_annovar_prefix

    @property
    def out_annotated_vcf(self):
        file_name = self.dataset_name + "_annotated.vcf"
        return join_path(self.data_out_folder, file_name)

    @property
    def table_annovar_cmd(self):
        cmd = "table_annovar.pl"
        cmd += " " + self.input_file
        cmd += " " + self.db_folder
        cmd += " -buildver " + self.buildver
        cmd += " -out " + self.tmp_annovar_prefix
        cmd += " -remove"
        cmd += " -protocol " + self.protocols
        cmd += " -operation " + self.operations
        cmd += " -nastring " + self.nastring
        cmd += " -vcfinput"
        return cmd

    def __gen_annotated_vcf(self):
        return exec_sh(self.table_annovar_cmd)

    def run_table_annovar(self):
        # run table_annovar.pl, the result is in local tmp folder
        self.__gen_annotated_vcf()
        # copy results into the expected output file
        annotated_file = self.tmp_annovar_prefix + ".hg19_multianno.vcf" 
        self.copy_file(annotated_file, self.out_annotated_vcf)
        # delete tmp file
        self.remove_dir(get_path(annotated_file))
        self.info(" > > > D O N E < < < ")
        self.info("The result is at "+self.out_annotated_vcf)


class PredictionInfo(pyCMMBase):
    """
    To encapsulate the effect prediction information
    """

    def __init__(self,
                 code,
                 description,
                 harmful,
                 ):
        self.__code = code
        self.__description = description
        self.__harmful = harmful

    @property
    def code(self):
        return self.__code

    @property
    def description(self):
        if self.__code == ".":
            return self.__code
        if self.__description == ".":
            return "NA(" + self.__code + ")"
        return self.__description

    @property
    def harmful(self):
        return self.__harmful


class PredictionTranslator(pyCMMBase):
    """
    A class to translate codes from effect predictors by using informaiton
    from http://annovar.openbioinformatics.org/en/latest/user-guide/filter/
    As 20 August 2015 the translation is as below

    SIFT (sift)                        D: Deleterious (sift<=0.05); T: tolerated (sift>0.tolerated05)
    PolyPhen 2 HDIV (pp2_hdiv)         D: Probably damaging (>=0.957),damaging P: possibly damaging (0.453<=pp2_hdiv<=0.956); B: benign (pp2_hdiv<=0.452)
    PolyPhen 2 HVar (pp2_hvar)         D: Probably damaging (>=0.909), P: possibly damaging (0.447<=pp2_hdiv<=0.909); B: benign (pp2_hdiv<=0.446)
    LRT (lrt)                          D: Deleterious; N: Neutral; U: Unknown
    MutationTaster (mt)                A" ("disease_causing_automatic"); "D" ("disease_causing"); "N" ("polymorphism"); "P" ("polymorphism_automatic")
    MutationAssessor (ma)              H: high; M: medium; L: low; N: neutral. H/M means functional and L/N   means non-functional
    FATHMM (fathmm)                    D: Deleterious; T: Tolerated
    MetaSVM (metasvm)                  D: Deleterious; T: Tolerated
    MetaLR (metasvmalr)                D: Deleterious; T: Tolerated
    GERP++ (gerp++)                    higher scores are more deleterious
    PhyloP (phylop)                    higher scores are more deleterious
    SiPhy (siphy)                      higher scores are more deleterious
    """

    def __init__(self):
        self.__set_prediction_info()
        self.__null_prediction = PredictionInfo(code='.',
                                                description='.',
                                                harmful=False,
                                                )

    def __set_prediction_info(self):
        self.__pred_info = {}

        self.__pred_info[SIFT_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[SIFT_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='T',
                                   description='Tolerated',
                                   harmful=False,
                                   )
        self.__pred_info[SIFT_PRED_COL]['T'] = pred_info

        self.__pred_info[POLYPHEN2_HDIV_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Probably damaging',
                                   harmful=True,
                                   )
        self.__pred_info[POLYPHEN2_HDIV_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='P',
                                   description='Possibly damaging',
                                   harmful=True,
                                   )
        self.__pred_info[POLYPHEN2_HDIV_PRED_COL]['P'] = pred_info
        pred_info = PredictionInfo(code='B',
                                   description='Benign',
                                   harmful=False,
                                   )
        self.__pred_info[POLYPHEN2_HDIV_PRED_COL]['B'] = pred_info

        self.__pred_info[POLYPHEN2_HVAR_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Probably damaging',
                                   harmful=True,
                                   )
        self.__pred_info[POLYPHEN2_HVAR_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='P',
                                   description='Possibly damaging',
                                   harmful=True,
                                   )
        self.__pred_info[POLYPHEN2_HVAR_PRED_COL]['P'] = pred_info
        pred_info = PredictionInfo(code='B',
                                   description='Benign',
                                   harmful=False,
                                   )
        self.__pred_info[POLYPHEN2_HVAR_PRED_COL]['B'] = pred_info

        self.__pred_info[LRT_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[LRT_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='N',
                                   description='Neutral',
                                   harmful=False,
                                   )
        self.__pred_info[LRT_PRED_COL]['N'] = pred_info
        pred_info = PredictionInfo(code='U',
                                   description='Unknown',
                                   harmful=False,
                                   )
        self.__pred_info[LRT_PRED_COL]['U'] = pred_info

        self.__pred_info[MUTATIONTASTER_PRED_COL] = {}
        pred_info = PredictionInfo(code='A',
                                   description='Disease causing automatic',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONTASTER_PRED_COL]['A'] = pred_info
        pred_info = PredictionInfo(code='D',
                                   description='Disease causing',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONTASTER_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='N',
                                   description='Polymorphism',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONTASTER_PRED_COL]['N'] = pred_info
        pred_info = PredictionInfo(code='P',
                                   description='Polymorphism automatic',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONTASTER_PRED_COL]['P'] = pred_info

        self.__pred_info[MUTATIONASSESSOR_PRED_COL] = {}
        pred_info = PredictionInfo(code='H',
                                   description='High',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['H'] = pred_info
        pred_info = PredictionInfo(code='M',
                                   description='Medium',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['M'] = pred_info
        pred_info = PredictionInfo(code='L',
                                   description='Low',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['L'] = pred_info
        pred_info = PredictionInfo(code='N',
                                   description='Neutral',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['N'] = pred_info
        pred_info = PredictionInfo(code='H/M',
                                   description='Functional',
                                   harmful=True,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['H/M'] = pred_info
        pred_info = PredictionInfo(code='L/N',
                                   description='Non-functional',
                                   harmful=False,
                                   )
        self.__pred_info[MUTATIONASSESSOR_PRED_COL]['L/N'] = pred_info

        self.__pred_info[FATHMM_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[FATHMM_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='T',
                                   description='Tolerated',
                                   harmful=False,
                                   )
        self.__pred_info[FATHMM_PRED_COL]['T'] = pred_info

        self.__pred_info[RADIALSVM_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[RADIALSVM_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='T',
                                   description='Tolerated',
                                   harmful=False,
                                   )
        self.__pred_info[RADIALSVM_PRED_COL]['T'] = pred_info

        self.__pred_info[LR_PRED_COL] = {}
        pred_info = PredictionInfo(code='D',
                                   description='Deleterious',
                                   harmful=True,
                                   )
        self.__pred_info[LR_PRED_COL]['D'] = pred_info
        pred_info = PredictionInfo(code='T',
                                   description='Tolerated',
                                   harmful=False,
                                   )
        self.__pred_info[LR_PRED_COL]['T'] = pred_info

    def get_prediction_info(self,
                            predictor_name,
                            code,
                            ):
        # if no prediction
        if code == ".":
            return self.__null_prediction
        # if predictor is not yet in the system
        if predictor_name not in self.__pred_info.keys():
            warning_msg = "!! Unknown predictor '" + predictor_name + "'"
            self.warning(warning_msg)
            # create a dummy for the new predictor entity in the system
            self.__pred_info[predictor_name] = {}
            pred_info = PredictionInfo(code=code,
                                       description='.',
                                       harmful=False,
                                       )
            self.__pred_info[predictor_name][code] = pred_info
        # if prediction code is not yet in the system
        elif code not in self.__pred_info[predictor_name].keys():
            warning_msg = "!! Unknown prediction code '" + code + "'"
            warning_msg += " for predictor '" + predictor_name + "'"
            self.warning(warning_msg)
            # create a dummy for the new prediction code entity in the system
            pred_info = PredictionInfo(code=code,
                                       description='.',
                                       harmful=False,
                                       )
            self.__pred_info[predictor_name][code] = pred_info
        return self.__pred_info[predictor_name][code]

    @property
    def predictor_list(self):
        return self.__pred_info.keys()
