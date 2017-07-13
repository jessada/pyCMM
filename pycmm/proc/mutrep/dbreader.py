import copy
import re
from collections import OrderedDict
from pycmm.settings import MAX_REF_MAF_COL_NAME
from pycmm.settings import REF_MAF_COL_NAMES
from pycmm.settings import EXAC03_CONSTRAINT_COL_NAMES
from pycmm.settings import EXAC03_CONSTRAINT_COL_NAME
from pycmm.settings import INTERVAR_AND_EVIDENCE_COL_NAME
from pycmm.settings import INTERVAR_CLASS_COL_NAME
from pycmm.settings import INTERVAR_EVIDENCE_COL_NAME
from pycmm.settings import PATHOGENIC_COUNT_COL_NAME
from pycmm.settings import PREDICTION_COLS
from pycmm.settings import LJB_METASVM_PREDICTION_COL_NAME
from pycmm.settings import LJB_METALR_PREDICTION_COL_NAME
from pycmm.settings import LJB_RADIALSVM_PREDICTION_COL_NAME
from pycmm.settings import LJB_LR_PREDICTION_COL_NAME
from pycmm.utils import is_number
from pycmm.template import pyCMMBase
from pycmm.proc.db.connector import SQLiteDB
from pycmm.proc.db.connector import VCF_PKEYS
from pycmm.proc.db.connector import TBL_NAME_ALL_GTZ_ANNOS
from pycmm.proc.db.connector import REF_MUTATED_COL_NAME
from pycmm.cmmlib.intervarlib import parse_intervar_class
from pycmm.cmmlib.intervarlib import parse_intervar_evidence
from pycmm.cmmlib.annovarlib import PredictionTranslator

MASTER_TABLE = TBL_NAME_ALL_GTZ_ANNOS

GT_WT = 'wt'
GT_HOM = 'hom'
GT_HET = 'het'
GT_OTH = 'oth'

FUNC_INTERGENIC = 'intergenic'
FUNC_INTRONIC = 'intronic'
FUNC_UPSTREAM = 'upstream'
FUNC_DOWNSTREAM = 'downstream'
FUNC_UTR = 'UTR'
EXONICFUNC_SYNONYMOUS = 'synonymous_SNV'

QRY_CHROM_COL_IDX = 0
QRY_POS_COL_IDX = 1
QRY_REF_COL_IDX = 2
QRY_ALT_COL_IDX = 3
#QRY_FILTER_COL_IDX = 4

EXAC_CONSTRAINT_PATTERN = re.compile(r'''(?P<var_name>.+?)=(?P<value>.+)''')

PREDICTION_COLS = list(PREDICTION_COLS)
PREDICTION_COLS.append(LJB_METASVM_PREDICTION_COL_NAME)
PREDICTION_COLS.append(LJB_METALR_PREDICTION_COL_NAME)
PREDICTION_COLS.remove(LJB_RADIALSVM_PREDICTION_COL_NAME)
PREDICTION_COLS.remove(LJB_LR_PREDICTION_COL_NAME)

class QryCall(pyCMMBase):
    """
    A structure to parse a calling info for a sample. Given raw genotype and
    alle frequency, it must be able to identify 
    - if a sample is actually mutated
    - what is the actual genotype compared to the whole population
    There is only one gt for one sample, contrary to how vcf store the data
    """

    def __init__(self,
                 *args,
                 **kwargs
                 ):
        super(QryCall, self).__init__(*args, **kwargs)
        self.__shared_mutation = False

    def get_raw_obj_str(self, *args, **kwargs):
        raw_str = OrderedDict()
        raw_str["raw zygosity"] = self.raw_zygo
        return raw_str

    @property
    def raw_zygo(self):
        return self.__raw_zygo

    @raw_zygo.setter
    def raw_zygo(self, val):
        self.__raw_zygo = val

    @property
    def actual_zygo(self):
        return self.__actual_zygo

    @actual_zygo.setter
    def actual_zygo(self, val):
        self.__actual_zygo = val

    @property
    def mutated(self):
        return self.__mutated

    @mutated.setter
    def mutated(self, val):
        self.__mutated = val

    @property
    def shared_mutation(self):
        return self.__shared_mutation

    @shared_mutation.setter
    def shared_mutation(self, val):
        self.__shared_mutation = val

class QryRecord(pyCMMBase):
    """
    A structure to parse a query record from SQLiteDB to be used by
    mutrep module. Three basic functions are
    - can get annotation value via .get_anno(anno_col)
    - can get call value call value via .get_call(sample_id)
    - can evaluate expression through .qry_eval(expr)
    """

    def __init__(self,
                 raw_qry,
                 anno_col_idxs,
                 sample_col_idxs,
                 family_infos,
                 pred_tran,
                 *args,
                 **kwargs
                 ):
        super(QryRecord, self).__init__(*args, **kwargs)
        self.__raw_qry = raw_qry
        self.__anno_col_idxs = anno_col_idxs
        self.__sample_col_idxs = sample_col_idxs
        self.__family_infos = copy.deepcopy(family_infos)
        self.__pred_tran = pred_tran
        self.__exac_constraints = -1
        self.__max_ref_maf = None
        self.__pathogenic_count = None
        self.__parse_genotype()
        self.__cal_shared()

    def get_raw_obj_str(self, *args, **kwargs):
        raw_str = OrderedDict()
        raw_str["CHROM"] = self.CHROM
        raw_str["POS"] = self.POS
        raw_str["REF"] = self.REF
        raw_str["ALT"] = self.ALT
        return raw_str

    @property
    def CHROM(self):
        chrom = self.__raw_qry[QRY_CHROM_COL_IDX]
        if is_number(chrom):
            chrom = int(chrom)
        return chrom

    @property
    def POS(self):
        return self.__raw_qry[QRY_POS_COL_IDX]

    @property
    def REF(self):
        return self.__raw_qry[QRY_REF_COL_IDX]

    @property
    def ALT(self):
        return str(self.__raw_qry[QRY_ALT_COL_IDX])

    @property
    def ref_mutated(self):
        return self.__get_anno(REF_MUTATED_COL_NAME) == 1

#    @property
#    def FILTER(self):
#        return self.__raw_qry[QRY_FILTER_COL_IDX]
#
    def __parse_genotype(self):
        raw_qry = self.__raw_qry
        self.gt = {}
        for sample_id, sample_col_idx in self.__sample_col_idxs.items():
            raw_zygo = raw_qry[sample_col_idx]
            if raw_zygo != GT_HOM and raw_zygo != GT_WT:
                actual_zygo = raw_zygo
            # if reference is not mutated the actual gt remain the same
            elif not self.ref_mutated:
                actual_zygo = raw_zygo
            # below should be only hom and wild type with allele freq >= 0.5
            elif raw_zygo == GT_HOM:
                actual_zygo = GT_WT
            else:
                actual_zygo = GT_HOM
            if actual_zygo == GT_HOM or actual_zygo == GT_HET:
                mutated = True
            else:
                mutated = False
            gt = QryCall()
            gt.raw_zygo = raw_zygo
            gt.actual_zygo = actual_zygo
            gt.mutated = mutated
            self.gt[sample_id] = gt

    def __cal_shared(self):
        """
        - with given family informations, this function will identify shared
        mutation for each family and for each genotype
        """
        if self.__family_infos is None:
            return
        for fam_id in self.__family_infos:
            family_info = self.__family_infos[fam_id]
            shared_mutation = True
            for member in family_info.members:
                if not self.gt[member.sample_id].mutated:
                    shared_mutation = False
                    break
            for member in family_info.members:
                self.gt[member.sample_id].shared_mutation = shared_mutation
            family_info.shared_mutation = shared_mutation

    def __parse_exac_constraint(self, exac_constraint_vals):
        exac_constaints = exac_constraint_vals.split('#')
        self.__exac_constraints = {}
        for exac_constaint in exac_constaints:
            match = EXAC_CONSTRAINT_PATTERN.match(exac_constaint)
            var_name = match.group('var_name')
            value = match.group('value')
            self.__exac_constraints[var_name] = value

    def __get_exac_constraint_val(self, var_name):
        if self.__exac_constraints is None:
            return None
        if self.__exac_constraints == -1:
            exac_constraint_vals = self.get_anno(EXAC03_CONSTRAINT_COL_NAME)
            if (exac_constraint_vals is None or
                exac_constraint_vals == "." or
                exac_constraint_vals == ""):
                self.__exac_constraints = None
                return None
            self.__parse_exac_constraint(exac_constraint_vals)
        if var_name in self.__exac_constraints:
            return self.__exac_constraints[var_name]
        raise IndexError("constraint: "+var_name+" is not found in exac03 constraint")

    def __get_max_ref_maf(self):
        if self.__max_ref_maf is None:
            max_ref_maf = 0
            for ref_maf_col_name in REF_MAF_COL_NAMES:
                ref_maf = self.get_anno(ref_maf_col_name)
                if ref_maf == "":
                    continue
                ref_maf = float(ref_maf)
                if ref_maf > 0.5:
                    ref_maf = 1 - ref_maf
                if ref_maf > max_ref_maf:
                    max_ref_maf = ref_maf
            self.__max_ref_maf = max_ref_maf
        return self.__max_ref_maf

    def __get_pathogenic_count(self):
        if self.__pathogenic_count is None:
            count = 0
            for col_name in PREDICTION_COLS:
                info = self.get_anno(col_name)
                if (info is None or info == ""):
                    count += 1
                    continue
                if not info.harmful:
                    continue
                count += 1
            self.__pathogenic_count = count
        return self.__pathogenic_count

    def __get_anno(self, anno_col):
        # - If info exist (maybe in an "encrypted" name),
        # return the info
        if anno_col in self.__anno_col_idxs:
            return self.__raw_qry[self.__anno_col_idxs[anno_col]]
        return None

    def get_anno(self, anno_col):
        # - If info exist (maybe in an "encrypted" name),
        # return the info
        # - If info calculate-able return the calculated result
        # - Otherwise, exception
        # - Before return simplify all possible NA into ""
        anno_val = self.__get_anno(anno_col)
        # Translate in-silico prediction if possible
        if (anno_col in self.__pred_tran.predictor_list and
            anno_val is not None):
            anno_val = self.__pred_tran.get_prediction_info(anno_col, anno_val)
        if anno_val is not None:
            pass
        elif anno_col in EXAC03_CONSTRAINT_COL_NAMES:
            anno_val = self.__get_exac_constraint_val(anno_col)
        elif anno_col == INTERVAR_CLASS_COL_NAME:
            anno_val = parse_intervar_class(self.get_anno(INTERVAR_AND_EVIDENCE_COL_NAME))
        elif anno_col == INTERVAR_EVIDENCE_COL_NAME:
            anno_val = parse_intervar_evidence(self.get_anno(INTERVAR_AND_EVIDENCE_COL_NAME))
        elif anno_col == MAX_REF_MAF_COL_NAME:
            anno_val = self.__get_max_ref_maf()
        elif anno_col == PATHOGENIC_COUNT_COL_NAME:
            anno_val = self.__get_pathogenic_count()
        if (anno_val == "" or
            anno_val is None or
            anno_val == [None] or
            anno_val == "."
            ):
            anno_val = ""
        return anno_val

    def has_mutation(self, samples_id):
        """
        to identify if a genotype mutation is found in any samples given
           - samples, list of samples to be compared if the list is
             empty, None, it will return False
        """
        if samples_id is None:
            return False
        if len(samples_id) == 0:
            return False
        for sample_id in samples_id:
            if self.gt[sample_id].mutated:
                return True
        return False

    def has_shared(self):
        for fam_id in self.__family_infos:
            family_info = self.__family_infos[fam_id]
            if family_info.shared_mutation:
                return True
        return False

    def qry_eval(self, expr):
        def anno_repl(match_obj):
            repl_txt = "self.get_anno("
            repl_txt += match_obj.group(0)
            repl_txt += ")"
            return repl_txt
        def is_str(val):
            return not eval("is_number('" + val + "')")

        # look for annotated fields
        anno_fields = {}
        for match in re.finditer(r'(\".+?\")', expr):
            anno_fields[match.group(1)] = 1
        # replace the annotated fields with the actual values
        repl_expr = expr
        for anno_field in anno_fields:
            anno_val = str(eval(re.sub(r'(\".+?\")',
                                       anno_repl,
                                       anno_field)))
            if is_str(anno_val):
                anno_val = "'" + anno_val + "'"
            repl_expr = re.sub(anno_field,
                               anno_val,
                               repl_expr)
        # then eval the expression
        return eval(repl_expr)

class SQLiteDBReader(SQLiteDB):
    """ A class to provide APIs for reading data CMM custrom SQLite DB """

    def __init__(self, db_path, family_infos=None, *args, **kwargs):
        kwargs['db_path'] = db_path
        super(SQLiteDBReader, self).__init__(*args, **kwargs)
        self.__qry_cols = None
        self.__family_infos = family_infos
        self.__pred_tran = PredictionTranslator()
        self.__init_properties()

    def __init_properties(self):
        self.__filter_non_intergenic = False
        self.__filter_non_intronic = False
        self.__filter_non_upstream = False
        self.__filter_non_downstream = False
        self.__filter_non_utr = False
        self.__filter_non_synonymous = False

    def get_raw_obj_str(self, *args, **kwargs):
        raw_str = super(SQLiteDBReader, self).get_raw_obj_str(*args, **kwargs)
        return raw_str

    def init_columns(self, anno_col_names):
        # - To tell the reader what columns it should expect
        # - To prepare the column names for later sql query
        # - Because the UI col name and sql col name might be different,
        # the function also maps UI col name and col name idx in the query
        # - The return value is a list of available columns
        self.__anno_col_idxs = OrderedDict()
        self.__qry_cols = []
        avail_cols = []
        anno_col_idx = len(VCF_PKEYS)
        for anno_col_name in anno_col_names:
            db_col_name = self.anno_col_to_db_col(anno_col_name)
            if db_col_name is not None:
                avail_cols.append(anno_col_name)
                self.__anno_col_idxs[anno_col_name] = anno_col_idx
                self.__qry_cols.append(db_col_name)
                anno_col_idx += 1
            elif (anno_col_name in EXAC03_CONSTRAINT_COL_NAMES and
                self.anno_col_to_db_col(EXAC03_CONSTRAINT_COL_NAME) is not None
                ):
                avail_cols.append(anno_col_name)
                if EXAC03_CONSTRAINT_COL_NAME not in self.__anno_col_idxs:
                    self.__anno_col_idxs[EXAC03_CONSTRAINT_COL_NAME] = anno_col_idx
                    db_col_name = self.anno_col_to_db_col(EXAC03_CONSTRAINT_COL_NAME)
                    self.__qry_cols.append(db_col_name)
                    anno_col_idx += 1
            elif ((anno_col_name == INTERVAR_CLASS_COL_NAME or
                anno_col_name == INTERVAR_EVIDENCE_COL_NAME) and
                self.anno_col_to_db_col(INTERVAR_AND_EVIDENCE_COL_NAME) is not None
                ):
                avail_cols.append(anno_col_name)
                if INTERVAR_AND_EVIDENCE_COL_NAME not in self.__anno_col_idxs:
                    self.__anno_col_idxs[INTERVAR_AND_EVIDENCE_COL_NAME] = anno_col_idx
                    db_col_name = self.anno_col_to_db_col(INTERVAR_AND_EVIDENCE_COL_NAME)
                    self.__qry_cols.append(db_col_name)
                    anno_col_idx += 1
            elif anno_col_name == MAX_REF_MAF_COL_NAME:
                avail_cols.append(anno_col_name)
            elif anno_col_name == PATHOGENIC_COUNT_COL_NAME:
                avail_cols.append(anno_col_name)
            else:
                self.warning("Column " + anno_col_name + " is missing")
        return avail_cols

    def init_samples(self, samples_id=None):
        # - To tell the reader what samples'id it should expect
        # - To prepare the column names for later sql query
        # - Because the UI sample id and sql sample id might be different,
        # the function also maps UI sample id and sample id idx in the query
        self.__sample_col_idxs = OrderedDict()
        self.__sample_tbls = OrderedDict()
        sample_col_idx = len(VCF_PKEYS) + len(self.__anno_col_idxs)
        if samples_id is None:
            samples_id = self.samples_id
        for sample_id in samples_id:
            db_col_name = self.sample_id_to_db_col(sample_id)
            if db_col_name is not None:
                self.__sample_col_idxs[sample_id] = sample_col_idx
                self.__qry_cols.append(db_col_name)
                sample_col_idx += 1
            else:
                raise IndexError("sample id: "+sample_id+" is not found in the database")
            tbl_name = self.sample_id_to_tbl(sample_id)
            if tbl_name is not None:
                self.__sample_tbls[tbl_name] = 1
            else:
                raise IndexError("sample id: "+sample_id+" is not found in any tables")
        self.__ref_mutated_col_idx = sample_col_idx
            
    @property
    def filter_non_intergenic(self):
        return self.__filter_non_intergenic

    @filter_non_intergenic.setter
    def filter_non_intergenic(self, value):
        self.__filter_non_intergenic = value

    @property
    def filter_non_intronic(self):
        return self.__filter_non_intronic

    @filter_non_intronic.setter
    def filter_non_intronic(self, value):
        self.__filter_non_intronic = value

    @property
    def filter_non_upstream(self):
        return self.__filter_non_upstream

    @filter_non_upstream.setter
    def filter_non_upstream(self, value):
        self.__filter_non_upstream = value

    @property
    def filter_non_downstream(self):
        return self.__filter_non_downstream

    @filter_non_downstream.setter
    def filter_non_downstream(self, value):
        self.__filter_non_downstream = value

    @property
    def filter_non_utr(self):
        return self.__filter_non_utr

    @filter_non_utr.setter
    def filter_non_utr(self, value):
        self.__filter_non_utr = value

    @property
    def filter_non_synonymous(self):
        return self.__filter_non_synonymous

    @filter_non_synonymous.setter
    def filter_non_synonymous(self, value):
        self.__filter_non_synonymous = value

    @property
    def filter_pass_vqsr(self):
        pass

    @filter_pass_vqsr.setter
    def filter_pass_vqsr(self, value):
        pass

    def get_qry_records(self, chrom=None, start_pos=None, end_pos=None):
        sql = "SELECT "
        if self.__qry_cols is None:
            anno_cols = self.get_avdb_info().values()
            anno_cols += self.get_annovar_info().values()
            anno_cols = map(lambda x: x.strip("_"),
                            reduce(lambda x, y: x+y,
                                   anno_cols))
            anno_cols.append(MAX_REF_MAF_COL_NAME)
            anno_cols.append(REF_MUTATED_COL_NAME)
            anno_cols.append(INTERVAR_AND_EVIDENCE_COL_NAME)
            self.init_columns(anno_cols)
            self.init_samples(self.samples_id)
        sql += ", ".join(VCF_PKEYS)
        sql += ", " + ", ".join(self.__qry_cols)
        sql += " FROM " + MASTER_TABLE
        where_clause = []
        if chrom is not None:
            where_clause.append("CHROM = '" + chrom + "'")
            if start_pos is not None:
                where_clause.append("POS >= " + start_pos)
                where_clause.append("POS <= " + end_pos)
        if self.filter_non_intergenic:
            where_clause.append("_Func_refGene != '" + FUNC_INTERGENIC + "'")
        if self.filter_non_intronic:
            where_clause.append("_Func_refGene != '" + FUNC_INTRONIC + "'")
        if self.filter_non_upstream:
            where_clause.append("_Func_refGene != '" + FUNC_UPSTREAM + "'")
        if self.filter_non_downstream:
            where_clause.append("_Func_refGene != '" + FUNC_DOWNSTREAM + "'")
        if self.filter_non_utr:
            where_clause.append("_Func_refGene != '" + FUNC_UTR + "3'")
            where_clause.append("_Func_refGene != '" + FUNC_UTR + "5'")
        if self.filter_non_synonymous:
            where_clause.append("_ExonicFunc_refGene != '" + EXONICFUNC_SYNONYMOUS + "'")
        or_clause = []
        for tbl_name in self.__sample_tbls:
            col_name = tbl_name.replace("gtz_", "") + "_FILTER"
            or_clause.append(col_name + " IS NOT NULL")
        where_clause.append("(" + " OR ".join(or_clause) + ")")
        if len(where_clause) > 0:
            sql += " WHERE " + " AND ".join(where_clause)
        for row in self.read_rows(sql=sql):
            yield QryRecord(row,
                            self.__anno_col_idxs,
                            self.__sample_col_idxs,
                            self.__family_infos,
                            self.__pred_tran,
                            )
