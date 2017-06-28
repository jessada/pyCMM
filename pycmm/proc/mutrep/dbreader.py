from collections import OrderedDict
from pycmm.utils import is_number
from pycmm.template import pyCMMBase
from pycmm.proc.db.connector import SQLiteDB
from pycmm.proc.db.connector import VCF_PKEYS
from pycmm.proc.db.connector import TBL_NAME_ALL_GTZ_ANNOS

MASTER_TABLE = TBL_NAME_ALL_GTZ_ANNOS

QRY_CHROM_COL_IDX = 0
QRY_POS_COL_IDX = 1
QRY_REF_COL_IDX = 2
QRY_ALT_COL_IDX = 3
#QRY_FILTER_COL_IDX = 4


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
                 *args,
                 **kwargs
                 ):
        super(QryRecord, self).__init__(*args, **kwargs)
        self.__raw_qry = raw_qry
        self.__anno_col_idxs = anno_col_idxs
        self.__sample_col_idxs = sample_col_idxs

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

#    @property
#    def FILTER(self):
#        return self.__raw_qry[QRY_FILTER_COL_IDX]
#
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
        return self.__get_anno(anno_col)

    def __get_gt(self, sample_id):
        # - If sample_id exist (maybe in an "encrypted" name),
        # return the gt data
        if sample_id in self.__sample_col_idxs:
            return self.__raw_qry[self.__sample_col_idxs[sample_id]]
        return None

    def get_call(self, sample_id):
        return self.__get_gt(sample_id)

# ********************************** Need testing ********************************
    def qry_eval(self, expr):
        def anno_repl(match_obj):
            repl_txt = "self.get_anno("
            repl_txt += match_obj.group(0)
            repl_txt += ")"
            return repl_txt

        # look for annotated fields
        anno_fields = {}
        for match in re.finditer(r'(\".+?\")', expr):
            anno_fields[match.group(1)] = 1
        # replace the annotated fields with the actual values
        repl_expr = expr
        for anno_field in ano_fields:
            anno_val = eval(re.sub(r'(\".+?\")',
                                   anno_repl,
                                   anno_field)) 
            repl_expr = re.sub(anno_field,
                               "'"+str(anno_val)+"'",
                               repl_expr)
        # then eval the expression
        return eval(repl_expr)
# ********************************** Need testing ********************************

class SQLiteDBReader(SQLiteDB):
    """ A class to provide APIs for reading data CMM custrom SQLite DB """

    def __init__(self, *args, **kwargs):
        super(SQLiteDBReader, self).__init__(*args, **kwargs)
        self.__samples_id = self.__parse_samples()
#        self.__connected_table = None

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(SQLiteDBReader, self).get_raw_obj_str(*args, **kwargs)
        return raw_repr


# obsolete
#    def create_mutrep_view(self, view_name, gtz_tbl_name):
#        # create view for further used in mutation report
#        # it's based on an assumption that all avdb columns in the database
#        # will be included into the view
#        samples_id = self.get_samples_id(gtz_tbl_name=gtz_tbl_name)
#        avdb_info = self.get_avdb_info()
#        annovar_info = self.get_annovar_info()
#        self.drop_view(view_name)
#        sql = 'CREATE VIEW IF NOT EXISTS ' + view_name + " AS "
#        sql += " SELECT gtz.CHROM, gtz.POS, gtz.REF, gtz.ALT, gtz.FILTER"
#        for annovar_tbl_name in annovar_info:
#            col_names = annovar_info[annovar_tbl_name]
#            for col_name in col_names:
#                sql += ", " + annovar_tbl_name + "." + col_name
#        for avdb_tbl_name in avdb_info:
#            col_names = avdb_info[avdb_tbl_name]
#            for col_name in col_names:
#                sql += ", " + avdb_tbl_name + "." + col_name
#        sql += ", " + ", ".join(map(lambda x: "gtz."+x, samples_id))
#        sql += " FROM " + gtz_tbl_name + " AS gtz"
#        for annovar_tbl_name in annovar_info:
#            sql += " LEFT OUTER JOIN " + annovar_tbl_name + " ON "
#            join_conditions = []
#            join_conditions.append("gtz.CHROM=" + annovar_tbl_name + ".CHROM")
#            join_conditions.append("gtz.POS=" + annovar_tbl_name + ".POS")
#            join_conditions.append("gtz.REF=" + annovar_tbl_name + ".REF")
#            join_conditions.append("gtz.ALT=" + annovar_tbl_name + ".ALT")
#            sql += " AND ".join(join_conditions)
#        for avdb_tbl_name in avdb_info:
#            sql += " LEFT OUTER JOIN " + avdb_tbl_name + " ON "
#            join_conditions = []
#            join_conditions.append("gtz.CHROM=" + avdb_tbl_name + ".Chr")
#            join_conditions.append("gtz.POS=" + avdb_tbl_name + ".Start")
#            join_conditions.append("gtz.REF=" + avdb_tbl_name + ".Ref")
#            join_conditions.append("gtz.ALT=" + avdb_tbl_name + ".ALt")
#            sql += " AND ".join(join_conditions)
#        self._exec_sql(sql)

# ********************************** test is required for the whole thing inside here ********************************
    def __parse_samples(self):
        master_col_names = self.get_col_names(MASTER_TABLE)
        db_samples_id = self.get_samples_id()
        samples_id = []
        for col_name in master_col_names:
            if col_name in db_samples_id:
                samples_id.append(col_name.strip("_").replace("_", "-"))
        return samples_id

    def init_columns(self, col_names):
        # - To tell the reader what columns it should expect
        # - To prepare the column names for later sql query
        # - Because the UI col name and sql col name might be different,
        # the function also maps UI col name and col name idx in the query
        # - The return value is a list of available columns
        master_col_names = self.get_col_names(MASTER_TABLE)
        mod_master_col_names = map(lambda x: x.strip("_"),
                                   master_col_names)
        self.__anno_col_idxs = OrderedDict()
        self.__qry_cols = []
        anno_col_idx = len(VCF_PKEYS)
        for col_name in col_names:
            if col_name in mod_master_col_names:
                self.__anno_col_idxs[col_name] = anno_col_idx
                master_col_idx = mod_master_col_names.index(col_name)
                self.__qry_cols.append(master_col_names[master_col_idx])
                anno_col_idx += 1
            else:
                self.warning("Column " + col_name + " is missing")
        avail_cols = self.__anno_col_idxs.keys()
        return avail_cols

    def init_samples(self, samples_id=None):
        # - To tell the reader what samples'id it should expect
        # - To prepare the column names for later sql query
        # - Because the UI sample id and sql sample id might be different,
        # the function also maps UI sample id and sample id idx in the query
        master_col_names = self.get_col_names(MASTER_TABLE)
        mod_master_col_names = map(lambda x: x.strip("_").replace("_", "-"),
                                   master_col_names)
        self.__sample_col_idxs = OrderedDict()
        sample_col_idx = len(VCF_PKEYS) + len(self.__anno_col_idxs)
        if samples_id is None:
            samples_id = self.samples_id
        for sample_id in samples_id:
            if sample_id in mod_master_col_names:
                self.__sample_col_idxs[sample_id] = sample_col_idx
                master_col_idx = mod_master_col_names.index(sample_id)
                self.__qry_cols.append(master_col_names[master_col_idx])
                sample_col_idx += 1
            else:
                raise IndexError("sample id: "+sample_id+" is not found in the database")
            
    @property
    def samples_id(self):
        return self.__samples_id

    @property
    def filter_pass_vqsr(self):
        pass

    @filter_pass_vqsr.setter
    def filter_pass_vqsr(self, value):
        pass

    def get_qry_records(self, chrom=None, start_pos=None, end_pos=None):
        sql = "SELECT "
        sql += ", ".join(VCF_PKEYS)
        sql += ", " + ", ".join(self.__qry_cols)
        sql += " FROM " + MASTER_TABLE
        if chrom is not None:
            sql += " WHERE CHROM = '" + chrom + "'"
            if start_pos is not None:
                sql += " AND POS >= " + start_pos 
                sql += " AND POS <= " + end_pos 
        for row in self.read_rows(sql=sql):
            yield QryRecord(row, self.__anno_col_idxs, self.__sample_col_idxs)

#    def __cal_qry_col_idxs(self, tbl_name, qry_col_names):
#        tbl_col_names = map(lambda x: x.strip("_"),
#                            self.get_col_names(tbl_name))
##        tbl_col_names = self.get_col_names(tbl_name)
#        qry_col_idxs = OrderedDict()
#        for qry_col_name in qry_col_names:
#            if qry_col_name in tbl_col_names:
#                qry_col_idxs[qry_col_name] = tbl_col_names.index(qry_col_name)
#            else:
#                self.warning("Column " + qry_col_name + " is missing")
##                qry_col_idxs[qry_col_name] = -1
#        return qry_col_idxs
#
#    def __cal_qry_sample_idxs(self, tbl_name):
#        # This function will only check if any columns in the tableis containing
#        # genotyping data as the table itself doesn't know anyting.
#        # The similar comparison with the samples id from job setup file is not 
#        # done here with the assumption that the samples id must have been used in
#        # generating table/view upstream so all the samples id are needed
#        tbl_col_names = self.get_col_names(tbl_name)
#        db_samples_id = self.get_samples_id()
#        qry_sample_idxs = OrderedDict()
#        for tbl_col_idx in xrange(len(tbl_col_names)):
#            tbl_col_name = tbl_col_names[tbl_col_idx]
#            if tbl_col_name in db_samples_id:
#                sample_id = tbl_col_name.strip("_").replace("_", "-")
#                qry_sample_idxs[sample_id] = tbl_col_idx
#        return qry_sample_idxs
#
#    def query_table(self, tbl_name, col_names):
#        self.__connected_table = tbl_name
#        self.__qry_col_idxs = self.__cal_qry_col_idxs(tbl_name, col_names)
#        self.__qry_sample_idxs = self.__cal_qry_sample_idxs(tbl_name)
#
#    @property
#    def connected_table(self):
#        return self.__connected_table
#
#    @property
#    def qry_col_idxs(self):
#        return self.__qry_col_idxs
#
#    @property
#    def qry_sample_idxs(self):
#        return self.__qry_sample_idxs
#
#    @property
#    def qry_records(self):
#        for row in self.read_rows(self.connected_table):
#            yield row
# ********************************** test is required for the whole thing inside here ********************************
