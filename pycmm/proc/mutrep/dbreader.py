from collections import OrderedDict
from pycmm.proc.db.dbms import SQLiteDB


class SQLiteDBReader(SQLiteDB):
    """ A class to provide APIs for reading data CMM custrom SQLite DB """

    def __init__(self, *args, **kwargs):
        super(SQLiteDBReader, self).__init__(*args, **kwargs)
        self.__connected_table = None

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(SQLiteDBReader, self).get_raw_obj_str(*args, **kwargs)
        return raw_repr

    def create_mutrep_view(self, view_name, gtz_tbl_name):
        # create view for further used in mutation report
        # it's based on an assumption that all avdb columns in the database
        # will be included into the view
        samples_id = self.get_samples_id(gtz_tbl_name=gtz_tbl_name)
        avdb_info = self.get_avdb_info()
        annovar_info = self.get_annovar_info()
        self.drop_view(view_name)
        sql = 'CREATE VIEW IF NOT EXISTS ' + view_name + " AS "
        sql += " SELECT gtz.CHROM, gtz.POS, gtz.REF, gtz.ALT, gtz.FILTER"
        for annovar_tbl_name in annovar_info:
            col_names = annovar_info[annovar_tbl_name]
            for col_name in col_names:
                sql += ", " + annovar_tbl_name + "." + col_name
        for avdb_tbl_name in avdb_info:
            col_names = avdb_info[avdb_tbl_name]
            for col_name in col_names:
                sql += ", " + avdb_tbl_name + "." + col_name
        sql += ", " + ", ".join(map(lambda x: "gtz."+x, samples_id))
        sql += " FROM " + gtz_tbl_name + " AS gtz"
        for annovar_tbl_name in annovar_info:
            sql += " LEFT OUTER JOIN " + annovar_tbl_name + " ON "
            join_conditions = []
            join_conditions.append("gtz.CHROM=" + annovar_tbl_name + ".CHROM")
            join_conditions.append("gtz.POS=" + annovar_tbl_name + ".POS")
            join_conditions.append("gtz.REF=" + annovar_tbl_name + ".REF")
            join_conditions.append("gtz.ALT=" + annovar_tbl_name + ".ALT")
            sql += " AND ".join(join_conditions)
        for avdb_tbl_name in avdb_info:
            sql += " LEFT OUTER JOIN " + avdb_tbl_name + " ON "
            join_conditions = []
            join_conditions.append("gtz.CHROM=" + avdb_tbl_name + ".Chr")
            join_conditions.append("gtz.POS=" + avdb_tbl_name + ".Start")
            join_conditions.append("gtz.REF=" + avdb_tbl_name + ".Ref")
            join_conditions.append("gtz.ALT=" + avdb_tbl_name + ".ALt")
            sql += " AND ".join(join_conditions)
        self._exec_sql(sql)

# ********************************** test is required for the whole thing inside here ********************************
    def __cal_qry_col_idxs(self, tbl_name, qry_col_names):
        tbl_col_names = map(lambda x: x.strip("_"),
                            self.get_col_names(tbl_name))
#        tbl_col_names = self.get_col_names(tbl_name)
        qry_col_idxs = OrderedDict()
        for qry_col_name in qry_col_names:
            if qry_col_name in tbl_col_names:
                qry_col_idxs[qry_col_name] = tbl_col_names.index(qry_col_name)
            else:
                self.warning("Column " + qry_col_name + " is missing")
#                qry_col_idxs[qry_col_name] = -1
        return qry_col_idxs

    def __cal_qry_sample_idxs(self, tbl_name):
        # This function will only check if any columns in the tableis containing
        # genotyping data as the table itself doesn't know anyting.
        # The similar comparison with the samples id from job setup file is not 
        # done here with the assumption that the samples id must have been used in
        # generating table/view upstream so all the samples id are needed
        tbl_col_names = self.get_col_names(tbl_name)
        db_samples_id = self.get_samples_id()
        qry_sample_idxs = OrderedDict()
        for tbl_col_idx in xrange(len(tbl_col_names)):
            tbl_col_name = tbl_col_names[tbl_col_idx]
            if tbl_col_name in db_samples_id:
                sample_id = tbl_col_name.strip("_").replace("_", "-")
                qry_sample_idxs[sample_id] = tbl_col_idx
        return qry_sample_idxs

    def query_table(self, tbl_name, col_names):
        self.__connected_table = tbl_name
        self.__qry_col_idxs = self.__cal_qry_col_idxs(tbl_name, col_names)
        self.__qry_sample_idxs = self.__cal_qry_sample_idxs(tbl_name)

    @property
    def connected_table(self):
        return self.__connected_table

    @property
    def qry_col_idxs(self):
        return self.__qry_col_idxs

    @property
    def qry_sample_idxs(self):
        return self.__qry_sample_idxs

    @property
    def qry_records(self):
        for row in self.read_rows(self.connected_table):
            yield row
# ********************************** test is required for the whole thing inside here ********************************
