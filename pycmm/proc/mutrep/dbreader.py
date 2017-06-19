from pycmm.proc.db.dbms import SQLiteDB


class SQLiteDBReader(SQLiteDB):
    """ A class to provide APIs for reading data CMM custrom SQLite DB """

    def __init__(self, *args, **kwargs):
        super(SQLiteDBReader, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(SQLiteDB, self).get_raw_obj_str(*args, **kwargs)
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
