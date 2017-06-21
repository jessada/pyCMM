import sqlite3
from pycmm.template import pyCMMBase

DATA_TYPE_AVDB_INFO = "data_type_avdb_info"
DATA_TYPE_ANNOVAR_INFO = "data_type_annovar_info"
DATA_TYPE_GTZ = "data_type_gtz"

# To keep the information about table created for CMM usage
# So far only keep data type and information columns
TBL_NAME_CMM_TBLS = "cmm_tables"
TBL_CMM_TBLS_INFO_COLS_COL_NAME = "info_cols"
TBL_CMM_TBLS_TBL_NAME_COL_NAME = "tbl_name"
TBL_CMM_TBLS_DATA_TYPE_COL_NAME = "data_type"


class SQLiteDB(pyCMMBase):
    """
    A class to provide APIs for conntecting to CMM custrom SQLite DB.
    The implementation include the handling of CMM-related-system table
    """

    def __init__(self, db_path, verbose=True, *args, **kwargs):
        super(SQLiteDB, self).__init__(*args, **kwargs)
        self.__db_path = db_path
        self.__verbose = verbose
        self.__init_db()

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(SQLiteDB, self).get_raw_obj_str(*args, **kwargs)
        return raw_repr

    @property
    def db_path(self):
        return self.__db_path

    @property
    def verbose(self):
        return self.__verbose

    def __init_db(self):
        self.info()
        self.info("connecting to database: " + self.db_path)
        self._create_table(tbl_name=TBL_NAME_CMM_TBLS,
                            col_names=[TBL_CMM_TBLS_TBL_NAME_COL_NAME,
                                       TBL_CMM_TBLS_INFO_COLS_COL_NAME,
                                       TBL_CMM_TBLS_DATA_TYPE_COL_NAME],
                            pkeys=[TBL_CMM_TBLS_TBL_NAME_COL_NAME,
                                   TBL_CMM_TBLS_INFO_COLS_COL_NAME],
                            )

    def __log_sql(self, sql):
        if self.verbose:
            self.info("executing sql: " + sql)

    def _fetch_all(self, sql):
        # to encapsulate sqlite sql fetchall execution
        # so far just for logging
        conn = self.__connect_db()
        self.__log_sql(sql)
        result = conn.execute(sql).fetchall()
        self.__close_connection()
        return result

    def _fetch_none(self, sql):
        # to encapsulate sqlite sql fetchnone execution
        # so far just for logging
        conn = self.__connect_db()
        self.__log_sql(sql)
        result = conn.execute(sql).fetchone()
        self.__close_connection()
        return result

    def _exec_sql(self, sql):
        # to encapsulate sqlite sql execution
        # so far just for logging
        conn = self.__connect_db() 
        self.__log_sql(sql)
        result = conn.execute(sql)
        self.__close_connection()
        return result

    def _exec_many_sql(self, sql, records):
        # to encapsulate sqlite sql execution
        # so far just for logging
        conn = self.__connect_db() 
        self.__log_sql(sql)
        result = conn.executemany(sql, records)
        self.__close_connection()
        return result

    def __connect_db(self):
        self.__conn = sqlite3.connect(self.db_path)
        return self.__conn

    def __close_connection(self):
        self.__conn.commit()
        self.__conn.close()

    def drop_table(self, tbl_name):
        sql = 'DROP TABLE IF EXISTS ' + tbl_name
        self._exec_sql(sql)

    def drop_view(self, view_name):
        sql = 'DROP VIEW IF EXISTS ' + view_name
        self._exec_sql(sql)

    def _create_table(self, tbl_name, col_names, pkeys=None):
        sql = 'CREATE TABLE IF NOT EXISTS ' + tbl_name
        sql += ' ( ' + ", ".join(col_names) 
        if pkeys is not None:
            sql += ', PRIMARY KEY (' + ",".join(pkeys) + ')' 
        sql += ' )'
        self._exec_sql(sql)

    def _update_sys_info(self,
                         data_type,
                         updating_info_cols,
                         updating_tbl_name,
                         ):
        # For adding system information
        # This has be called when AVDB, Annovar Info and genotyping data
        # are loaded into the database
        # The process start by checking if the "updating_tbl_name
        # is alreay there
        # if yes, insert else update
        where_clause = TBL_CMM_TBLS_TBL_NAME_COL_NAME + "='" + updating_tbl_name + "'"
        col_names_txt = ",".join(updating_info_cols)
        if self.count_rows(TBL_NAME_CMM_TBLS, where_clause) > 0:
            sql = 'UPDATE ' + TBL_NAME_CMM_TBLS
            sql += ' SET ' + TBL_CMM_TBLS_INFO_COLS_COL_NAME + "='" + col_names_txt + "'"
            sql += ' WHERE ' + TBL_CMM_TBLS_TBL_NAME_COL_NAME + "='" + updating_tbl_name + "'"
        else:
            sql = 'INSERT INTO ' + TBL_NAME_CMM_TBLS
            sql += " ( " + TBL_CMM_TBLS_TBL_NAME_COL_NAME
            sql += ", " + TBL_CMM_TBLS_INFO_COLS_COL_NAME
            sql += ", " + TBL_CMM_TBLS_DATA_TYPE_COL_NAME + " ) "
            sql += " select '" + updating_tbl_name + "'"
            sql += ", '" + col_names_txt + "'"
            sql += ", '" + data_type + "'"
        self.dbg(sql)
        self._exec_sql(sql)

    def table_exist(self, tbl_name):
        sql = "SELECT name FROM sqlite_master WHERE type='table'"
        sql += "AND name='" + tbl_name + "'"
        tbl_exist = self._fetch_none(sql)
        if tbl_exist:
            return True
        return False

    def count_rows(self, tbl_name, where_clause=None):
        sql = "SELECT COUNT(*) FROM " + tbl_name
        if where_clause is not None:
            sql += " WHERE " + where_clause
        (rows_count,) = self._fetch_none(sql)
        return rows_count

    def __get_sys_info(self, datatype):
        sys_info = {}
        sql = "SELECT " + TBL_CMM_TBLS_TBL_NAME_COL_NAME
        sql += ", " + TBL_CMM_TBLS_INFO_COLS_COL_NAME
        sql += " FROM " + TBL_NAME_CMM_TBLS
        sql += " WHERE " + TBL_CMM_TBLS_DATA_TYPE_COL_NAME
        sql += " = '" + datatype + "'"
        rows = self._fetch_all(sql)
        for (tbl_name, info_cols_txt) in rows:
            sys_info[tbl_name] = info_cols_txt.split(",")
        return sys_info

    def get_avdb_info(self):
        return self.__get_sys_info(DATA_TYPE_AVDB_INFO)

    def get_annovar_info(self):
        return self.__get_sys_info(DATA_TYPE_ANNOVAR_INFO)

    def get_samples_id(self, gtz_tbl_name=None):
        samples_info = self.__get_sys_info(DATA_TYPE_GTZ)
        if gtz_tbl_name is not None:
            return samples_info[gtz_tbl_name]
        return reduce(lambda x, y: x+y,
                      samples_info.values())

    def get_col_names(self, tbl_name):
        cursor = self._exec_sql('SELECT * FROM ' + tbl_name + " LIMIT 1")
        col_names = map(lambda x: x[0], cursor.description)
        return col_names

    def read_rows(self, tbl_name):
        conn = self.__connect_db() 
        sql = "SELECT * FROM " + tbl_name
        rows = conn.execute(sql)
        for row in rows:
            yield(row)
        self.__close_connection()
