import sqlite3
from pycmm.settings import MAX_REF_MAF_COL_NAME
from pycmm.settings import REF_MAF_COL_NAMES
#from pycmm.settings import INTERVAR_AND_EVIDENCE_COL_NAME
from pycmm.settings import GNOMAD_GENOME_ALL_COL_NAME
from pycmm.settings import ALL_MUTREP_ANNO_COLS
from pycmm.template import pyCMMBase

DATA_TYPE_AVDB_INFO = "data_type_avdb_info"
DATA_TYPE_ANNOVAR_INFO = "data_type_annovar_info"
DATA_TYPE_GTZ = "data_type_gtz"

VCF_PKEYS = []
VCF_PKEYS.append("CHROM")
VCF_PKEYS.append("POS")
VCF_PKEYS.append("REF")
VCF_PKEYS.append("ALT")

AVDB_PKEYS = []
AVDB_PKEYS.append("Chr")
AVDB_PKEYS.append("Start")
AVDB_PKEYS.append("Ref")
AVDB_PKEYS.append("Alt")

# A table to keep track of information columsn that can be used
# for generating reports later on
TBL_NAME_CMM_TBLS = "cmm_tables"
TBL_CMM_TBLS_INFO_COLS_COL_NAME = "info_cols"
TBL_CMM_TBLS_TBL_NAME_COL_NAME = "tbl_name"
TBL_CMM_TBLS_DATA_TYPE_COL_NAME = "data_type"

# A table to keep all the variants coordinate (CHROM, POS, REF, ALT)
# The columns in the table are als CHROM, POS, REF, ALT
TBL_NAME_GTZ_COORS = "gtz_coors"

# A table that have all vcf coordinates of all samples,
# all annotataion columns,
# and genotyping data from all samples
TBL_NAME_ALL_GTZ_ANNOS = "all_gtz_annos"

REF_MUTATED_COL_NAME = "ref_mutated"

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
        raw_str = super(SQLiteDB, self).get_raw_obj_str(*args, **kwargs)
        return raw_str

    @property
    def db_path(self):
        return self.__db_path

    @property
    def verbose(self):
        return self.__verbose

    def __init_db(self):
        self.info()
        self.info("connecting to database: " + self.db_path)
        self.__create_system_tables()
        self.__avdb_info = None
        self.__annovar_info = None
        self.__col_mapping = None
        self.__sample_id_mapping = None
        self.__sample_tbl_mapping = None
        self.__samples_id = None

    def __create_system_tables(self):
        self._create_table(tbl_name=TBL_NAME_CMM_TBLS,
                            col_names=[TBL_CMM_TBLS_TBL_NAME_COL_NAME,
                                       TBL_CMM_TBLS_INFO_COLS_COL_NAME,
                                       TBL_CMM_TBLS_DATA_TYPE_COL_NAME],
                            pkeys=[TBL_CMM_TBLS_TBL_NAME_COL_NAME,
                                   TBL_CMM_TBLS_INFO_COLS_COL_NAME],
                            )
        self._create_table(tbl_name=TBL_NAME_GTZ_COORS,
                            col_names=VCF_PKEYS,
                            pkeys=VCF_PKEYS,
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
        self.__conn.text_factory = str
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

    def vacuum(self):
        sql = 'VACUUM'
        self._exec_sql(sql)

    def _create_table(self, tbl_name, col_names, pkeys=None):
        sql = 'CREATE TABLE IF NOT EXISTS ' + tbl_name
        sql += ' ( ' + ", ".join(col_names) 
        if pkeys is not None:
            sql += ', PRIMARY KEY (' + ",".join(pkeys) + ')' 
        sql += ' )'
        self._exec_sql(sql)

    def _add_column(self, tbl_name, col_name, col_def=""):
        sql = "ALTER TABLE " + tbl_name
        sql += " ADD " + col_name
        sql += " " + col_def
        self._exec_sql(sql)

    def _update_gtz_coors(self,
                          tbl_name,
                          ):
        # insert all the input coors into the table (only the new ones)
        sql = 'INSERT OR IGNORE INTO ' + TBL_NAME_GTZ_COORS
        sql += ' (' + ", ".join(VCF_PKEYS) + ')'
        sql += ' SELECT ' + ", ".join(VCF_PKEYS)
        sql += ' FROM ' + tbl_name
        row_count = self._exec_sql(sql).rowcount
        return row_count

    def _update_tbl_info(self,
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
            sql += " SELECT '" + updating_tbl_name + "'"
            sql += ", '" + col_names_txt + "'"
            sql += ", '" + data_type + "'"
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
        (row_count,) = self._fetch_none(sql)
        return row_count

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
        if self.__avdb_info is None:
            self.__avdb_info = self.__get_sys_info(DATA_TYPE_AVDB_INFO)
        return self.__avdb_info

    def get_annovar_info(self):
        return self.__get_sys_info(DATA_TYPE_ANNOVAR_INFO)

    def get_samples_id(self, gtz_tbl_name=None):
        samples_info = self.__get_sys_info(DATA_TYPE_GTZ)
        if gtz_tbl_name is not None:
            return samples_info[gtz_tbl_name]
        return reduce(lambda x, y: x+y,
                      samples_info.values())

    def get_gtz_tbls(self):
        samples_info = self.__get_sys_info(DATA_TYPE_GTZ)
        return samples_info

    def get_col_names(self, tbl_name):
        cursor = self._exec_sql('SELECT * FROM ' + tbl_name + " LIMIT 1")
        col_names = map(lambda x: x[0], cursor.description)
        return col_names

    def read_rows(self, tbl_name=None, sql=None):
        conn = self.__connect_db() 
        if tbl_name is not None:
            sql = "SELECT * FROM " + tbl_name
        rows = conn.execute(sql)
        for row in rows:
            yield(row)
        self.__close_connection()

    def dbg_query(self, sql):
        for row in self.read_rows(sql=sql):
            self.dbg(row)

    @property
    def samples_id(self):
        if self.__sample_id_mapping is None:
            self.__map_samples()
        if self.__samples_id is None:
            self.__samples_id = self.__sample_id_mapping.keys()
        return self.__samples_id

    def __map_columns(self):
        # Due to the different in annotation column names and the column names
        # in the database, and there are some risks in changing the table
        # columns schema, this funtion will map them in dictionary by using the
        # annotation column name as the key and the table column name as value.
        self.__col_mapping = {}
        db_col_names = self.get_col_names(TBL_NAME_ALL_GTZ_ANNOS)
        mod_db_col_names = map(lambda x: x.strip("_"),
                               db_col_names)
        anno_col_names = self.get_avdb_info().values()
        anno_col_names += self.get_annovar_info().values()
        anno_col_names = map(lambda x: x.strip("_"),
                             reduce(lambda x, y: x+y,
                                    anno_col_names))
        anno_col_names += ALL_MUTREP_ANNO_COLS
        anno_col_names.append(REF_MUTATED_COL_NAME)
#        anno_col_names.append(INTERVAR_AND_EVIDENCE_COL_NAME)
        for anno_col_name in anno_col_names:
            if anno_col_name in mod_db_col_names:
                db_col_idx = mod_db_col_names.index(anno_col_name)
                db_col_name = db_col_names[db_col_idx]
                self.__col_mapping[anno_col_name] = db_col_name
            else:
                mod_anno_col_name = anno_col_name.replace(".", "_")
                if mod_anno_col_name in mod_db_col_names:
                    db_col_idx = mod_db_col_names.index(mod_anno_col_name)
                    db_col_name = db_col_names[db_col_idx]
                    self.__col_mapping[anno_col_name] = db_col_name
                    self.__col_mapping[mod_anno_col_name] = db_col_name
                mod_anno_col_name = anno_col_name.replace(":", "_")
                if mod_anno_col_name in mod_db_col_names:
                    db_col_idx = mod_db_col_names.index(mod_anno_col_name)
                    db_col_name = db_col_names[db_col_idx]
                    self.__col_mapping[anno_col_name] = db_col_name
                    self.__col_mapping[mod_anno_col_name] = db_col_name

    def anno_col_to_db_col(self, anno_col_name):
        if self.__col_mapping is None:
            self.__map_columns()
        if anno_col_name in self.__col_mapping:
            return self.__col_mapping[anno_col_name]
        return None

    def __map_samples(self):
        # Due to the different in the actual sample ids and the column names
        # in the database, this funtion will map them in dictionary by using the
        # actual sample ids as the key and the table column name as value.
        self.__sample_id_mapping = {}
        self.__sample_tbl_mapping = {}
        gtz_tbls = self.get_gtz_tbls()
        db_samples_id = reduce(lambda x, y: x+y,
                               self.get_gtz_tbls().values())
        for gtz_tbl in gtz_tbls:
            for db_sample_id in gtz_tbls[gtz_tbl]:
                actual_id = db_sample_id.strip("_").replace("_", "-")
                self.__sample_id_mapping[actual_id] = db_sample_id
                self.__sample_tbl_mapping[actual_id] = gtz_tbl

    def sample_id_to_db_col(self, sample_id):
        if self.__sample_id_mapping is None:
            self.__map_samples()
        if sample_id in self.__sample_id_mapping:
            return self.__sample_id_mapping[sample_id]
        return None

    def sample_id_to_tbl(self, sample_id):
        if self.__sample_tbl_mapping is None:
            self.__map_samples()
        if sample_id in self.__sample_tbl_mapping:
            return self.__sample_tbl_mapping[sample_id]
        return None

    def join_all_gtz_anno(self):
        # This function will create a very huge table
        # having all gtz coors, all annotation, and gtz from all samples
        # The rational behind this big join is because report region selection
        # If this method is too time consuming, the other option is to have 
        # views dedicated for reports (one view per one report region),
        # which might be very difficult to maintain
        annovar_info = self.get_annovar_info()
        avdb_info = self.get_avdb_info()
        gtz_tbls = self.get_gtz_tbls()
        # create the table
        tbl_cols = list(VCF_PKEYS)
        tbl_cols += map(lambda x: x.replace("gtz_", "")+"_FILTER",
                        gtz_tbls.keys())
        tbl_cols += reduce(lambda x, y: x+y, annovar_info.values()) 
        tbl_cols += reduce(lambda x, y: x+y, avdb_info.values()) 
        tbl_cols += reduce(lambda x, y: x+y, gtz_tbls.values()) 
        self.drop_table(TBL_NAME_ALL_GTZ_ANNOS)
        self._create_table(tbl_name=TBL_NAME_ALL_GTZ_ANNOS,
                           col_names=tbl_cols,
                           pkeys=VCF_PKEYS,
                           )
        # join and insert
        sql = "INSERT INTO " + TBL_NAME_ALL_GTZ_ANNOS
        sql += " ( " + ", ".join(tbl_cols) + " ) "
        sql += " SELECT"
        sql += " " + ", ".join(map(lambda x: TBL_NAME_GTZ_COORS+"."+x,
                                   VCF_PKEYS))
        sql += ", " + ", ".join(map(lambda x: x+".FILTER", gtz_tbls.keys()))
        for annovar_tbl_name, col_names in annovar_info.items():
            sql += ", " + ", ".join(map(lambda x: annovar_tbl_name+"."+x,
                                        col_names))
        for avdb_tbl_name, col_names in avdb_info.items():
            sql += ", " + ", ".join(map(lambda x: avdb_tbl_name+"."+x,
                                        col_names))
        for gtz_tbl_name, samples_id in gtz_tbls.items():
            sql += ", " + ", ".join(map(lambda x: gtz_tbl_name+"."+x,
                                        samples_id))
        sql += " FROM " + TBL_NAME_GTZ_COORS
        for gtz_tbl_name in gtz_tbls:
            sql += " LEFT OUTER JOIN " + gtz_tbl_name + " ON "
            sql += " AND ".join(map(lambda x: TBL_NAME_GTZ_COORS+"."+x+" = "+gtz_tbl_name+"."+x,
                                    VCF_PKEYS))
        for annovar_tbl_name in annovar_info:
            sql += " LEFT OUTER JOIN " + annovar_tbl_name + " ON "
            sql += " AND ".join(map(lambda x: TBL_NAME_GTZ_COORS+"."+x+" = "+annovar_tbl_name+"."+x,
                                    VCF_PKEYS))
        for avdb_tbl_name in avdb_info:
            sql += " LEFT OUTER JOIN " + avdb_tbl_name + " ON "
            sql += " AND ".join(map(lambda idx: TBL_NAME_GTZ_COORS+"."+VCF_PKEYS[idx]+" = "+avdb_tbl_name+"."+AVDB_PKEYS[idx],
                                    xrange(len(VCF_PKEYS))))
        sql += " ORDER BY"
        sql += " " + ", ".join(map(lambda x: TBL_NAME_GTZ_COORS+"."+x,
                                   VCF_PKEYS))
        return self._exec_sql(sql).rowcount

    def cal_max_ref_maf(self):
        tbl_name = TBL_NAME_ALL_GTZ_ANNOS
        self._add_column(tbl_name, MAX_REF_MAF_COL_NAME)
        sql = "UPDATE " + tbl_name
        sql += " SET " + MAX_REF_MAF_COL_NAME + " = 0"
        self._exec_sql(sql)
        tbl_ref_maf_col_names = []
        for col_name in self.get_col_names(tbl_name):
            if col_name.strip("_") in REF_MAF_COL_NAMES:
                sql = "UPDATE " + tbl_name
                sql += " SET " + MAX_REF_MAF_COL_NAME + " = " + col_name
                sql += " WHERE"
                sql += " " + col_name + " IS NOT NULL"
                sql += " AND " + col_name + " <> ''"
                sql += " AND " + col_name + " > " + MAX_REF_MAF_COL_NAME
                sql += " AND " + col_name + " < 0.5"
                self._exec_sql(sql) 
                sql = "UPDATE " + tbl_name
                sql += " SET " + MAX_REF_MAF_COL_NAME + " = 1 - " + col_name
                sql += " WHERE"
                sql += " " + col_name + " IS NOT NULL"
                sql += " AND " + col_name + " <> ''"
                sql += " AND 1 - " + col_name + " > " + MAX_REF_MAF_COL_NAME
                sql += " AND " + col_name + " > 0.5"
                self._exec_sql(sql) 
                tbl_ref_maf_col_names.append(col_name)

    def set_ref_mutated(self):
        tbl_name = TBL_NAME_ALL_GTZ_ANNOS
        self._add_column(tbl_name, REF_MUTATED_COL_NAME)
        sql = "UPDATE " + tbl_name
        sql += " SET " + REF_MUTATED_COL_NAME + " = 0"
        self._exec_sql(sql)
        sql = "UPDATE " + tbl_name
        sql += " SET " + REF_MUTATED_COL_NAME + " = 1"
        sql += " WHERE " + self.anno_col_to_db_col(GNOMAD_GENOME_ALL_COL_NAME) + " > 0.5"
        self._exec_sql(sql)
