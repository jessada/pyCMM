import sys
from xlrd import open_workbook
from xlrd.sheet import ctype_text 
from pycmm.template import pyCMMBase


def check_equal(var1, var2):
    return var1 == var2

def check_in(var1, var2):
    return var1 in var2

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class XlsUtils(pyCMMBase):
    """ general pyCMM template for testing """

    def __init__(self, xls_file):
        pyCMMBase.__init__(self)
        self.__xls_file = xls_file

    @property
    def xls_file(self):
        return self.__xls_file

    @property
    def nsheets(self):
        wb = open_workbook(self.xls_file)
        return len(wb.sheet_names())
    
    def __get_sheet(self, sheet_idx=0):
        wb = open_workbook(self.xls_file)
        return wb.sheet_by_index(sheet_idx)

    def __get_col_idx(self, col_name, ws=None, sheet_idx=0):
        if ws is None:
            ws = self.__get_sheet(sheet_idx)
        header = ws.row(0)
        for idx, cell_obj in enumerate(header):
            if cell_obj.value == col_name:
                return idx
        return -1

    def get_sheet_idx(self, sheet_name):
        wb = open_workbook(self.xls_file)
        sheet_names = wb.sheet_names()
        for sheet_idx in xrange(len(sheet_names)):
            if str(sheet_names[sheet_idx]) == sheet_name:
                return sheet_idx
        return -1
    
    def compare_vals(self,
                     col_name,
                     row_idx1,
                     row_idx2,
                     sheet_idx=0):
        ws = self.__get_sheet(sheet_idx=sheet_idx)
        col_idx = self.__get_col_idx(col_name, ws=ws)
        val1 = ws.cell(row_idx1, col_idx).value
        val2 = ws.cell(row_idx2, col_idx).value
        return val1 == val2
    
    def count_rows(self, sheet_idx=0):
        ws = self.__get_sheet(sheet_idx=sheet_idx)
        return ws.nrows
    
    def count_cols(self, col_name1=None, col_name2=None, sheet_idx=0):
        ws = self.__get_sheet(sheet_idx=sheet_idx)
        if col_name1 is not None:
            col_idx1 = self.__get_col_idx(col_name1, ws=ws)
            col_idx2 = self.__get_col_idx(col_name2, ws=ws)
            return col_idx2 - col_idx1 + 1
        return ws.ncols
