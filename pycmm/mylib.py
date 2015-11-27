import sys
import xlrd


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

def count_xls_rows(xls_file_name, sheet_idx=0):
    wb = xlrd.open_workbook(xls_file_name)
    ws = wb.sheet_by_index(sheet_idx)
    return ws.nrows
