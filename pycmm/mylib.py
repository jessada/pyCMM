import sys
import xlrd


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
