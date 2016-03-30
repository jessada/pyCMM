from collections import OrderedDict
from xlsxwriter.worksheet import Worksheet
from xlsxwriter import Workbook
from pycmm.template import pyCMMBase
from pycmm.cmmlib.colorlib import COLORS_RGB

NO_COLOR = 'no_color'


class CMMWorksheet(Worksheet, pyCMMBase):
    """ a customized inherited from xlsxwriter.worksheet.Worksheet """

    def __init__(self, 
                 dflt_fmt=None,
                 **kwargs
                 ):
        self.__dflt_fmt = dflt_fmt
        super(CMMWorksheet, self).__init__(**kwargs)
        pyCMMBase.__init__(self)

    def write(self, row, col, val, cell_format=None):
        if (cell_format is None) and (self.__dflt_fmt is not None):
            cell_format = self.__dflt_fmt
        super(CMMWorksheet, self).write(row, col, val, cell_format) 

class CMMWorkbook(Workbook, pyCMMBase):
    """ a customized workbook inherited from xlsxwriter.Workbook """

    def __init__(self, 
                 **kwargs
                 ):
        super(CMMWorkbook, self).__init__(**kwargs)
        pyCMMBase.__init__(self)
        self.__dflt_hash_fmt = {'font_name': 'Arial', 'font_size': 9}
        self.__dflt_fmt = self.add_format(self.__dflt_hash_fmt)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
#        raw_repr["file name"] = self.__file_name
#        raw_repr["parser"] = self.__parser
        return raw_repr

    def _add_sheet(self, name, is_chartsheet):
        # Utility for shared code in add_worksheet() and add_chartsheet().
        sheet_index = len(self.worksheets_objs)
        name = self._check_sheetname(name, is_chartsheet)
        # Initialisation data to pass to the worksheet.
        init_data = {
            'name': name,
            'index': sheet_index,
            'str_table': self.str_table,
            'worksheet_meta': self.worksheet_meta,
            'optimization': self.optimization,
            'tmpdir': self.tmpdir,
            'date_1904': self.date_1904,
            'strings_to_numbers': self.strings_to_numbers,
            'strings_to_formulas': self.strings_to_formulas,
            'strings_to_urls': self.strings_to_urls,
            'default_date_format': self.default_date_format,
            'default_url_format': self.default_url_format,
            'excel2003_style': self.excel2003_style,
        }
        if is_chartsheet:
            worksheet = Chartsheet()
        else:
            worksheet = CMMWorksheet(dflt_fmt=self.__dflt_fmt)
        worksheet._initialize(init_data)
        self.worksheets_objs.append(worksheet)
        self.sheetnames.append(name)
        return worksheet

    def add_colors_format(self, added_hash_fmt):
        fmts = OrderedDict()
        base_hash_fmt = self.__dflt_hash_fmt.copy()
        for key in added_hash_fmt:
            base_hash_fmt[key] = added_hash_fmt[key]
        fmts[NO_COLOR] = self.add_format(base_hash_fmt)
        for color in COLORS_RGB:
            color_hash_fmt = base_hash_fmt.copy()
            color_hash_fmt['bg_color'] = COLORS_RGB[color]
            fmts[color] = self.add_format(color_hash_fmt)
        return fmts

#    def init_colors_formats(self):
#        dflt_cell_hash_fmt = self.__dflt_hash_fmt.copy()
#        self.__cell_fmts = self.__init_colors_format(dflt_cell_hash_fmt)

#    @property
#    def reader(self):
#        return self.__reader
#
#    def __iter__(self):
#        return self
#
#    def next(self):
#        return self.__parser(self.__reader.next())

#class CellFormatManager(PlinkBase):
#    """ A class to define cell format property """
#
#    def __init__(self, work_book, color_dict):
#        self.__wb = work_book
#        self.__color_dict = color_dict
#        self.__dflt_hash_fmt = {'font_name': 'Arial', 'font_size': 9}
#        self.__dflt_fmt = self.__add_fmt(self.__dflt_hash_fmt)
#        self.__init_colors_formats()
#
#    def get_raw_repr(self):
#        return {"color dict": self.__color_dict,
#                "number of color": self.n_colors,
#                }
#
#    def __add_fmt(self, fmt_dict):
#        return self.__wb.add_format(fmt_dict)
#
#    @property
#    def default_format(self):
#        return self.__dflt_fmt
#
#    def __init_colors_format(self, default_hash_format):
#        fmts = OrderedDict()
#        fmts[DFLT_FMT] = self.__add_fmt(default_hash_format)
#        colors = self.__color_dict.keys()
#        for color_idx in xrange(len(colors)):
#            color = colors[color_idx]
#            color_hash_fmt = default_hash_format.copy()
#            color_hash_fmt['bg_color'] = self.__color_dict[color]
#            fmts[color_idx] = self.__add_fmt(color_hash_fmt)
#            fmts[color] = fmts[color_idx]
#        return fmts
#
#    def __init_colors_formats(self):
#        dflt_bp_hash_fmt = self.__dflt_hash_fmt.copy()
#        dflt_bp_hash_fmt['align'] = 'center'
#        self.__bp_fmts = self.__init_colors_format(dflt_bp_hash_fmt)
#        dflt_stat_hash_fmt = self.__dflt_hash_fmt.copy()
#        dflt_stat_hash_fmt['rotation'] = 90
#        self.__stat_fmts = self.__init_colors_format(dflt_stat_hash_fmt)
#        dflt_snp_hash_fmt = self.__dflt_hash_fmt.copy()
#        self.__snp_fmts = self.__init_colors_format(dflt_snp_hash_fmt)
#
#    @property
#    def n_colors(self):
#        return len(self.__color_dict)
#
#    @property
#    def bp_fmts(self):
#        return self.__bp_fmts
#
#    @property
#    def stat_fmts(self):
#        return self.__stat_fmts
#
#    @property
#    def snp_fmts(self):
#        return self.__snp_fmts

