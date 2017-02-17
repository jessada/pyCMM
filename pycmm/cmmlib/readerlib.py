import gzip
from collections import OrderedDict
from pycmm.template import pyCMMBase


class Reader(pyCMMBase):
    """ a base class for read and parse file content """

    def __init__(self, 
                 file_name,
                 parser,
                 **kwargs
                 ):
        self.__file_name = file_name
        # open file to read in text mode
        if file_name.endswith(".gz"):
            self.__reader = gzip.open(file_name, 'rt')
        else:
            self.__reader = open(file_name, 'rt')
        self.__row_count = 0
        self.__parser = parser
        super(Reader, self).__init__(**kwargs)

    def get_raw_obj_str(self):
        raw_repr = OrderedDict()
        raw_repr["file name"] = self.__file_name
        raw_repr["parser"] = self.__parser
        return raw_repr

    # to solve self recalling issues
    @property
    def reader(self):
        return self.__reader

    def __iter__(self):
        return self

    def next(self, mod=1):
        rec = self.__parser(self.__reader.next())
        self.__row_count += 1
        while self.__row_count % mod != 0:
            rec = self.__parser(self.__reader.next())
            self.__row_count += 1
        return rec
