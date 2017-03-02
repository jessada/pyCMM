import re
from pycmm.template import pyCMMBase

ALL_CHROMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "Y"]
CHROM_POS_PATTERN = re.compile(r'''(?P<chrom>.+?):(?P<start_pos>.+?)-(?P<end_pos>.+)''')

class DNARegion(pyCMMBase):
    """  To parse DNA region  """

    def __init__(self, raw_region, **kwargs):
        super(DNARegion, self).__init__(**kwargs)
        self.__parse_region(raw_region)
        self.__raw_region = raw_region

    def get_raw_obj_str(self):
        return {"chromosome": self.chrom,
                "start position": self.start_pos,
                "end postion": self.end_pos,
                }

    @property
    def chrom(self):
        return self.__chrom

    @property
    def start_pos(self):
        return self.__start_pos

    @property
    def end_pos(self):
        return self.__end_pos

    @property
    def region_key(self):
        key = "chr"
        key += self.__chrom
        if self.start_pos is not None:
            key += "_" + self.start_pos
        return key
    
    @property
    def raw_region(self):
        return self.__raw_region

    def __parse_region(self, raw_region):
        match = CHROM_POS_PATTERN.match(raw_region)
        if match is not None:
            self.__chrom = match.group('chrom')
            self.__start_pos = match.group('start_pos')
            self.__end_pos = match.group('end_pos')
        else:
            self.__chrom = raw_region
            self.__start_pos = None
            self.__end_pos = None
