from pycmm.template import pyCMMBase

ALL_CHROMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "Y"]

class DNARegion(pyCMMBase):
    """  To parse DNA region  """

    def __init__(self, region):
        pyCMMBase.__init__(self)
        self.__parse_region(region)

    def get_raw_repr(self):
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
    
    def __parse_region(self, region):
        region_items = region.split(":")
        self.__chrom = region_items[0]
        if len(region_items) == 1:
            self.__start_pos = None
            self.__end_pos = None
        else:
            self.__start_pos, self.__end_pos = region_items[1].split("-")
