from pycmm.template import pyCMMBase


class HapAssocUtils(pyCMMBase):
    """ Encapsulated util to handle excel sheet """

    def __init__(self, hap_assoc_file):
        pyCMMBase.__init__(self)
        self.__hap_assoc_file = hap_assoc_file

    @property
    def hap_assoc_file(self):
        return self.__hap_assoc_file

    @property
    def nlines(self):
        with open(self.hap_assoc_file, 'r') as f:
            nlines = 0
            for line in f:
                nlines += 1
        return nlines-1

    @property
    def ncols(self):
        with open(self.hap_assoc_file, 'r') as f:
            return len(f.next().strip().split())
