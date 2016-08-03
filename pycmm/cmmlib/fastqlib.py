import sys
from pycmm.template import pyCMMBase
from pycmm.cmmlib.readerlib import Reader

ENCODING_SOLEXA = "Solexa"
ENCODING_ILLUMINA_1_3 = "Illumina-1.3"
ENCODING_ILLUMINA_1_5 = "Illumina-1.5"
ENCODING_ILLUMINA_1_8 = "Illumina-1.8"

ENCODING_RANGES = {
    ENCODING_SOLEXA: (59, 104),
    ENCODING_ILLUMINA_1_3: (64, 104),
    ENCODING_ILLUMINA_1_5: (66, 105),
    ENCODING_ILLUMINA_1_8: (33, 74),
}

class FastqParser(pyCMMBase):
    """ To parse a record in fastq format """

    def __init__(self, rec, **kwargs):
        self.__rec = rec
        super(FastqParser, self).__init__(**kwargs)

    @property
    def phred_score(self):
        if len(self.__rec) < 100:
            return ""
        return self.__rec.strip()

class FastqReader(Reader):
    """ To read and parse fastq file """

    def __init__(self, 
                 **kwargs
                 ):
        kwargs['parser'] = FastqParser
        super(FastqReader, self).__init__(**kwargs)

    def next_phred_score(self):
        return self.next(mod=4)

class Fastq(pyCMMBase):
    """ To read and parse fastq file """

    def __init__(self, 
                 file_name,
                 **kwargs
                 ):
        self.__file_name = file_name
        self.__encoding = None
        super(Fastq, self).__init__(**kwargs)

    @property
    def encoding(self):
        if self.__encoding == None:
            f = FastqReader(file_name=self.__file_name)
            # define negative range
            gmin, gmax  = 99, 0
            while True:
                q_scores = map(lambda x: ord(x),
                               f.next_phred_score().phred_score)
                qmin = min(q_scores)
                qmax = max(q_scores)
                if qmin < gmin or qmax > gmax:
                    # wider the ranges
                    gmin, gmax = min(qmin, gmin), max(qmax, gmax)
                    # get possible encoding ranges
                    valid = []
                    for encoding, (emin, emax) in ENCODING_RANGES.items():
                        if gmin >= emin and gmax <= emax:
                            valid.append(encoding)
                    if len(valid) == 0:
                        self.warning("no encodings for range: %s for file %s" % (str((gmin, gmax)), self.__file_name) )
                        return ""
                    if len(valid) == 1:
                        self.__encoding = valid[0]
                        return self.__encoding
            self.throw("no encodings for range: %s" % str((gmin, gmax)))
            sys.exit()
        return self.__encoding
