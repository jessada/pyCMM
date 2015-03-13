import sys
from collections import OrderedDict
from pycmm.template import pyCMMBase
from pycmm.settings import VCF_MUT_KEY_INT_FMT
from pycmm.settings import VCF_MUT_KEY_STR_FMT
from pycmm.utils import mylogger
from pycmm.utils import exec_sh

RAW_AVDB_0_IDX_AVDB_CHROM = 0
RAW_AVDB_0_IDX_AVDB_START_POS = 1
RAW_AVDB_0_IDX_AVDB_END_POS = 2
RAW_AVDB_0_IDX_AVDB_REF = 3
RAW_AVDB_0_IDX_AVDB_ALT = 4
RAW_AVDB_0_IDX_CHROM = 5
RAW_AVDB_0_IDX_POS = 6
RAW_AVDB_0_IDX_REF = 8
RAW_AVDB_0_IDX_ALT = 9

KEY_AVDB_0_IDX_VCF_KEY = 5
KEY_AVDB_0_IDX_CHROM = 6
KEY_AVDB_0_IDX_POS = 7
KEY_AVDB_0_IDX_REF = 8
KEY_AVDB_0_IDX_ALT = 9

class RawAVDBRecord(pyCMMBase):
    """ to automatically parse content from raw AVDB (avinputs) file """

    def __init__(self, rec):
        if isinstance(rec, list):
            self.rec = rec
        else:
            self.rec = rec.split('\t')

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr['AVDB CHROM'] = self.avdb_chrom
        raw_repr['AVDB START POS'] = self.avdb_start_pos
        raw_repr['AVDB END POS'] = self.avdb_end_pos
        raw_repr['AVDB REF'] = self.avdb_ref
        raw_repr['AVDB ALT'] = self.avdb_alt
        raw_repr['VCF KEY'] = self.vcf_key
        raw_repr['CHROM'] = self.chrom
        raw_repr['POS'] = self.pos
        raw_repr['REF'] = self.ref
        raw_repr['ALT'] = self.alt
        return raw_repr

    @property
    def avdb_chrom(self):
        """ AVDB CHROM """
        return self.rec[RAW_AVDB_0_IDX_AVDB_CHROM]

    @property
    def avdb_start_pos (self):
        """ AVDB START POS """
        return self.rec[RAW_AVDB_0_IDX_AVDB_START_POS]

    @property
    def avdb_end_pos(self):
        """ AVDB END POS """
        return self.rec[RAW_AVDB_0_IDX_AVDB_END_POS]

    @property
    def avdb_ref(self):
        """ AVDB REF """
        return self.rec[RAW_AVDB_0_IDX_AVDB_REF]

    @property
    def avdb_alt(self):
        """ AVDB ALT """
        return self.rec[RAW_AVDB_0_IDX_AVDB_ALT]

    @property
    def vcf_key(self):
        """ VCF KEY """
        if self.chrom.isdigit():
            return VCF_MUT_KEY_INT_FMT.format(chrom=int(self.chrom),
                                              pos=int(self.pos),
                                              ref=self.ref,
                                              alt=self.alt,
                                              )
        else:
            return VCF_MUT_KEY_STR_FMT.format(chrom=self.chrom,
                                              pos=int(self.pos),
                                              ref=self.ref,
                                              alt=self.alt,
                                              )

    @property
    def chrom(self):
        """ CHROM """
        return self.rec[RAW_AVDB_0_IDX_CHROM]

    @property
    def pos(self):
        """ POS """
        return self.rec[RAW_AVDB_0_IDX_POS]

    @property
    def ref(self):
        """ REF """
        return self.rec[RAW_AVDB_0_IDX_REF]

    @property
    def alt(self):
        """ ALT """
        return self.rec[RAW_AVDB_0_IDX_ALT]


class KeyAVDBRecord(RawAVDBRecord):
    """ to automatically parse content from KeyAVDB file """

    def __init__(self, *args):
        RawAVDBRecord.__init__(self, *args)

    @property
    def vcf_key(self):
        """ VCF KEY """
        return self.rec[KEY_AVDB_0_IDX_VCF_KEY]

    @property
    def chrom(self):
        """ CHROM """
        return self.rec[KEY_AVDB_0_IDX_CHROM]

    @property
    def pos(self):
        """ POS """
        return self.rec[KEY_AVDB_0_IDX_POS]

    @property
    def ref(self):
        """ REF """
        return self.rec[KEY_AVDB_0_IDX_REF]

    @property
    def alt(self):
        """ ALT """
        return self.rec[KEY_AVDB_0_IDX_ALT]


def uniq_avdb(avdb_in_file,
              uniq_avdb_out_file,
              ):
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    #mylogger.debug("inside proc.uniq_avdb")
    cmd = "cut -f1-11 " + avdb_in_file
    cmd += " | sort"
    cmd += " | uniq"
    cmd += " | sort -k6,6"
    cmd += " > " + uniq_avdb_out_file
    exec_sh(cmd)
