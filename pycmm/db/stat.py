import sys
from pycmm.utils import mylogger
import vcf

def vcf_query(vcf_in_file,
              vcf_chr=None,
              vcf_start_pos=None,
              vcf_end_pos=None,
              ):
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    if vcf_chr is None:
        vcf_reader = vcf.Reader(filename=vcf_in_file)
    elif vcf_start_pos is None:
        vcf_reader = vcf.Reader(filename=vcf_in_file).fetch(vcf_chr)
    else:
        vcf_reader = vcf.Reader(filename=vcf_in_file).fetch(vcf_chr,
                                                            int(vcf_start_pos),
                                                            int(vcf_end_pos))

    return vcf_reader
