import sys
import pysam

def create_empty_vcf(vcf_out_file,
                     ):
    f = open(vcf_out_file, "w")
    print >> f, "##fileformat=VCFv4.2"
    print >> f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    f.close()
    pysam.tabix_index(vcf_out_file, preset='vcf')
