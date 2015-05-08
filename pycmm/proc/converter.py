import sys
import glob
from os.path import basename
from os.path import join as join_path
from pycmm.utils import mylogger
from pycmm.utils import exec_sh
from pycmm.utils import concat_files
from pycmm.proc.avdb import uniq_avdb
from pycmm.proc.avdb import RawAVDBRecord
from pycmm.proc.avdb import KeyAVDBRecord

def vcf2pavdb(vcf_in_file,
              pavdb_out_file,
              ):
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    f = open(pavdb_out_file, "w")
    print >> f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdummy1\tdummy2"
    f.close()
    cmd = "vcf-query"
    cmd += " -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT\t1/2\t3/4\n' "
    cmd += " " + vcf_in_file
    cmd += " >> " + pavdb_out_file
    exec_sh(cmd)

def pavdb2avinputs(pavdb_in_file,
                   avinput_file_prefix,
                   ):
    """
    using convert2annovar.pl to convert mutation from vcf format into raw
    avdb (avinput) format

    P L S   B E   A W A R E   T H A T   K E Y S   M I G H T   B E   W R O N G
    """
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    cmd = "convert2annovar.pl"
    cmd += " --format vcf4"
    cmd += " " + pavdb_in_file
    cmd += " --include"
    cmd += " --allsample"
    cmd += " --outfile " + avinput_file_prefix
    exec_sh(cmd)

def raw_avdb2key_avdb(raw_avdb_in_file,
                      key_avdb_out_file,
                      ):
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    f_in = open(raw_avdb_in_file, "r")
    f_out = open(key_avdb_out_file, "w")
    for line in f_in:
        avdb_rec = RawAVDBRecord(line)
        buffer = avdb_rec.avdb_chrom
        buffer += "\t" + avdb_rec.avdb_start_pos
        buffer += "\t" + avdb_rec.avdb_end_pos
        buffer += "\t" + avdb_rec.avdb_ref
        buffer += "\t" + avdb_rec.avdb_alt
        buffer += "\t" + avdb_rec.vcf_key
        buffer += "\t" + avdb_rec.chrom
        buffer += "\t" + avdb_rec.pos
        buffer += "\t" + avdb_rec.ref
        buffer += "\t" + avdb_rec.alt
        f_out.write(buffer + '\n')
    f_in.close()
    f_out.close()

def vcf2avdb(vcf_in_file,
             avdb_out_file,
             working_dir,
# The following params are not being used in this version
# They are just here for reminding that I might need to implement them
             log_file=None,
             ):
    """
    the expected format of output avdb file is entirely based on "-include"
    parameter of "convert2annovar.pl" command. It basically has the first
    five columns in avdb format and the rest are the raw content from vcf
    (chrom, pos, id, ref, alt, qual, filter, info, and the rest)
    """
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    tmp_file_prefix = join_path(working_dir,
                                basename(vcf_in_file))
    tmp_pavdb_file = tmp_file_prefix + '.tmp.avdb'
    vcf2pavdb(vcf_in_file,
              tmp_pavdb_file
              )
    tmp_avinput_file_prefix = tmp_file_prefix
    pavdb2avinputs(tmp_pavdb_file,
                   tmp_avinput_file_prefix
                   )
    avinput_count = 0
    key_avdb_file_fmt = '{prefix}.key_{file_no:03d}.avdb'
    key_avdb_files = []
    for raw_avdb_file in glob.glob(tmp_avinput_file_prefix+'*.avinput'):
        avinput_count += 1
        tmp_key_avdb_file = key_avdb_file_fmt.format(prefix=tmp_file_prefix,
                                                     file_no=avinput_count,
                                                     )
        raw_avdb2key_avdb(raw_avdb_file, tmp_key_avdb_file)
        key_avdb_files.append(tmp_key_avdb_file)
    tmp_concat_file = tmp_file_prefix + '.concat.avdb'
    concat_files(key_avdb_files,
                 tmp_concat_file
                 )
    uniq_avdb(tmp_concat_file, avdb_out_file)

def avdb2bed(avdb_in_file,
             bed_out_file,
             working_dir,
# The following params are not being used in this version
# They are just here for reminding that I might need to implement them
             log_file=None,
             ):
    """
    the expected format of output bed file is
    - The first three columns are converted vcf coordinate to be used with
      UCSC liftover
    - The forth column is vcf mutation key (pos_chrom_ref_alt)
    - The fifth to eighth columns are the each component of the key
    """
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    f_in = open(avdb_in_file, "r")
    f_out = open(bed_out_file, "w")
    for line in f_in:
        avdb_rec = KeyAVDBRecord(line)
        if avdb_rec.chrom.startswith("chr"):
            buffer = avdb_rec.avdb_chrom
        else:
            buffer = "chr" + avdb_rec.avdb_chrom
        buffer += "\t" + avdb_rec.avdb_start_pos
        buffer += "\t" + str(int(avdb_rec.avdb_start_pos)+1)
        buffer += "\t" + avdb_rec.vcf_key 
        f_out.write(buffer + '\n')
    f_in.close()
    f_out.close()

