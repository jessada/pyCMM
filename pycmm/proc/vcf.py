from pycmm.utils import mylogger
from pycmm.utils import exec_sh
from pycmm.utils import concat_files
from pycmm.proc.avdb import uniq_avdb

def vcf2pavdb(vcf_in_file,
              pavdb_out_file,
              ):
    mylogger.info("inside proc.vcf2pavdb")
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
    mylogger.info("inside proc.pavdb2avinputs")
    cmd = "convert2annovar.pl"
    cmd += " -format vcf4"
    cmd += " " + pavdb_in_file
    cmd += " -include"
    cmd += " --allsample"
    cmd += " --outfile " + avinput_file_prefix
    exec_sh(cmd)

def vcf2avdb(vcf_in_file,
             avdb_out_file,
             working_dir,
# The following params are not being used in this version
# They are just here for reminding that I might need to implement them
             log_file=None,
             ):
    mylogger.info("inside proc.vcf2avdb")
    tmp_pavdb_file = vcf_in_file + 'tmp_pavdb'
    vcf2pavdb(vcf_in_file,
              tmp_pavdb_file
              )
    tmp_avinput_file_prefix = vcf_in_file
    pavdb2avinputs(tmp_pavdb_file,
                   tmp_avinput_file_prefix
                   )
    tmp_concat_file = vcf_in_file + 'raw.avdb'
    concat_files(tmp_avinput_file_prefix+'*.avinput',
                 tmp_concat_file
                 )
    uniq_avdb(tmp_concat_file, avdb_out_file)

