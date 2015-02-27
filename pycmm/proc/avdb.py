from pycmm.utils import mylogger
from pycmm.utils import exec_sh

def uniq_avdb(avdb_in_file,
              uniq_avdb_out_file,
              ):
    mylogger.info("inside proc.uniq_avdb")
    cmd = "cut -f1-11 " + avdb_in_file
    cmd += " | sort"
    cmd += " | uniq"
    cmd += " >> " + uniq_avdb_out_file
    exec_sh(cmd)
