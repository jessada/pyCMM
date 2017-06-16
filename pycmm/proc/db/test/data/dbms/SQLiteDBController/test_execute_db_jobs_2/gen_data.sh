#!/bin/bash

source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function extract_avdb_data {
    annovar_file_name=$1
    head -1 "$ANNOVAR_HUMAN_DB_DIR/$annovar_file_name" > "$script_dir/test_$annovar_file_name"
    cmd="grep -P \"^10\\t8968\" $ANNOVAR_HUMAN_DB_DIR/$annovar_file_name >> $script_dir/test_$annovar_file_name"
    eval_cmd "$cmd"
}


#extract_avdb_data "hg19_gnomad_genome.txt"
#extract_avdb_data "hg19_clinvar_20150330.txt"
#extract_avdb_data "hg19_clinvar_20170130.txt"
#extract_avdb_data "hg19_CMM_Axeq_chr5_19.txt"
#extract_avdb_data "hg19_CMM_OAF_BrC_CRC_prostate.txt"
#extract_avdb_data "hg19_CMM_OAF_familial_CRCs.txt"
#extract_avdb_data "hg19_CMM_OAF_CHEK2.txt"
#extract_avdb_data "hg19_CMM_Swegen_20161223.txt"
#extract_avdb_data "hg19_exac03constraint.txt"
#extract_avdb_data "hg19_spidex.txt"
#extract_avdb_data "hg19_intervar.txt"
#extract_avdb_data "hg19_AFR.sites.2014_10.txt"
#extract_avdb_data "hg19_ljb26_all.txt"
#extract_avdb_data "hg19_dbnsfp30a.txt"
#extract_avdb_data "hg19_avsnp144.txt"
#


