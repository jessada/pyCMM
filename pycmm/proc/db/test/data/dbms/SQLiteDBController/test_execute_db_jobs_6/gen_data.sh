#!/bin/bash

module load annovar
source $PYCMM/bash/cmm_functions.sh
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function extract_avdb_data {
    annovar_file_name=$1
    head -1 "$ANNOVAR_HUMAN_DB_DIR/$annovar_file_name" > "$script_dir/test_$annovar_file_name"
    cmd="grep -P \"^5\\t112154\" $ANNOVAR_HUMAN_DB_DIR/$annovar_file_name >> $script_dir/test_$annovar_file_name"
    eval_cmd "$cmd"
}

extract_avdb_data "hg19_gnomad_genome.txt"
extract_avdb_data "hg19_clinvar_20170130.txt"
extract_avdb_data "hg19_249_Swedes.txt"
extract_avdb_data "hg19_intervar_20170202.txt"
extract_avdb_data "hg19_CMM_Swegen_20161223.txt"

function generate_empty_vcf {
    input_vcf=$1
    output_vcf=$2
    cmd="tabix $input_vcf 5:112154000-112155000 | cut -f1-5 | awk -F'\t' '{ printf \"%s\t.\t.\t.\n\", \$0 }' >> $output_vcf"
    eval_cmd "$cmd"
}

tmp_empty_vcf="$script_dir/tmp_empty.vcf"
> $tmp_empty_vcf

generate_empty_vcf $VCF_CO309_GTZ $tmp_empty_vcf
generate_empty_vcf $VCF_WES294 $tmp_empty_vcf
generate_empty_vcf $VCF_THYRCA $tmp_empty_vcf

tmp_concat_vcf="$script_dir/tmp_concat_vcf.vcf"
echo "##fileformat=VCFv4.2" > $tmp_concat_vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> $tmp_concat_vcf

cmd="grep -v \"^X\" $tmp_empty_vcf | sort | uniq | sort -k1,1n -k2,2n >> $tmp_concat_vcf"
eval_cmd "$cmd"
cmd="grep \"^X\" $tmp_empty_vcf | sort | uniq | sort -k1,1n -k2,2n >> $tmp_concat_vcf"
eval_cmd "$cmd"
cmd="bgzip -f $tmp_concat_vcf"
eval_cmd "$cmd"
cmd="tabix -f -p vcf $tmp_concat_vcf.gz"
eval_cmd "$cmd"

tmp_annovar_prefix="$script_dir/tmp_annovar_prefix"

table_annovar_cmd="table_annovar.pl"
table_annovar_cmd+=" $tmp_concat_vcf.gz"
table_annovar_cmd+=" $ANNOVAR_HUMAN_DB_DIR"
table_annovar_cmd+=" -buildver hg19"
table_annovar_cmd+=" -out $tmp_annovar_prefix"
table_annovar_cmd+=" -protocol refGene,cytoBand,targetScanS,wgRna,tfbsConsSites,exac03constraint,1000g2014oct_all,1000g2014oct_eur"
table_annovar_cmd+=" -operation g,r,r,r,r,r,f,f"
table_annovar_cmd+=" -nastring ."
table_annovar_cmd+=" -vcfinput"
eval_cmd "$table_annovar_cmd"

test_annovar_vcf="$script_dir/test_annovar_vcf.vcf"
cp $tmp_annovar_prefix.hg19_multianno.vcf $test_annovar_vcf
cmd="bgzip -f $test_annovar_vcf"
eval_cmd "$cmd"
cmd="tabix -f -p vcf $test_annovar_vcf.gz"
eval_cmd "$cmd"
rm tmp*

function generate_gtz_vcf {
    input_vcf=$1
    output_vcf=$2

    cmd="tabix -h $input_vcf 5:112154000-112155000 > $output_vcf"
    eval_cmd "$cmd"
    
    bgzip -f $output_vcf
    tabix -f -p vcf $output_vcf.gz
}


test_gtz_Co_309="$script_dir/test_gtz_Co_309.vcf"
generate_gtz_vcf $VCF_CO309_GTZ $test_gtz_Co_309
test_gtz_WES294="$script_dir/test_gtz_WES294.vcf"
generate_gtz_vcf $VCF_WES294 $test_gtz_WES294
test_gtz_THYRCA="$script_dir/test_gtz_THYRCA.vcf"
generate_gtz_vcf $VCF_THYRCA $test_gtz_THYRCA

