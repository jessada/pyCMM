#!/bin/bash
source $PYCMM/bash/cmm_functions.sh
module load annovar

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-k {name}          dataset name (required)
-i {file}          specify tabix file (required)
-o {file}          output file name (required)
EOF
)

while getopts ":k:i:t:o:" OPTION; do
  case "$OPTION" in
    k)
      dataset_name="$OPTARG"
      ;;
    i)
      tabix_file="$OPTARG"
      ;;
    o)
      out_file="$OPTARG"
      ;;
    *)
      die "unrecognized option from executing: $0 $@"
      ;;
  esac
done

[ ! -z $dataset_name ] || die "Please specify dataset name (-k)"
[ ! -z $tabix_file ] || die "Please specify tabix file (-i)"
[ ! -z $out_file ] || die "Plesae specify output file name (-o)"
[ -f $tabix_file ] || die "$tabix_file is not a valid file name"

working_dir=`mktemp -d`

## ****************************************  display configuration  ****************************************
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "parameters"
info_msg "  $params"
info_msg
info_msg "description"
info_msg "  A script to key table to bridge coordinates between vcf and avdb formats:"
info_msg
info_msg
## display required configuration
info_msg "overall configuration"
display_param "dataset name (-k)" "$dataset_name"
display_param "tabix file (-i)" "$tabix_file"
display_param "output file (-o)" "$out_file"

## display misc configuration
info_msg
info_msg "misc configuration"
display_param "working direcotry" "$working_dir"

# ****************************************  executing  ****************************************

new_section_txt "E X E C U T I N G"

tmp_vcf_query=$working_dir/$dataset_name"_tmp_vcf_query"
avinput_file_prefix=$working_dir/$dataset_name"_avdb_individual"
tmp_avdb_uniq=$working_dir/$dataset_name"_tmp_uniq.avdb"
avdb_uniq=$working_dir/$dataset_name".uniq.avdb"

#---------- vcf2avdb --------------
#generate query header
query_header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
query_header+="\tdummy1\tdummy2\tdummy3"
echo -e "$query_header" > $tmp_vcf_query

vcf_query_cmd="vcf-query "
vcf_query_cmd+=" -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT\t1/2\t3/4\t5/6\n'"
vcf_query_cmd+=" $tabix_file"
vcf_query_cmd+=" >> $tmp_vcf_query"
eval_cmd "$vcf_query_cmd"

rm $avinput_file_prefix*
convert2annovar="convert2annovar.pl -format vcf4 $tmp_vcf_query -include --allsample --outfile $avinput_file_prefix"
eval_cmd "$convert2annovar" 

:>$tmp_avdb_uniq
for f in $avinput_file_prefix*
do
    concat_avdb_cmd="cut -f1-11 $f >> $tmp_avdb_uniq"
    eval_cmd "$concat_avdb_cmd" 
done

sort $tmp_avdb_uniq | uniq > $avdb_uniq
#---------- vcf2avdb --------------

#---------- rearrange avdb and add key --------------
:>$out_file
add_key_to_avdb="grep -P \"^[0-9]\" $avdb_uniq"
add_key_to_avdb+=" | awk -F'\t' '{ printf \"%02d_%012d_%s_%s\t%s\t%s\t%s\t%s\t%s\n\", \$6, \$7, \$9, \$10, \$1, \$2, \$3, \$4, \$5 }'"
add_key_to_avdb+=" >> $out_file"
eval_cmd "$add_key_to_avdb"
add_key_to_avdb="grep -vP \"^[0-9]\" $avdb_uniq"
add_key_to_avdb+=" | awk -F'\t' '{ printf \"%s_%012d_%s_%s\t%s\t%s\t%s\t%s\t%s\n\", \$6, \$7, \$9, \$10, \$1, \$2, \$3, \$4, \$5 }'"
add_key_to_avdb+=" >> $out_file"
eval_cmd "$add_key_to_avdb"
#---------- rearrange avdb and add key --------------

new_section_txt "F I N I S H <$script_name>"
