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
-i {file}          specify tabix file (required)
-o {file}          output file name (required)
EOF
)

while getopts ":i:t:o:" OPTION; do
  case "$OPTION" in
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
display_param "tabix file (-i)" "$tabix_file"
display_param "output file (-o)" "$out_file"

## display misc configuration
info_msg
info_msg "misc configuration"
display_param "working direcotry" "$working_dir"

# ****************************************  executing  ****************************************

new_section_txt "E X E C U T I N G"

cmd="tabix -h $tabix_file 1:1-1 > $out_file" 
eval_cmd "$cmd"

printf_phrase="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
param_phrase="\$1, \$2, \$8, \$4, \$5, \$11, \$12, \$13, \$14"

n_samples=`vcf-query -l $tabix_file | wc -l `
for ((n_sample=1; n_sample<=$n_samples; n_sample++));
do
    printf_phrase+="\t%s"
    param_phrase+=", \$$((n_sample+14))"
done;

cmd="convert2annovar.pl"
cmd+=" -format vcf4old"
cmd+=" $tabix_file"
cmd+=" --includeinfo"
cmd+=" | awk -F '\t' '{ printf \"$printf_phrase\n\", $param_phrase}'"
cmd+=" | vcf-sort"
cmd+=" >> $out_file"
eval_cmd "$cmd"

cmd="bgzip -f $out_file"
eval_cmd "$cmd"
cmd="tabix -f -p vcf $out_file.gz"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
