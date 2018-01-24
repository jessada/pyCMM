#!/bin/bash 
set -o pipefail

module load annovar
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}         input VCF file (required)
-o {file}         output file (required)
-h                this help
EOF
)

# parse option
while getopts ":i:o:h" OPTION; do
  case "$OPTION" in
    i)
      input_vcf="$OPTARG"
      ;;
    o)
      out_file="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $input_vcf ] || die "an input VCF file is required (-i)"
[ ! -z $out_file ] || die "please indicate output file name (-o)"
[ -f "$input_vcf" ] || die "$input_vcf is not a valid file name."

time_stamp=$( date )

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

working_dir=`mktemp -d`

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  decompose and left align vcf file"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "input VCF file (-i)" "$input_vcf"
display_param "output file (-o)" "$out_file"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************
 
new_section_txt "E X E C U T I N G"

function decompose_n_leftalign {
    chrom=$1
    # decompose and leftalign input vcf
    cmd="tabix -h"
    cmd+=" $input_vcf"
    cmd+=" $chrom"
    cmd+=" | vt decompose -s -"
    cmd+=" | vt normalize -q -r $GRCH37_REF -"
    cmd+=" | grep -v \"^#\""
    cmd+=" | grep -Pv \"\t\*\t\""
    cmd+=" | grep -P \"\tPASS\t\""
    cmd+=" >> $tmp_left_align"
    eval_cmd "$cmd"
}

tmp_left_align="$working_dir/tmp_la.vcf"
tabix -h $input_vcf 1:1-1 > $tmp_left_align
decompose_n_leftalign 1
decompose_n_leftalign 2
decompose_n_leftalign 3
decompose_n_leftalign 4
decompose_n_leftalign 5
decompose_n_leftalign 6
decompose_n_leftalign 7
decompose_n_leftalign 8
decompose_n_leftalign 9
decompose_n_leftalign 10
decompose_n_leftalign 11
decompose_n_leftalign 12
decompose_n_leftalign 13 
decompose_n_leftalign 14
decompose_n_leftalign 15
decompose_n_leftalign 16
decompose_n_leftalign 17
decompose_n_leftalign 18
decompose_n_leftalign 19
decompose_n_leftalign 20
decompose_n_leftalign 21
decompose_n_leftalign 22
decompose_n_leftalign "X"
decompose_n_leftalign "Y"

cmd="bgzip -f $tmp_left_align"
eval_cmd "$cmd"

cmd="tabix -p vcf $tmp_left_align.gz"
eval_cmd "$cmd"

cmd="tabix -h $tmp_left_align.gz 1:1-1 > $out_file" 
eval_cmd "$cmd"

# retain all the original information
#printf_phrase="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
#param_phrase="\$1, \$2, \$8, \$4, \$5, \$11, \$12, \$13, \$14"
# suppress QUAL and INFO column
printf_phrase="%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s"
param_phrase="\$1, \$2, \$8, \$4, \$5, \$12, \$14"

n_samples=`vcf-query -l $tmp_left_align.gz | wc -l `
for ((n_sample=1; n_sample<=$n_samples; n_sample++));
do
    printf_phrase+="\t%s"
    param_phrase+=", \$$((n_sample+14))"
done;

cmd="convert2annovar.pl"
cmd+=" -format vcf4old"
cmd+=" $tmp_left_align.gz"
cmd+=" --includeinfo"
cmd+=" | awk -F '\t' '{ printf \"$printf_phrase\n\", $param_phrase}'"
cmd+=" | vcf-sort -c"
cmd+=" >> $out_file"
eval_cmd "$cmd"

cmd="bgzip -f $out_file"
eval_cmd "$cmd"
cmd="tabix -f -p vcf $out_file.gz"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
