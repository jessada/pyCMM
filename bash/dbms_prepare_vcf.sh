#!/bin/bash 
set -e
#set -o pipefail

source $PYCMM/bash/cmm_functions.sh

module load annovar
script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}           list of vcf input files (comma separated) (required)
-o {dir}            project output root directory (required)
-h                  this help
EOF
)

# parse option
while getopts ":i:o:h" OPTION; do
  case "$OPTION" in
    i)
      comma_separated_vcf_files="$OPTARG"
      ;;
    o)
      project_output_root_dir="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $comma_separated_vcf_files ] || die "list of VCF input files(-i)"
[ ! -z $project_output_root_dir ] || die "project output root dir is required (-o)"

IFS=',' read -ra vcf_files <<< "$comma_separated_vcf_files"

time_stamp=$( date )

working_dir=`mktemp -d`
ta_vcf_out_dir="$project_output_root_dir/ta_vcf"
if [ ! -d "$ta_vcf_out_dir" ]
then
    mkdir -p "$ta_vcf_out_dir"
fi
output_vcf="$ta_vcf_out_dir/all_coors_annotated.vcf"

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

## ****************************************  display configuration  ****************************************
# display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  performing gene based and region based annotation of the variants by doing"
info_msg "  - empty the vcf files (keep only field 1-5)"
info_msg "  - concat the files"
info_msg "  - annotate with ANNOVAR"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "input vcf files (-i)"
for (( n=0; n<$((${#vcf_files[@]})); n++ ))
do
#for vcf_file in "${vcf_files[@]}"; do
#    display_param "  " "$vcf_file"
display_param "  #$((n+1))" "${vcf_files[$n]}"
done

display_param "working directory" "$working_dir"
display_param "output annotated VCF file" "$output_vcf.gz"

# ****************************************  executing  ****************************************

new_section_txt "E X E C U T I N G"

info_msg
info_msg
info_msg ">>> creating empty concat vcf <<<"


tmp_concated_vcf="$working_dir/tmp_concated.vcf"
echo "##fileformat=VCFv4.2" > $tmp_concated_vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> $tmp_concated_vcf

tmp_concated="$working_dir/tmp_concated"
#for chrom in 5
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
    > $tmp_concated
    for vcf_file in "${vcf_files[@]}"
    do
        cmd="tabix -h $vcf_file $chrom"
        cmd+=" | vt decompose -s -"
        cmd+=" | vt normalize -q -r $GRCH37_REF -"
        cmd+=" | grep -P \"\tPASS\t\""
        cmd+=" | grep -Pv \"\t\*\t\""
        cmd+=" | cut -f1-5"
        cmd+=" | awk -F'\t' '{ printf \"%s\t.\t.\t.\n\", \$0 }'"
        cmd+=" >> $tmp_concated"
        eval_cmd "$cmd"
    done
    cmd="sort -k2,2n $tmp_concated | uniq >> $tmp_concated.vcf"
    eval_cmd "$cmd"
done

cmd="bgzip -f $tmp_concated_vcf"
eval_cmd "$cmd"
cmd="tabix -f -p vcf $tmp_concated_vcf.gz"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> running ANNOVAR <<<"

tmp_ta_prefix="$working_dir/tmp_ta"
table_annovar_cmd="table_annovar.pl"
table_annovar_cmd+=" $tmp_concated_vcf.gz"
table_annovar_cmd+=" $ANNOVAR_HUMAN_DB_DIR"
table_annovar_cmd+=" -buildver hg19"
table_annovar_cmd+=" -out $tmp_ta_prefix"
table_annovar_cmd+=" -remove"
table_annovar_cmd+=" -protocol refGene,cytoBand,targetScanS,wgRna,exac03constraint,1000g2014oct_all,1000g2014oct_eur"
table_annovar_cmd+=" -operation g,r,r,r,r,f,f"
table_annovar_cmd+=" -nastring ."
table_annovar_cmd+=" -vcfinput"
eval_cmd "$table_annovar_cmd"

tmp_ta_vcf="$tmp_ta_prefix.hg19_multianno.vcf"
cmd="bgzip -f $tmp_ta_vcf"
eval_cmd "$cmd"
cmd="tabix -f -p vcf $tmp_ta_vcf.gz"
eval_cmd "$cmd"

cmd="tabix -h $tmp_ta_vcf.gz 1:1-1 > $output_vcf" 
eval_cmd "$cmd"

printf_phrase="%s\t%s\t%s\t%s\t%s\t.\t%s\t%s"
param_phrase="\$1, \$2, \$8, \$4, \$5, \$12, \$13"

cmd="convert2annovar.pl"
cmd+=" -format vcf4old"
cmd+=" $tmp_ta_vcf.gz"
cmd+=" --includeinfo"
cmd+=" | awk -F '\t' '{ printf \"$printf_phrase\n\", $param_phrase}'"
cmd+=" | vcf-sort -c"
cmd+=" >> $output_vcf"
eval_cmd "$cmd"

cmd="bgzip -f $output_vcf"
eval_cmd "$cmd"
cmd="tabix -f -p vcf $output_vcf.gz"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
