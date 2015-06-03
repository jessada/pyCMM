#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-I {file}     input bam file to be realigned(required)
-R {file}     reference file (required)
-t {file}     indels interval (required)
-k {file}     known indels (required)
-o {file}     output re-aligned reads file (required)
-h            this help
EOF
)

# parse option
while getopts ":I:R:t:k:o:h" OPTION; do
  case "$OPTION" in
    I)
      bam_file="$OPTARG"
      ;;
    R)
      ref="$OPTARG"
      ;;
    t)
      target_intervals="$OPTARG"
      ;;
    k)
      known_indels="$OPTARG"
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

[ ! -z $bam_file ] || die "a bam file to be re-aligned is required (-I)"
[ ! -z $known_indels ] || die "list of known indels in vcf format is required (-k)"
[ ! -z $target_intervals ] || die "a list of target intervals is required (-t)"
[ ! -z $out_file ] || die "output file name is required (-o)"
[ ! -z $ref ] || die "reference file is required (-R)"
[ -f "$bam_file" ] || die "$bam_file is not found"
[ -f "$ref" ] || die "$ref is not found"

IFS=',' read -ra known_indels_files <<< "$known_indels"
time_stamp=$( date )

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  Perform realignment of the target intervals"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "bam file (-I)" "$bam_file"
display_param "known indels (-k)"
known_indels_param=""
for (( indel_idx=0; indel_idx<$((${#known_indels_files[@]})); indel_idx++ ))
do
    display_param "  known indels file #$(( indel_idx+1 )) " "${known_indels_files[$indel_idx]}"
    known_indels_param+=" -known ${known_indels_files[$indel_idx]}"
done
display_param "reference sequence (-R)" "$ref"
display_param "target intervals (-t)" "$target_intervals"
display_param "output bam file with indel realigned (-o)" "$out_file"
# ****************************************  executing  ****************************************
#Perform realignment of the target intervals
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T IndelRealigner"
cmd+=" -R $ref"
cmd+=" -I $bam_file"
cmd+=" -targetIntervals $target_intervals"
cmd+=" $known_indels_param"
cmd+=" -o $out_file"
eval_cmd "$cmd"
new_section_txt "F I N I S H <$script_name>"

