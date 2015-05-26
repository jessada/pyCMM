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
-k {file}     known indels (required)
-R {file}     reference file (required)
-o {file}     output target indels interval (required)
-h            this help
EOF
)

# parse option
while getopts ":I:k:R:o:h" OPTION; do
  case "$OPTION" in
    I)
      bam_file="$OPTARG"
      ;;
    k)
      known_indels="$OPTARG"
      ;;
    R)
      ref="$OPTARG"
      ;;
    o)
      out_target_intervals="$OPTARG"
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
[ ! -z $out_target_intervals ] || die "a file to keep list of target intervals is required (-o)"
[ ! -z $ref ] || die "reference file is required (-R)"
[ -f "$bam_file" ] || die "$bam_file is not found"
[ -f "$ref" ] || die "$ref is not found"

IFS=',' read -ra known_indels_files <<< "$known_indels"
time_stamp=$( date )

cd $PYCMM_DIR
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  Create a target list of intervals to be realigned"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM_DIR"
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
display_param "reference sequence (-r)" "$ref"
display_param "ourput target intervals (-o)" "$out_target_intervals"
# ****************************************  executing  ****************************************
# Create a target list of intervals to be realigned
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T RealignerTargetCreator"
cmd+=" -I $bam_file"
cmd+=" -R $ref"
cmd+=" $known_indels_param"
cmd+=" -o $out_target_intervals"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
