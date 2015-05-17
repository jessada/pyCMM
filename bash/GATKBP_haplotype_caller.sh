#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-B {BAM file}       BAM file (required)
-o {file}           output GVCF file (required)
-r {file}           reference file
-h                  this help
EOF
)

# parse option
while getopts ":B:o:r:h" OPTION; do
  case "$OPTION" in
    B)
      bam_file="$OPTARG"
      ;;
    o)
      out_file="$OPTARG"
      ;;
    r)
      ref="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $bam_file ] || die "a BAM input file is required (-B)"
[ ! -z $out_file ] || die "please indicate output file name (-o)"
[ ! -z $ref ] || die "reference file is required (-r)"
[ -f "$bam_file" ] || die "$bam_file is not found"
[ -f "$ref" ] || die "$ref is not found"

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
info_msg "  Mark duplicates and index the BAM file"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM_DIR"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "BAM input file (-B)" "$bam_file"
display_param "output GVCF file (-o)" "$out_file"
display_param "reference file (-r)" "$ref"

# ****************************************  executing  ****************************************
 
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T HaplotypeCaller"
cmd+=" -R $ref"
cmd+=" -I $bam_file"
cmd+=" --emitRefConfidence GVCF"
cmd+=" --variant_index_type LINEAR"
cmd+=" --variant_index_parameter 128000"
cmd+=" -o $out_file"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
