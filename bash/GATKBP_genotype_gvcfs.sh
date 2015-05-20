#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-G {file}         a file containing list of gvcf files to be genotyped (required)
-o {file}         genotyped GVCF file (required)
-r {file}         reference file
-h                this help
EOF
)

# parse option
while getopts ":G:o:r:h" OPTION; do
  case "$OPTION" in
    G)
      gvcf_list="$OPTARG"
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

[ ! -z $gvcf_list ] || die "a list of gvcf files to be genotyped is required (-G)"
[ ! -z $out_file ] || die "please indicate output genotyped gvcf file name (-o)"
[ ! -z $ref ] || die "reference file is required (-r)"
[ -f "$ref" ] || die "reference $ref is not found"
[ -f "$gvcf_list" ] || die "$gvcf_list is not a valid file name. List of gvcf files must be put togetther in a file"

time_stamp=$( date )
IFS=',' read -ra gvcf_files <<< "$gvcf_list"

cd $PYCMM_DIR
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  Perform joint genotyping on gVCF files produced by HaplotypeCaller"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM_DIR"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "input gvcf files (-G):"
variant_cmd=""
gvcf_idx=1
while read gvcf; do
    display_param "  gvcf file #$gvcf_idx " "$gvcf"
    variant_cmd+=" --variant $gvcf"
    ((gvcf_idx++))
done < "$gvcf_list"
display_param "output genotyped GVCF file (-o)" "$out_file"
display_param "reference file (-r)" "$ref"

# ****************************************  executing  ****************************************
 
cmd="java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T GenotypeGVCFs"
cmd+=" $variant_cmd"
cmd+=" -R $ref" 
cmd+=" -o $out_file"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
