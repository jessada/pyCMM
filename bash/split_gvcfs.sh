#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {dir}            gvcf files directory (required)
-r {file}           list of chromosome regions (required)
-o {dir}            output split gvcf files directory (required)
-h                  this help
EOF
)

# parse option
while getopts ":i:r:o:h" OPTION; do
  case "$OPTION" in
    i)
      input_dir="$OPTARG"
      ;;
    r)
      regions_list="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $input_dir ] || die "gvcf files directory is required (-i)"
[ ! -z $regions_list ] || die "chromosome regions list is required (-r)"
[ ! -z $output_dir ] || die "output split gvcf files directory is required (-o)"
[ -d "$input_dir" ] || die "$input_dir is not found"
[ -f "$regions_list" ] || die "$regions_list is not found"

if [ ! -d "$output_dir" ]
then
    mkdir -p "$output_dir"
fi

time_stamp=$( date )
if [ -z $working_dir ]
then
    working_dir=`mktemp -d`
fi

tmp_bwa_mem_out="$working_dir/sample.sam"
tmp_sorted_bam="$working_dir/sample.bam"
tmp_dedup_bam="$working_dir/sample.dedup.bam"
tmp_metrics_file="$working_dir/sample.metics.txt"
tmp_realign_target_intevals="$working_dir/sample.target.intervals"
tmp_realigned_bam="$working_dir/sample.dedup.realigned.bam"
tmp_recal_table="$working_dir/sample.recal.table"

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  Split gvcf files by region"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "input gvcf files directory (-i)" "$input_dir"
display_param "chromosome regions list (-r)" "$regions_list"
display_param "output split gvcf files direcotry (-o)" "$output_dir"

# ****************************************  executing  ****************************************

new_section_txt "E X E C U T I N G"

while read region
do
    region_txt=`echo $region | sed 's/:/_/' | sed 's/-/_/'`
    region_output_dir="$output_dir/$region_txt"
    if [ ! -d "$region_output_dir" ]
    then
        mkdir -p $region_output_dir
    fi
    for gvcf_gz_full_file in $input_dir/*.vcf.gz
    do
        gvcf_gz_file_name=$(basename "$gvcf_gz_full_file")
        sample_prefix="${gvcf_gz_file_name%.g.vcf.gz}"
        output_vcf_full_file="$region_output_dir/$sample_prefix.$region_txt.g.vcf"
    
        cmd="tabix -h"
        cmd+=" $gvcf_gz_full_file"
        cmd+=" $region"
        cmd+=" > $output_vcf_full_file"
        eval_cmd "$cmd"
    
        cmd="bgzip -f"
        cmd+=" $output_vcf_full_file"
        eval_cmd "$cmd"
    
        cmd="tabix"
        cmd+=" -p vcf"
        cmd+=" $output_vcf_full_file.gz"
        eval_cmd "$cmd"
    done
done < $regions_list

new_section_txt "F I N I S H <$script_name>"
