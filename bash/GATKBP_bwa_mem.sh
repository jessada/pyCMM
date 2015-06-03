#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-I {fastq list}     list of fastq files (comma separated) (required)
-R {file}           reference file
-g {sample group}   sample group (required) 
-o {file}           output aligned reads file (required)
-h                  this help
EOF
)

# parse option
while getopts ":I:R:g:o:h" OPTION; do
  case "$OPTION" in
    I)
      sample_list="$OPTARG"
      ;;
    R)
      ref="$OPTARG"
      ;;
    g)
      sample_group="$OPTARG"
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

[ ! -z $sample_list ] || die "list of fastq files is required (-I)"
[ ! -z $sample_group ] || die "sample group is required (-g)"
[ ! -z $out_file ] || die "output file name is required (-o)"
[ ! -z $ref ] || die "reference file is required (-R)"
IFS=',' read -ra fastq_files <<< "$sample_list"
sample1="${fastq_files[0]}"
sample2="${fastq_files[1]}"
[ -f "$sample1" ] || die "$sample1 is not found"
[ -f "$sample2" ] || die "$sample2 is not found"
[ -f "$ref" ] || die "$ref is not found"

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
info_msg "  Generate a SAM file containing aligned reads"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "sample prefix (-I)" "$sample_prefix"
display_param "  sample 1" "$sample1"
display_param "  sample 2" "$sample2"
display_param "sample group (-g)" "$sample_group"
display_param "output directory (-o)" "$out_file"
display_param "reference file (-R)" "$ref"

# ****************************************  executing  ****************************************
sample_file=$(basename "$sample1")
sample_name=${sample_file%_1.fastq.gz}
read_group="$sample_group.$sample_name"

cmd=" bwa mem -M -R '@RG\tID:$read_group\tSM:$sample_name\tPL:illumina\tLB:lib1\tPU:unit1'"
cmd+=" $ref"
cmd+=" $sample1"
cmd+=" $sample2"
cmd+=" > $out_file" 
eval_cmd "$cmd"
new_section_txt "F I N I S H <$script_name>"

