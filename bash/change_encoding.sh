#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-I {file}     input fastq file in gz format (required)
-o {file}     output directory (required)
-h            this help
EOF
)

# parse option
while getopts ":I:o:h" OPTION; do
  case "$OPTION" in
    I)
      fastq_gz_file="$OPTARG"
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

[ ! -z $fastq_gz_file ] || die "a fastq.gz file to be re-encoded is required (-I)"
[ ! -z $output_dir ] || die "a output directory is required (-w)"
[ -f "$fastq_gz_file" ] || die "$fastq_gz_file is not found"
[ -d "$output_dir" ] || die "$output_dir is not found"

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

short_fastq_gz_file=$(basename "$fastq_gz_file")
fastq_file="$output_dir/${short_fastq_gz_file%.gz}"
encoded_file="${fastq_file%.fastq}.fix.fastq"

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  a script to re-encode fastq file from Illumina 1.3+ to Illumina 1.8+"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "input file (-I)" "$fastq_gz_file"
display_param "output direcotry (-w)" "$output_dir"
display_param "  output file" "$encoded_file".gz
# ****************************************  executing  ****************************************
cmd="cp $fastq_gz_file $output_dir"
eval_cmd "$cmd"

cmd="gunzip $output_dir/$short_fastq_gz_file"
eval_cmd "$cmd"

cmd="seqret fastq-illumina::$fastq_file fastq::$encoded_file"
#info_msg "$cmd"
eval_cmd "$cmd"

rm $fastq_file
cmd="gzip $encoded_file"
#info_msg "$cmd"
eval_cmd "$cmd"
new_section_txt "F I N I S H <$script_name>"

