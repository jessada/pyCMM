#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}            file containg list of VCF files to be concatenated (required)
-o {file}            output concatenated VCF file (required)
-h                  this help
EOF
)

# parse option
while getopts ":i:o:h" OPTION; do
  case "$OPTION" in
    i)
      vcfs_list_file="$OPTARG"
      ;;
    o)
      output_vcf_gz="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $vcfs_list_file ] || die "list of VCF files is required (-i)"
[ ! -z $output_vcf_gz ] || die "output VCF file is required (-o)"
[ -f "$vcfs_list_file" ] || die "$vcfs_list_file is not found"

time_stamp=$( date )

working_dir=`mktemp -d`
output_vcf="$working_dir/vcf_concated.vcf"

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  Concatenate VCF files"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "list of VCF files to be concatenated (-i)" "$vcfs_list_file"
display_param "output concatenated VCF file (-o)" "$output_vcf_gz"

# ****************************************  executing  ****************************************

new_section_txt "E X E C U T I N G"

# create vcf header
output_vcf=${output_vcf_gz%.gz}
vcf_file=`head -1 $vcfs_list_file`
cmd="tabix -h "
cmd+=" $vcf_file"
cmd+=" 1:0-0"
cmd+=" > $output_vcf"
eval_cmd "$cmd"

while read vcf_file
do
#    info_msg $vcf_file
    cmd="zgrep -v \"^#\""
    cmd+=" $vcf_file"
    cmd+=" >> $output_vcf"
    eval_cmd "$cmd"
done < $vcfs_list_file

#cmd=" mv $output_vcf"
#cmd+=" $output_vcf"
#eval_cmd "$cmd"
#
cmd=" bgzip"
cmd+=" -f $output_vcf"
eval_cmd "$cmd"

cmd=" tabix"
cmd+=" -p vcf"
cmd+=" $output_vcf_gz"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
