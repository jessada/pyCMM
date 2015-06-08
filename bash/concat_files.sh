#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-I {file}     input sam files to be concatenated, comma separated(required)
-o {file}     output concatenated file (required)
-h            this help
EOF
)

# parse option
while getopts ":I:o:h" OPTION; do
  case "$OPTION" in
    I)
      sam_files="$OPTARG"
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

[ ! -z $sam_files ] || die "sam files to be concatenated is required (-I)"
[ ! -z $out_file ] || die "output file name is required (-o)"

IFS=',' read -ra sam_files_list <<< "$sam_files"
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
info_msg "  SAM files concatenation"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "sam files (-I)"
for (( file_idx=0; file_idx<$((${#sam_files_list[@]})); file_idx++ ))
do
    display_param "  file #$(( file_idx+1 )) " "${sam_files_list[$file_idx]}"
done
display_param "output file (-o)" "$out_file"
# ****************************************  executing  ****************************************
> $out_file
for (( file_idx=0; file_idx<$((${#sam_files_list[@]})); file_idx++ ))
do
    cat "${sam_files_list[$file_idx]}" >> $out_file
done

new_section_txt "F I N I S H <$script_name>"

