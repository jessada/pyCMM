#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-S {SAM file}       SAM file (required)
-o {file}           output sorted reads file (required)
-h                  this help
EOF
)

# parse option
while getopts ":S:o:h" OPTION; do
  case "$OPTION" in
    S)
      sam_file="$OPTARG"
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

[ ! -z $sam_file ] || die "a SAM input file is required (-S)"
[ ! -z $out_file ] || die "please indicate output file name (-o)"
[ -f "$sam_file" ] || die "$sam_file is not found"

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
info_msg "  Sort the SAM file and convert it to BAM"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM_DIR"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "SAM input file (-S)" "$sam_file"
display_param "output file (-o)" "$out_file"

# ****************************************  executing  ****************************************
cmd=" java -jar $CLASSPATH/SortSam.jar"
cmd+=" INPUT=$sam_file"
cmd+=" OUTPUT=$out_file"
cmd+=" SORT_ORDER=coordinate"
eval_cmd "$cmd"
new_section_txt "F I N I S H <$script_name>"

