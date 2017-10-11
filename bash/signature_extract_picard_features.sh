#!/bin/bash 
set -o pipefail

module load annovar
module load picard/2.0.1
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}         input VCF file (required)
-o {file}         output file prefix (required)
-h                this help
EOF
)

# parse option
while getopts ":i:o:h" OPTION; do
  case "$OPTION" in
    i)
      input_vcf="$OPTARG"
      ;;
    o)
      out_file_prefix="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $input_vcf ] || die "an input VCF file is required (-i)"
[ ! -z $out_file_prefix ] || die "please indicate output file prefix (-o)"
[ -f "$input_vcf" ] || die "$input_vcf is not a valid file name."

time_stamp=$( date )

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

working_dir=`mktemp -d`

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  decompose and left align vcf file"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "input VCF file (-i)" "$input_vcf"
display_param "output file prefix (-o)" "$out_file_prefix"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************
 
new_section_txt "Annovar annotation"

tmp_ta="$working_dir/ta"
table_annovar_cmd="table_annovar.pl"
table_annovar_cmd+=" $input_vcf"
table_annovar_cmd+=" $ANNOVAR_HUMAN_DB_DIR"
table_annovar_cmd+=" -buildver hg19"
table_annovar_cmd+=" -out $tmp_ta"
table_annovar_cmd+=" -remove"
table_annovar_cmd+=" -protocol refGene"
table_annovar_cmd+=" -operation g"
table_annovar_cmd+=" -nastring ."
table_annovar_cmd+=" -vcfinput"
eval_cmd "$table_annovar_cmd"


new_section_txt "Extract only coding and splice regions"

tmp_picard_input="$working_dir/tmp_picard_input.vcf"
picard_input="$working_dir/picard_input.vcf"

cmd="head -10000 $tmp_ta.hg19_multianno.vcf"
cmd+=" | grep \"^#\""
cmd+=" > $tmp_picard_input"
eval_cmd "$cmd"

cmd="grep -v \"^#\" $tmp_ta.hg19_multianno.vcf"
cmd+=" | grep \"=exonic\""
cmd+=" >> $tmp_picard_input"
eval_cmd "$cmd"

cmd="grep -v \"^#\" $tmp_ta.hg19_multianno.vcf"
cmd+=" | grep \"=splicing\""
cmd+=" >> $tmp_picard_input"
eval_cmd "$cmd"

cmd="vcf-sort $tmp_picard_input > $picard_input"
eval_cmd "$cmd"

cmd="bgzip -f $picard_input"
eval_cmd "$cmd"

cmd="tabix -p vcf $picard_input.gz"
eval_cmd "$cmd"


new_section_txt "Generate Picard features"

>&2 date
# compute variant calling metrics using picard
cmd="java -jar"
cmd+=" $CLASSPATH/picard.jar"
cmd+=" CollectVariantCallingMetrics"
cmd+=" INPUT=$picard_input.gz"
cmd+=" OUTPUT=$out_file_prefix"
cmd+=" DBSNP=$REFERENCE_DATA_DIR/dbsnp_138.b37.vcf"
eval_cmd "$cmd"
>&2 date

new_section_txt "F I N I S H <$script_name>"
