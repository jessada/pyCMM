#!/bin/bash 

#set -e
#set -u
set -o pipefail

module load annovar
source $PYCMM/bash/cmm_functions.sh

REFERENCE_GENE_DEFAULT="$ANNOVAR_HUMAN_DB_DIR/hg19_refGene.txt"

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}         input VCF file (required)
-d {database}     specify reference gene database (default:$REFERENCE_GENE_DEFAULT)
-g {database}     gene list for searching (required)
-o {file}         output file (required)
-h                this help
EOF
)

# parse option
while getopts ":i:d:g:o:h" OPTION; do
  case "$OPTION" in
    i)
      input_vcf="$OPTARG"
      ;;
    d)
      reference_gene="$OPTARG"
      ;;
    g)
      suggested_gene_list="$OPTARG"
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

[ ! -z $input_vcf ] || die "an input VCF file is required (-i)"
#[ ! -z $suggested_gene_list ] || die "please indicate gene list (-g)"
[ ! -z $out_file ] || die "please indicate output file name (-o)"
[ -f "$input_vcf" ] || die "$input_vcf is not a valid file name."

#setting default values:
: ${reference_gene=$REFERENCE_GENE_DEFAULT}

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
#display_param "reference gene database (-d)" "$reference_gene"
#display_param "gene list (-g)" "$suggested_gene_list"
display_param "output file (-o)" "$out_file"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************
 
new_section_txt "Extract only coding and splice regions"

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

tmp_output="$working_dir/tmp_output.vcf"
picard_input="$working_dir/picard_input.vcf"

cmd="head -10000 $tmp_ta.hg19_multianno.vcf"
cmd+=" | grep \"^#\""
cmd+=" > $tmp_output"
eval_cmd "$cmd"

cmd="grep -v \"^#\" $tmp_ta.hg19_multianno.vcf"
cmd+=" | grep \"=exonic\""
cmd+=" | grep -P \"\tPASS\t\""
cmd+=" >> $tmp_output"
eval_cmd "$cmd"

cmd="grep -v \"^#\" $tmp_ta.hg19_multianno.vcf"
cmd+=" | grep \"=splicing\""
cmd+=" | grep -P \"\tPASS\t\""
cmd+=" >> $tmp_output"
eval_cmd "$cmd"

cmd="vcf-sort -c $tmp_output > $out_file"
eval_cmd "$cmd"

cmd="bgzip -f $out_file"
eval_cmd "$cmd"

cmd="tabix -p vcf $out_file.gz"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
