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
[ ! -z $suggested_gene_list ] || die "please indicate gene list (-g)"
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
display_param "reference gene database (-d)" "$reference_gene"
display_param "gene list (-g)" "$suggested_gene_list"
display_param "output file (-o)" "$out_file"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************
 
new_section_txt "Extract list of coding regions"

tmp_refgene="$working_dir/refgene"
cmd="sort -k3,3 -k5,5n $reference_gene"
cmd+=" | cut -f3,5,6,13"
cmd+=" > $tmp_refgene"
eval_cmd "$cmd"

tmp_gene_list="$working_dir/gene_list"
# look for gene that only present in one place (not more than 1 chromosome)
cmd="cut -f1,4 $tmp_refgene"
cmd+=" | sort"
cmd+=" | uniq"
cmd+=" | sort -k2,2"
cmd+=" | uniq -u -f1"
cmd+=" | cut -f2"
#cmd+=" | grep -v \"^MIR6\""
#cmd+=" | grep -v \"^MIR4\""
#cmd+=" | grep -v \"^MIR3\""
#cmd+=" | grep -v \"^SNAR-A\""
##cmd+=" | grep -v \"^FAM74\""
#cmd+=" | grep -vP \"^LINC[0-9]\""
#cmd+=" | grep -vP \"^L0C[0-9]\""
#cmd+=" | grep -vP \"^NBPF[0-9]\""
#cmd+=" | grep -v \"^ULK4P3$\""
#cmd+=" | grep -v \"^TTTY17A$\""
#cmd+=" | grep -v \"^TTTY17B$\""
#cmd+=" | grep -v \"^TTTY17C$\""
#cmd+=" | grep -v \"^LOC103908605$\""
#cmd+=" | grep -v \"^SPDYE8P$\""
#cmd+=" | grep -v \"^PMS2P5$\""
#cmd+=" | grep -v \"^PMS2P7$\""
#cmd+=" | grep -v \"^STAG3L1$\""
#cmd+=" | grep -v \"^RGPD5$\""
#cmd+=" | grep -v \"^TTTY8$\""
#cmd+=" | grep -v \"^TTTY8B$\""
#cmd+=" | grep -v \"^TTTY7$\""
#cmd+=" | grep -v \"^TTTY7B$\""
#cmd+=" | grep -v \"^TTTY21$\""
#cmd+=" | grep -v \"^TTTY21B$\""
#cmd+=" | grep -v \"^TTTY2$\""
#cmd+=" | grep -v \"^TTTY2B$\""
#cmd+=" | grep -v \"^TTTY1$\""
#cmd+=" | grep -v \"^TTTY1B$\""
#cmd+=" | grep -v \"^LOC100288162$\""
#cmd+=" | grep -v \"^TTTY23$\""
#cmd+=" | grep -v \"^TTTY23B$\""
#cmd+=" | grep -v \"^RNU1-13P$\""
#cmd+=" | grep -v \"^LOC100132057$\""
#cmd+=" | grep -v \"^LOC554249$\""
#cmd+=" | grep -v \"^LOC102724238$\""
#cmd+=" | grep -v \"^PPIAL4B$\""
cmd+=" | grep -v \"^NBPF20$\""
cmd+=" | grep -v \"^RGPD5$\""
cmd+=" | grep -v \"^NBPF9$\""
cmd+=" | grep -v \"^NBPF8$\""
cmd+=" | grep -v \"^PPIAL4C$\""
#cmd+=" | grep -v \"^LOC102723709$\""
#cmd+=" | grep -v \"^FAM74A6$\""
#cmd+=" | grep -v \"^LOC286297$\""
cmd+=" | grep -v \"^SPATA31A7$\""
cmd+=" | grep -v \"^SPATA31A5$\""
cmd+=" | grep -v \"^ANKRD20A3$\""
#cmd+=" | grep -v \"^FAM74A1$\""
#cmd+=" | grep -v \"^LOC101928381$\""
cmd+=" | grep -v \"^FOXD4L4$\""
cmd+=" | grep -v \"^FAM72C$\""
#cmd+=" | grep -v \"^SRGAP2D$\""
cmd+=" > $tmp_gene_list"
eval_cmd "$cmd"

tmp_filtered_gene_list="$working_dir/filtered_gene_list"
cmd="join"
cmd+=" -1 1"
cmd+=" -2 1"
cmd+=" -t $'\t'"
cmd+=" -o 0"
cmd+=" $tmp_gene_list"
cmd+=" <( sort $suggested_gene_list )"
cmd+=" > $tmp_filtered_gene_list"
eval_cmd "$cmd"

tmp_coding_region="$working_dir/coding_region"
tmp_debug="$working_dir/debug"
:>$tmp_coding_region
:>$tmp_debug
for gene_name in `cat $tmp_filtered_gene_list`
do
#    info_msg "Processing information of gene: $gene_name"
    query_cmd="grep -P '\t$gene_name$' $tmp_refgene"
    query_cmd+=" | cut -f1"
    query_cmd+=" | head -1"
    tmp=`eval "$query_cmd"`
    chromosome="${tmp/chr/}"

    # get start position (the minimum one)
    query_cmd="grep -P '\t$gene_name$' $tmp_refgene"
    query_cmd+=" | cut -f2"
    query_cmd+=" | sort -n"
    query_cmd+=" | head -1"
    start_pos=`eval "$query_cmd"`

    # get end position (the maximum one)
    query_cmd="grep -P '\t$gene_name$' $tmp_refgene"
    query_cmd+=" | cut -f3"
    query_cmd+=" | sort -n"
    query_cmd+=" | tail -1"
    end_pos=`eval "$query_cmd"`

    echo "$chromosome:$start_pos-$end_pos" >> $tmp_coding_region
    echo -e "$chromosome:$start_pos-$end_pos\t$(expr $end_pos - $start_pos)\t$gene_name" >> $tmp_debug
done 

new_section_txt "Extract variants from coding regions that pass the QC"

tmp_exome="$working_dir/exome"
:> $tmp_exome
for region in `cat $tmp_coding_region`
do
    cmd="tabix $input_vcf $region"
    cmd+=" | grep -P \"\tPASS\t\""
    cmd+=" >> $tmp_exome"
    eval_cmd "$cmd"
done

tmp_uniq_exome="$working_dir/uniq_exome.vcf"
tmp_pass_qc="$working_dir/pass_qc.vcf"
tabix -h $input_vcf 1:1-1 > $tmp_uniq_exome
cmd="sort $tmp_exome"
cmd+=" | uniq"
cmd+=" >> $tmp_uniq_exome"
eval_cmd "$cmd"
cmd="vcf-sort -c $tmp_uniq_exome > $tmp_pass_qc"
eval_cmd "$cmd"
cmd="bgzip -f $tmp_pass_qc"
eval_cmd "$cmd"
cmd="tabix -f -p vcf $tmp_pass_qc.gz"
eval_cmd "$cmd"

new_section_txt "Extract only coding and splice regions"

tmp_ta="$working_dir/ta"
table_annovar_cmd="table_annovar.pl"
table_annovar_cmd+=" $tmp_pass_qc.gz"
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
cmd+=" >> $tmp_output"
eval_cmd "$cmd"

cmd="grep -v \"^#\" $tmp_ta.hg19_multianno.vcf"
cmd+=" | grep \"=splicing\""
cmd+=" >> $tmp_output"
eval_cmd "$cmd"

cmd="vcf-sort -c $tmp_output > $out_file"
eval_cmd "$cmd"

cmd="bgzip -f $out_file"
eval_cmd "$cmd"

cmd="tabix -p vcf $out_file.gz"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
