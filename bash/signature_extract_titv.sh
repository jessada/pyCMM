#!/bin/bash 
set -e
set -u
set -o pipefail

source $PYCMM/bash/cmm_functions.sh

module load samtools
module load annovar

#define default values
GT_FORMAT_DEFAULT="GTR"

gt_format=$GT_FORMAT_DEFAULT
frequency_ratio="1"

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}           input VCF file (required)
-g {file}           genotype format (default=$GT_FORMAT_DEFAULT)
-r {percent}        filtering frequency ratio (detault=1)
-e {file}           output events file (required)
-o {file}           output titv file (required)
-h                  this help
EOF
)

# parse option
while getopts ":i:g:e:r:o:h" OPTION; do
  case "$OPTION" in
    i)
      input_vcf="$OPTARG"
      ;;
    g)
      gt_format="$OPTARG"
      ;;
    r)
      frequency_ratio="$OPTARG"
      ;;
    e)
      out_event_file="$OPTARG"
      ;;
    o)
      out_titv_file="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $out_event_file ] || die "output events file is required (-a)"
[ ! -z $out_titv_file ] || die "output titv count is required (-o)"
[ -f "$input_vcf" ] || die "$input_vcf is not found"


working_dir=`mktemp -d`

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

time_stamp=$( date )

working_dir=`mktemp -d`

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  This application will count the following changes in all possible 3'5'"
info_msg "    - C > A"
info_msg "    - C > G"
info_msg "    - C > T"
info_msg "    - T > A"
info_msg "    - T > C"
info_msg "    - T > G"
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
display_param "genotype format (-g)" "$gt_format"
display_param "filtering frequency ratiio (-r)" "$frequency_ratio"
display_param "output event file (-e)" "$out_event_file"
display_param "output titv file (-o)" "$out_titv_file"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************

new_section_txt "Decompose variants"

tmp_decompose="$working_dir/tmp_decompose"
cmd="gunzip -c"
cmd+=" $input_vcf"
cmd+=" | vt decompose -s -"
cmd+=" | grep -Pv \"\t\*\t\""
cmd+=" | grep -v \"\\x3b\""
cmd+=" > $tmp_decompose"
eval_cmd "$cmd"

# suppress QUAL and INFO column
printf_phrase="%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s"
param_phrase="\$1, \$2, \$3, \$4, \$5, \$7, \$9"

n_samples=`vcf-query -l $input_vcf | wc -l `
for ((n_sample=1; n_sample<=$n_samples; n_sample++));
do
    printf_phrase+="\t%s"
    param_phrase+=", \$$((n_sample+9))"
done;

new_section_txt "Removing non-titv variants and empty QUAL and INFO columns"

tmp_removed_non_titv="$working_dir/tmp_removed_non_titv"
cmd="awk -F\$'\t'"
cmd+=" '{ if (length(\$4) == 1 && length(\$5) == 1) print \$0 }'"
cmd+=" $tmp_decompose"
cmd+=" | awk -F '\t' '{ printf \"$printf_phrase\n\", $param_phrase}'"
cmd+=" > $tmp_removed_non_titv"
eval_cmd "$cmd"

tmp_removed_non_titv_vcf="$working_dir/tmp_removed_non_titv.vcf"
cmd="tabix -h $input_vcf 1:1-1 > $tmp_removed_non_titv_vcf" 
eval_cmd "$cmd"
cmd="cat $tmp_removed_non_titv >> $tmp_removed_non_titv_vcf" 
eval_cmd "$cmd"

cmd="bgzip -f $tmp_removed_non_titv_vcf"
eval_cmd "$cmd"

cmd="tabix -p vcf $tmp_removed_non_titv_vcf.gz"
eval_cmd "$cmd"

new_section_txt "Annotate variants with refGene information"

tmp_ta="$working_dir/ta"
table_annovar_cmd="table_annovar.pl"
table_annovar_cmd+=" $tmp_removed_non_titv_vcf.gz"
table_annovar_cmd+=" $ANNOVAR_HUMAN_DB_DIR"
table_annovar_cmd+=" -buildver hg19"
table_annovar_cmd+=" -out $tmp_ta"
table_annovar_cmd+=" -remove"
table_annovar_cmd+=" -protocol refGene,gnomad_genome"
table_annovar_cmd+=" -operation g,f"
table_annovar_cmd+=" -nastring ."
table_annovar_cmd+=" -vcfinput"
eval_cmd "$table_annovar_cmd"

multianno_vcf="$tmp_ta.hg19_multianno.vcf"
tmp_ta_vcf="$tmp_ta.vcf"
cmd="sed 's/Func.refGene/Func_refGene/g'"
cmd+=" $multianno_vcf"
cmd+=" | sed 's/Gene.refGene/Gene_refGene/g'"
cmd+=" > $tmp_ta_vcf"
eval_cmd "$cmd"

cmd="bgzip -f $tmp_ta_vcf"
eval_cmd "$cmd"

cmd="tabix -p vcf $tmp_ta_vcf.gz"
eval_cmd "$cmd"

new_section_txt "generate list of all titv events in all samples"

info_msg
info_msg ">>> extracting variants information <<<"

TMP_VARIANTS_LIST_CHROM_COL_IDX=1
TMP_VARIANTS_LIST_POS_COL_IDX=2
TMP_VARIANTS_LIST_REF_COL_IDX=3
TMP_VARIANTS_LIST_ALT_COL_IDX=4
TMP_VARIANTS_LIST_GENE_REFGENE_COL_IDX=5
TMP_VARIANTS_LIST_GNOMAD_GENOME_ALL_COL_IDX=6
TMP_VARIANTS_LIST_GNOMAD_GENOME_AFR_COL_IDX=7
TMP_VARIANTS_LIST_GNOMAD_GENOME_AMR_COL_IDX=8
TMP_VARIANTS_LIST_GNOMAD_GENOME_ASJ_COL_IDX=9
TMP_VARIANTS_LIST_GNOMAD_GENOME_EAS_COL_IDX=10
TMP_VARIANTS_LIST_GNOMAD_GENOME_FIN_COL_IDX=11
TMP_VARIANTS_LIST_GNOMAD_GENOME_NFE_COL_IDX=12
TMP_VARIANTS_LIST_GNOMAD_GENOME_OTH_COL_IDX=13
TMP_VARIANTS_LIST_FIRST_GT_COL_IDX=14

vcf_query_format="'"
vcf_query_format+="%CHROM"
vcf_query_format+="\t%POS"
vcf_query_format+="\t%REF"
vcf_query_format+="\t%ALT"
vcf_query_format+="\t%INFO/Gene_refGene"
vcf_query_format+="\t%INFO/gnomAD_genome_ALL"
vcf_query_format+="\t%INFO/gnomAD_genome_AFR"
vcf_query_format+="\t%INFO/gnomAD_genome_AMR"
vcf_query_format+="\t%INFO/gnomAD_genome_ASJ"
vcf_query_format+="\t%INFO/gnomAD_genome_EAS"
vcf_query_format+="\t%INFO/gnomAD_genome_FIN"
vcf_query_format+="\t%INFO/gnomAD_genome_NFE"
vcf_query_format+="\t%INFO/gnomAD_genome_OTH"
vcf_query_format+="[\t%$gt_format]"
vcf_query_format+="\n"
vcf_query_format+="'"

awk_filtering_condition="\$$TMP_VARIANTS_LIST_GNOMAD_GENOME_ALL_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_GENOME_AFR_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_GENOME_AMR_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_GENOME_ASJ_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_GENOME_EAS_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_GENOME_FIN_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_GENOME_NFE_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_GENOME_OTH_COL_IDX < $frequency_ratio"

tmp_filtered_variants_list="$working_dir/tmp_filtered_variants_list"
cmd="vcf-query"
cmd+=" -f $vcf_query_format"
cmd+=" $tmp_ta_vcf.gz"
cmd+=" | awk '{ if ($awk_filtering_condition) print \$0 }'"
cmd+=" > $tmp_filtered_variants_list"
eval_cmd "$cmd"

tmp_gtz="$working_dir/tmp_gtz"
info_msg
info_msg
info_msg ">>> extracting zygosities <<<"
raw_cmd=" cut"
raw_cmd+=" -f"
for col_idx in $(seq $TMP_VARIANTS_LIST_FIRST_GT_COL_IDX $((TMP_VARIANTS_LIST_FIRST_GT_COL_IDX+n_samples-1)))
do
    raw_cmd+=",$col_idx"
done
raw_cmd+=" $tmp_filtered_variants_list"
raw_cmd+=" > $tmp_gtz"
cmd="$( echo $raw_cmd | sed 's/-f,/-f/g' )"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> extracting codons <<<"

tmp_get_codon_cmds="$working_dir/tmp_get_codon_cmds"
cmd="awk -F\$'\t'"
cmd+=" '{ printf \"samtools faidx \$GRCH37_REF %s:%s-%s | tail -1\n\", \$$TMP_VARIANTS_LIST_CHROM_COL_IDX, \$$TMP_VARIANTS_LIST_POS_COL_IDX-1, \$$TMP_VARIANTS_LIST_POS_COL_IDX+1 }'"
cmd+=" $tmp_filtered_variants_list"
cmd+=" > $tmp_get_codon_cmds"
eval_cmd "$cmd"

tmp_codons="$working_dir/tmp_codons"
:>$tmp_codons
while read cmd; do
    eval "$cmd >> $tmp_codons"
done <  $tmp_get_codon_cmds

info_msg
info_msg
info_msg ">>> extracting coordinates <<<"

tmp_coors="$working_dir/tmp_coors"
cmd="awk -F\$'\t'"
cmd+=" '{ printf \"%s\t%s\t%s\t%s\n\", \$$TMP_VARIANTS_LIST_CHROM_COL_IDX, \$$TMP_VARIANTS_LIST_POS_COL_IDX, \$$TMP_VARIANTS_LIST_REF_COL_IDX, \$$TMP_VARIANTS_LIST_ALT_COL_IDX }'"
cmd+=" $tmp_filtered_variants_list"
cmd+=" > $tmp_coors"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> extracting genes <<<"

tmp_genes="$working_dir/tmp_genes"
cmd="awk -F\$'\t'"
cmd+=" '{ printf \"%s\n\", \$$TMP_VARIANTS_LIST_GENE_REFGENE_COL_IDX }'"
cmd+=" $tmp_filtered_variants_list"
cmd+=" > $tmp_genes"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> merging genes, coordinates, codons, and zygosities together <<<"

tmp_genes_coors_codons_gtz="$working_dir/tmp_genes_coors_codons_gtz"
cmd="paste"
cmd+=" $tmp_genes"
cmd+=" $tmp_coors"
cmd+=" $tmp_codons"
cmd+=" $tmp_gtz"
cmd+=" > $tmp_genes_coors_codons_gtz"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> extracting transcription strands <<<"

tmp_strand_genes="$working_dir/tmp_strand_genes"
cmd="join"
cmd+=" -1 1"
cmd+=" -2 2"
cmd+=" -t $'\t'"
cmd+=" -o 2.1,1.1"
cmd+=" <( sort $tmp_genes )"
cmd+=" <( cut -f4,13 $ANNOVAR_HUMAN_DB_DIR/hg19_refGene.txt | sort -k2,2 | uniq )"
cmd+=" | sort -k2,2"
cmd+=" | uniq -f1"
cmd+=" > $tmp_strand_genes"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> merging strands, genes, coordinates, codons, and zygosities together <<<"

cmd="vcf-query -l"
cmd+=" $input_vcf"
cmd+=" | tr \"\n\" \"\t\""
samples_list=`eval $cmd`
echo -e "#Strand\tGene\tChr\tPos\tRef\tAlt\tcodon\t$samples_list" > $out_event_file

out_join_phrase="1.1,0,2.2,2.3,2.4,2.5,2.6"

for ((n_sample=1; n_sample<=$n_samples; n_sample++));
do
    out_join_phrase+=",2.$((n_sample+6))"
done

cmd="join"
cmd+=" -1 2"
cmd+=" -2 1"
cmd+=" -t $'\t'"
cmd+=" -o $out_join_phrase"
cmd+=" <( sort -k2,2 $tmp_strand_genes )"
cmd+=" <( sort -k1,1 $tmp_genes_coors_codons_gtz )"
cmd+=" >> $out_event_file"
eval_cmd "$cmd"

new_section_txt "extracting Ti/Tv signatures"

first_sample_idx=8
n_samples=`echo $samples_list | wc -w`

function count_subevents {
    events=$1
    strand=$2

    event_out=""
    for col_idx in $(seq $first_sample_idx $((n_samples+first_sample_idx-1)))
    do
        cmd="cut"
        cmd+=" -f$col_idx"
        cmd+=" $events"
        cmd+=" | grep -o 1"
        cmd+=" | wc -l"
        event_count=$( eval $cmd  || true ) 
        event_out+="\n$event_count"
    done
    echo -e "$event_out" | grep -v "^$" >> $tmp_events_count
}

function count_events {
    ref=$1
    alt=$2
    codon=$3

    case "$ref" in
      C)
        rev_ref="G"
        ;;
      T)
        rev_ref="A"
        ;;
      *)
        die "unrecognized ref: $ref"
        ;;
    esac

    case "$alt" in
      A)
        rev_alt="T"
        ;;
      C)
        rev_alt="G"
        ;;
      G)
        rev_alt="C"
        ;;
      T)
        rev_alt="A"
        ;;
      *)
        die "unrecognized alt: $alt"
        ;;
    esac

    case "$codon" in
      ACA)
        rev_cpl_codon="TGT"
        ;;
      ACC)
        rev_cpl_codon="GGT"
        ;;
      ACG)
        rev_cpl_codon="CGT"
        ;;
      ACT)
        rev_cpl_codon="AGT"
        ;;
      CCA)
        rev_cpl_codon="TGG"
        ;;
      CCC)
        rev_cpl_codon="GGG"
        ;;
      CCG)
        rev_cpl_codon="CGG"
        ;;
      CCT)
        rev_cpl_codon="AGG"
        ;;
      GCA)
        rev_cpl_codon="TGC"
        ;;
      GCC)
        rev_cpl_codon="GGC"
        ;;
      GCG)
        rev_cpl_codon="CGC"
        ;;
      GCT)
        rev_cpl_codon="AGC"
        ;;
      TCA)
        rev_cpl_codon="TGA"
        ;;
      TCC)
        rev_cpl_codon="GGA"
        ;;
      TCG)
        rev_cpl_codon="CGA"
        ;;
      TCT)
        rev_cpl_codon="AGA"
        ;;
      ATA)
        rev_cpl_codon="TAT"
        ;;
      ATC)
        rev_cpl_codon="GAT"
        ;;
      ATG)
        rev_cpl_codon="CAT"
        ;;
      ATT)
        rev_cpl_codon="AAT"
        ;;
      CTA)
        rev_cpl_codon="TAG"
        ;;
      CTC)
        rev_cpl_codon="GAG"
        ;;
      CTG)
        rev_cpl_codon="CAG"
        ;;
      CTT)
        rev_cpl_codon="AAG"
        ;;
      GTA)
        rev_cpl_codon="TAC"
        ;;
      GTC)
        rev_cpl_codon="GAC"
        ;;
      GTG)
        rev_cpl_codon="CAC"
        ;;
      GTT)
        rev_cpl_codon="AAC"
        ;;
      TTA)
        rev_cpl_codon="TAA"
        ;;
      TTC)
        rev_cpl_codon="GAA"
        ;;
      TTG)
        rev_cpl_codon="CAA"
        ;;
      TTT)
        rev_cpl_codon="AAA"
        ;;
      *)
        die "unrecognized codon: $codon"
        ;;
    esac

    tmp_events="$working_dir/tmp_events"

    ts_fwd_ref=$ref
    ts_fwd_alt=$alt
    ts_fwd_codon=$codon
    ts_rev_ref=$rev_ref
    ts_rev_alt=$rev_alt
    ts_rev_codon=$rev_cpl_codon

    cmd="awk -F\$'\t'"
    cmd+=" '{ if ((\$1 == \"+\" && \$5 == \"$ts_fwd_ref\" && \$6 == \"$ts_fwd_alt\" && \$7 == \"$ts_fwd_codon\") || (\$1 == \"-\" && \$5 == \"$ts_rev_ref\" && \$6 == \"$ts_rev_alt\" && \$7 == \"$ts_rev_codon\")) print \$0 }'"
    cmd+=" $out_event_file"
    cmd+=" > $tmp_events"
    eval_cmd "$cmd"

    count_subevents $tmp_events "transcribed"

    uts_fwd_ref=$rev_ref
    uts_fwd_alt=$rev_alt
    uts_fwd_codon=$rev_cpl_codon
    uts_rev_ref=$ref
    uts_rev_alt=$alt
    uts_rev_codon=$codon

    cmd="awk -F\$'\t'"
    cmd+=" '{ if ((\$1 == \"+\" && \$5 == \"$uts_fwd_ref\" && \$6 == \"$uts_fwd_alt\" && \$7 == \"$uts_fwd_codon\") || (\$1 == \"-\" && \$5 == \"$uts_rev_ref\" && \$6 == \"$uts_rev_alt\" && \$7 == \"$uts_rev_codon\")) print \$0 }'"
    cmd+=" $out_event_file"
    cmd+=" > $tmp_events"
    eval_cmd "$cmd"

    count_subevents $tmp_events "untranscribed"
}

out_header="sample_id\tref\talt\tif_transcribed"
tmp_reshape_siguatures_info="$working_dir/tmp_reshape_siguatures_info"
paste_buffer="$working_dir/paste_buffer"
:>$paste_buffer
tmp_paste_events_count="$working_dir/tmp_paste_events_count"
:>$tmp_paste_events_count
array_samples_list=( $samples_list )  

# create "reshape" signatures info
for ref in C T
do
    for alt in A C G T
    do
        if [ "$ref" != "$alt" ]
        then
            for col_idx in $(seq 0 $((n_samples-1)))
            do
                echo -e "${array_samples_list[$col_idx]}\t$ref\t$alt\tTranscribed" >> $tmp_reshape_siguatures_info
            done
            for col_idx in $(seq 0 $((n_samples-1)))
            do
                echo -e "${array_samples_list[$col_idx]}\t$ref\t$alt\tUnTranscribed" >> $tmp_reshape_siguatures_info
            done

        fi
    done
done

tmp_events_count="$working_dir/tmp_events_count"
for prime3 in A C G T
do
    for prime5 in A C G T
    do
        out_header+="\t$prime3"_"$prime5"
        # each column in event counting represent a pair of 3' and 5'
        :>$tmp_events_count
        for ref in C T
        do
            for alt in A C G T
            do
                if [ "$ref" != "$alt" ]
                then
                    codon="$prime3$ref$prime5"
                    count_events $ref $alt $codon
                fi
            done
        done
        paste "$paste_buffer" "$tmp_events_count" > "$tmp_paste_events_count"
        cp "$tmp_paste_events_count" "$paste_buffer"
    done
done
echo -e "$out_header" > $out_titv_file

paste -d '' "$tmp_reshape_siguatures_info" "$tmp_paste_events_count" >> $out_titv_file

new_section_txt "F I N I S H <$script_name>"

