#!/bin/bash 
set -e
set -u
set -o pipefail

source $PYCMM/bash/cmm_functions.sh

module load samtools

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}           input VCF file (required)
-o {file}           output file (required)
-h                  this help
EOF
)

# parse option
while getopts ":i:o:h" OPTION; do
  case "$OPTION" in
    i)
      input_vcf="$OPTARG"
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

[ ! -z $out_file ] || die "output titv count is required (-o)"
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
display_param "output file (-o)" "$out_file"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************

new_section_txt "Decompose and leftalign"

tmp_left_align="$working_dir/tmp_la"
cmd="gunzip -c"
cmd+=" $input_vcf"
cmd+=" | vt decompose -s -"
cmd+=" | vt normalize -q -r $GRCH37_REF -"
cmd+=" | grep -Pv \"\t\*\t\""
cmd+=" | grep -P \"\tPASS\t\""
cmd+=" > $tmp_left_align"
eval_cmd "$cmd"

new_section_txt "Removing non-titv variants"

tmp_removed_non_titv="$working_dir/tmp_removed_non_titv"
cmd="awk -F\$'\t'"
cmd+=" '{ if (length(\$4) == 1 && length(\$5) == 1) print \$0 }'"
#cmd+=" '{ if ((\$4 == \"C\" || \$4 == \"T\") && (length(\$5) == 1)) print \$0 }'"
cmd+=" $tmp_left_align"
cmd+=" > $tmp_removed_non_titv"
eval_cmd "$cmd"

new_section_txt "generate list of all titv events in all samples"

info_msg
info_msg ">>> extracting zygosities <<<"

tmp_removed_non_titv_vcf="$working_dir/tmp_removed_non_titv.vcf"
cmd="tabix -h $input_vcf 1:1-1 > $tmp_removed_non_titv_vcf" 
eval_cmd "$cmd"
cmd="cat $tmp_removed_non_titv >> $tmp_removed_non_titv_vcf" 
eval_cmd "$cmd"

cmd="bgzip -f $tmp_removed_non_titv_vcf"
eval_cmd "$cmd"

cmd="tabix -p vcf $tmp_removed_non_titv_vcf.gz"
eval_cmd "$cmd"

tmp_gtz="$working_dir/tmp_gtz"
cmd="vcf-query"
cmd+=" -f '[%GTR\t]\n'"
cmd+=" $tmp_removed_non_titv_vcf.gz"
cmd+=" > $tmp_gtz"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> extracting codons <<<"

tmp_get_codon_cmds="$working_dir/tmp_get_codon_cmds"
cmd="awk -F\$'\t'"
cmd+=" '{ printf \"samtools faidx \$GRCH37_REF %s:%s-%s | tail -1\n\", \$1, \$2-1, \$2+1 }'"
cmd+=" $tmp_removed_non_titv"
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
cmd+=" '{ printf \"%s\t%s\t%s\t%s\n\", \$1, \$2, \$4, \$5 }'"
cmd+=" $tmp_removed_non_titv"
cmd+=" > $tmp_coors"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> extracting genes <<<"

tmp_genes="$working_dir/tmp_genes"
cmd="grep -o 'Gene.refGene=[a-zA-Z0-9]*;'"
cmd+=" $tmp_removed_non_titv"
cmd+=" | grep -oP '(?<=Gene.refGene=)[a-zA-Z0-9]*'"
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
cmd+=" -o 1.1,2.1"
cmd+=" <( sort $tmp_genes )"
cmd+=" <( cut -f4,13 $ANNOVAR_HUMAN_DB_DIR/hg19_refGene.txt | sort -k2,2 | uniq )"
cmd+=" > $tmp_strand_genes"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> merging strands, genes, coordinates, codons, and zygosities together <<<"

tmp_strand_genes_coors_codons_gtz="$working_dir/tmp_strand_genes_coors_codons_gtz"
cmd="paste"
cmd+=" <( cut -f2 $tmp_strand_genes )"
cmd+=" <( sort -k1,1 $tmp_genes_coors_codons_gtz )"
cmd+=" > $tmp_strand_genes_coors_codons_gtz"
eval_cmd "$cmd"

cmd="vcf-query -l"
cmd+=" $tmp_removed_non_titv_vcf.gz"
cmd+=" | tr \"\n\" \"\t\""
samples_list=`eval $cmd`
tmp_info_ready="$working_dir/tmp_info_ready"
echo -e "#Strand\tGene\tChr\tPos\tRef\tAlt\tcodon\t$samples_list" > $tmp_info_ready
cat $tmp_strand_genes_coors_codons_gtz >> $tmp_info_ready

new_section_txt "extracting Ti/Tv signatures"

first_sample_idx=8
n_samples=`echo $samples_list | wc -w`
#n_samples=1

function count_subevents {
    events=$1
    strand=$2

    event_out=""
    for col_idx in $(seq $first_sample_idx $((n_samples+first_sample_idx-1)))
    do
        cmd="cut"
        cmd+=" -f$col_idx"
        cmd+=" $events"
        cmd+=" | awk -F\$'\t'"
        cmd+=" '{ if (\$0 == \"0/1\" || \$0 == \"1/1\") print \$0 }'"
        cmd+=" | wc -l"
        event_count=`eval $cmd`
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
#    info_msg "Counting transcribed event for $ts_fwd_ref>$ts_fwd_alt with codon $ts_fwd_codon in '+' gene and $ts_rev_ref>$ts_rev_alt with codon $ts_rev_codon in '-' gene"

    cmd="awk -F\$'\t'"
    cmd+=" '{ if ((\$1 == \"+\" && \$5 == \"$ts_fwd_ref\" && \$6 == \"$ts_fwd_alt\" && \$7 == \"$ts_fwd_codon\") || (\$1 == \"-\" && \$5 == \"$ts_rev_ref\" && \$6 == \"$ts_rev_alt\" && \$7 == \"$ts_rev_codon\")) print \$0 }'"
    cmd+=" $tmp_info_ready"
    cmd+=" > $tmp_events"
    eval_cmd "$cmd"

    count_subevents $tmp_events "transcribed"

    uts_fwd_ref=$rev_ref
    uts_fwd_alt=$rev_alt
    uts_fwd_codon=$rev_cpl_codon
    uts_rev_ref=$ref
    uts_rev_alt=$alt
    uts_rev_codon=$codon
#    info_msg
#    info_msg "Counting untranscribed event for $uts_fwd_ref>$uts_fwd_alt with codon $uts_fwd_codon in '+' gene and $uts_rev_ref>$uts_rev_alt with codon $uts_rev_codon in '-' gene"

    cmd="awk -F\$'\t'"
    cmd+=" '{ if ((\$1 == \"+\" && \$5 == \"$uts_fwd_ref\" && \$6 == \"$uts_fwd_alt\" && \$7 == \"$uts_fwd_codon\") || (\$1 == \"-\" && \$5 == \"$uts_rev_ref\" && \$6 == \"$uts_rev_alt\" && \$7 == \"$uts_rev_codon\")) print \$0 }'"
    cmd+=" $tmp_info_ready"
    cmd+=" > $tmp_events"
    eval_cmd "$cmd"

    count_subevents $tmp_events "untranscribed"
}

#echo -e "sample_id\tevent\tif_transcribed" > $out_file
out_header="sample_id\tref\talt\tif_transcribed"
#for prime3 in A C G T
#do
#    for prime5 in A C G T
#    do
#        out_header+="\t$prime3"_"$prime5"
#    done
#done
#echo -e "$out_header" > $out_file
#echo -e "#Ref\tAlt\tcodon\tif_transcribed\t$samples_list" > $out_file
tmp_reshape_siguatures_info="$working_dir/tmp_reshape_siguatures_info"
paste_buffer="$working_dir/paste_buffer"
:>$paste_buffer
tmp_paste_events_count="$working_dir/tmp_paste_events_count"
:>$tmp_paste_events_count
array_samples_list=( $samples_list )  

# create "reshape" signatures info
for ref in C T
do
#    for alt in A C
    for alt in A C G T
    do
        if [ "$ref" != "$alt" ]
        then
            for col_idx in $(seq 0 $((n_samples+-1)))
            do
                echo -e "${array_samples_list[$col_idx]}\t$ref\t$alt\tTranscribed" >> $tmp_reshape_siguatures_info
#                echo -e "${array_samples_list[$col_idx]}\t$ref\t$alt\tUntranscribed" >> $tmp_reshape_siguatures_info
            done
            for col_idx in $(seq 0 $((n_samples+-1)))
            do
                echo -e "${array_samples_list[$col_idx]}\t$ref\t$alt\tUnTranscribed" >> $tmp_reshape_siguatures_info
#                echo -e "${array_samples_list[$col_idx]}\t$ref\t$alt\tUntranscribed" >> $tmp_reshape_siguatures_info
            done

        fi
    done
done
#count_events T C ATA
#for ref in T
tmp_events_count="$working_dir/tmp_events_count"
#for prime3 in A
for prime3 in A C G T
do
#    for prime5 in A
    for prime5 in A C G T
    do
        out_header+="\t$prime3"_"$prime5"
        # each column in event counting represent a pair of 3' and 5'
        :>$tmp_events_count
        for ref in C T
        do
#            for alt in A C
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
##count_events T C ATA
##for ref in T
#for ref in C T
#do
#    for alt in A G
##    for alt in A C G T
#    do
#        if [ "$ref" != "$alt" ]
#        then
#            # add "reshape" signatures info
#            for col_idx in $(seq 1 $((n_samples)))
#            do
#                echo -e "${array_samples_list[$col_idx]}\t$ref\t$alt\tTranscribed" >> $tmp_reshape_siguatures_info
##                echo -e "${array_samples_list[$col_idx]}\t$ref\t$alt\tUntranscribed" >> $tmp_reshape_siguatures_info
#            done
#
#            for prime3 in A C
##            for prime3 in A C G T
#            do
#                for prime5 in A C
##                for prime5 in A C G T
#                do
#                    codon="$prime3$ref$prime5"
#                    out_header+="\t$codon"
#                    count_events $ref $alt $codon
#                    paste "$paste_buffer" "$tmp_events_count" > "$tmp_paste_events_count"
#                    cp "$tmp_paste_events_count" "$paste_buffer"
#                done
#            done
#        fi
#    done
#done
echo -e "$out_header" > $out_file

paste "$tmp_reshape_siguatures_info" "$tmp_paste_events_count" >> $out_file

new_section_txt "F I N I S H <$script_name>"

