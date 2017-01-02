#!/bin/bash
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

#define default values
COL_CONFIG_DEFAULT="ALL"
VCF_REGION_DEFAULT=""

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-k {name}          dataset name (required)
-i {file}          specify tabix file (required)
-r {region}        specify vcf region to be exported (default:None)
-o {file}          output file name (required)
EOF
)

while getopts ":k:i:r:c:o:" OPTION; do
  case "$OPTION" in
    k)
      dataset_name="$OPTARG"
      ;;
    i)
      tabix_file="$OPTARG"
      ;;
    r)
      vcf_region="$OPTARG"
      ;;
    o)
      out_file="$OPTARG"
      ;;
    *)
      die "unrecognized option from executing: $0 $@"
      ;;
  esac
done

[ ! -z $dataset_name ] || die "Please specify dataset name (-k)"
[ ! -z $tabix_file ] || die "Please specify tabix file (-i)"
[ ! -z $out_file ] || die "Plesae specify output file name (-o)"
[ -f $tabix_file ] || die "$tabix_file is not a valid file name"

#setting default values:
: ${vcf_region=$VCF_REGION_DEFAULT}

## ****************************************  display configuration  ****************************************
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "parameters"
info_msg "  $params"
info_msg
info_msg "description"
info_msg "  A script to count/calculate mutation statistics. There are three kind of frequencies calculated:"
info_msg "    - genotyping frequency: It's the ratio of \"number of samples being genotyped\"/\"total number samples\""
info_msg "    - allelic frequency: It's the ratio of \"number of that particular allele in the samples\"/(\"number of genotyped samples\"*2)"
info_msg "    - population frequency: It's the ratio of \"number of that particular allele in the samples\"/(\"total number of samples\"*2)"
info_msg
info_msg
## display required configuration
info_msg "overall configuration"
display_param "dataset name (-k)" "$dataset_name"
display_param "tabix file (-i)" "$tabix_file"
display_param "statistics output file (-o)" "$out_file"

## display optional configuration
info_msg
info_msg "optional configuration"
if [ ! -z "$vcf_region" ]; then
    display_param "vcf region (-r)" "$vcf_region"
    IFS=$',' read -ra vcf_region_list <<< "$vcf_region"
    if [ $((${#vcf_region_list[@]})) -gt 1 ]; then
    for (( i=0; i<$((${#vcf_region_list[@]})); i++ ))
    do
        display_param "      region $(( i+1 ))" "${vcf_region_list[$i]}"
    done
    fi
else
    display_param "vcf region" "ALL"
fi

# ****************************************  executing  ****************************************
VCF_QUERY_FORMAT="'%CHROM\t%POS\t%REF\t%ALT\t%INFO/Het_KGA\t%INFO/Hom_KGA\t%INFO/Hemi_KGA\t%INFO/AC_KGA\t%INFO/AN_KGA\n'"
TOTAL_SAMPLES=1000
TOTAL_ALLELES=`echo $TOTAL_SAMPLES*2 | bc -l`
IDX_CHR_COL=1
IDX_POS_COL=2
IDX_REF_COL=3
IDX_ALT_COL=4
IDX_HET_COL=5
IDX_HOM_COL=6
IDX_HEMI_COL=7
IDX_AC_COL=8
IDX_AN_COL=9

function query_vcf {
    query_region=$1
    
    query_cmd="tabix"
    query_cmd+=" -h"
    query_cmd+=" $tabix_file"
    if [ ! -z "$query_region" ]; then
        query_cmd+=" $query_region"
    fi
    query_cmd+=" | vt decompose -s -"
    query_cmd+=" | vcf-query -f $VCF_QUERY_FORMAT"
    query_cmd+=" | grep -v \"*\""
    query_cmd+=" | awk -F'\t' '{printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\n\", \$$IDX_CHR_COL, \$$IDX_POS_COL, \$$IDX_POS_COL, \$$IDX_REF_COL, \$$IDX_ALT_COL, \$$IDX_HET_COL, \$$IDX_HOM_COL/2, \$$IDX_HEMI_COL, \$$IDX_AC_COL, \$$IDX_AN_COL, \$$IDX_AN_COL/$TOTAL_ALLELES, \$$IDX_AC_COL/\$$IDX_AN_COL, \$$IDX_AC_COL/$TOTAL_ALLELES}'"
    info_msg
    info_msg "generating vcf genotyping using data from $query_cmd"
    eval "$query_cmd" 
}

function count_frequency {
    region="$1"

    # calculate statistics
    query_vcf "$1" 
}

col_prefix=$( echo $dataset_name | awk '{print toupper($0)}' )
# create header
header="#Chr"
header+="\tStart"
header+="\tEnd"
header+="\tRef"
header+="\tAlt"

header+="\t$col_prefix"_HET
header+="\t$col_prefix"_HOM
header+="\t$col_prefix"_HEMI
header+="\t$col_prefix"_AC
header+="\t$col_prefix"_GT
header+="\t$col_prefix"_GF
header+="\t$col_prefix"_AF
header+="\t$col_prefix"_PF
echo -e "$header" > "$out_file"
        
if [ ! -z "$vcf_region" ]; then
    for (( n=0; n<$((${#vcf_region_list[@]})); n++ ))
    do
        count_frequency "${vcf_region_list[$n]}" >> "$out_file"
    done
else
    count_frequency "" >> "$out_file"
fi

new_section_txt "F I N I S H <$script_name>"
