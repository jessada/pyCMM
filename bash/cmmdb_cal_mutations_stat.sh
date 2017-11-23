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
-t {file}          specify vcf2avdb key table (required)
-r {region}        specify vcf region to be exported (default:None)
-c {patient list}  specify vcf columns to exported. This can be either in comma-separated format or it can be a file name (default:$COL_CONFIG_DEFAULT)
-o {file}          output file name (required)
EOF
)

while getopts ":k:i:t:r:c:o:" OPTION; do
  case "$OPTION" in
    k)
      dataset_name="$OPTARG"
      ;;
    i)
      tabix_file="$OPTARG"
      ;;
    t)
      vcf2avdb_key_table="$OPTARG"
      ;;
    r)
      vcf_region="$OPTARG"
      ;;
    c)
      col_config="$OPTARG"
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
[ ! -z $vcf2avdb_key_table ] || die "Please specify vcf2avdb key table (-t)"
[ ! -z $out_file ] || die "Plesae specify output file name (-o)"
[ -f $tabix_file ] || die "$tabix_file is not a valid file name"
[ -f $vcf2avdb_key_table ] || die "$vcf2avdb_key_table is not a valid file name"

#setting default values:
: ${vcf_region=$VCF_REGION_DEFAULT}
: ${col_config=$COL_CONFIG_DEFAULT}

if [ "$col_config" == "$COL_CONFIG_DEFAULT" ]
then
    col_count=$( vcf-query -l $tabix_file | wc -l)
    parsed_col_names=""
else
    if [ -f "$col_config" ]
    then
        parsed_col_names=`paste -sd, $col_config`
    else
        parsed_col_names="$col_config"
    fi

    IFS=',' read -ra col_list <<< "$parsed_col_names"
    for (( i=0; i<$((${#col_list[@]})); i++ ))
    do
        col_exist=$( vcf_col_exist $tabix_file ${col_list[$i]} )
	if [ "$col_exist" -ne 1 ]
	then
	    die "column ${col_list[$i]} is not exist"
	fi
    done
    col_count=${#col_list[@]}
fi

working_dir=`mktemp -d`
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
display_param "vcf2avdb key table (-t)" "$vcf2avdb_key_table"
display_param "statistics output file (-o)" "$out_file"

## display optional configuration
info_msg
info_msg "optional configuration"
if [ ! -z "$parsed_col_names" ]; then
    display_param "column names (-c)" "$parsed_col_names"
else
    display_param "column names" "ALL"
fi
display_param "column count" "$col_count"
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

## display misc configuration
info_msg
info_msg "misc configuration"
display_param "working direcotry" "$working_dir"

# ****************************************  executing  ****************************************

new_section_txt "E X E C U T I N G"

VCF_QUERY_FORMAT="'%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'"
COL_KEY_COUNT=4
IDX_0_CHR_COL=0
IDX_0_POS_COL=1
IDX_0_REF_COL=2
IDX_0_ALT_COL=3
IDX_0_GT_COL=4

function query_vcf {
    query_region=$1
    
    vcf_query_cmd="vcf-query "
    if [ ! -z "$query_region" ]; then
        vcf_query_cmd+=" -r $query_region"
    fi
    if [ ! -z "$parsed_col_names" ]; then
        vcf_query_cmd+=" -c $parsed_col_names"
    fi
    vcf_query_cmd+=" -f "$VCF_QUERY_FORMAT" $tabix_file "
    info_msg
    info_msg "generating vcf genotyping using data from $vcf_query_cmd"
    eval "$vcf_query_cmd" 
}

function count_frequency {
    region="$1"

    # calculate statistics
    query_vcf "$1" | 
    while read rec_in; do
        #parse input vcf record into vcf columns
        IFS=$'\t' read -ra rec_col <<< "$rec_in"
        chr=${rec_col[$IDX_0_CHR_COL]}
        pos=${rec_col[$IDX_0_POS_COL]}
        ref=${rec_col[$IDX_0_REF_COL]}
        alt_list=${rec_col[$IDX_0_ALT_COL]}

        # split ALT field in case that there are more than one alternate alleles
        # for all ALT
        IFS=',' read -ra alt <<< "$alt_list"
        for (( i=0; i<$((${#alt[@]})); i++ ))
        do
            rec_out=$( printf "%s\t%s\t%s\t%s\t%s" $chr $pos $pos $ref "${alt[$i]}" )
            # for all GT fields
            wt_count=0
            het_count=0
            hom_count=0
            oth_count=0
            na_count=0
            gt_count=0
            al_count=0
            for (( j=$IDX_0_GT_COL; j<$((${#rec_col[@]})); j++ ))
            do
                # count genotypes
            	if [ "${rec_col[$j]}" != "./." ] && [ "${rec_col[$j]}" != "." ]
                then
                    let gt_count++
                else
                    let na_count++
                fi
                # count alleles
                IFS='/' read -ra gt <<< "${rec_col[$j]}"
                # for both chromosomes
            	for (( k=0; k<$((${#gt[@]})); k++ ))
            	do
                    if [ "${gt[$k]}" = "${alt[$i]}" ]
        	        then
                        let al_count++
        	            if [ "${gt[0]}" = "${gt[1]}" ]
        	            then
                            let hom_count++
                            let al_count++
                            break
                        else
                            let het_count++
                        fi
                    elif [ "${gt[0]}" != "${alt[$i]}" ] && [ "${gt[1]}" != "${alt[$i]}" ] && [ "${gt[$k]}" != "$ref" ] && [ "${gt[0]}" != "." ] && [ "${gt[1]}" != "." ]
                    then
                        let oth_count++
                        break
                    fi
            	done
                if [ "${gt[0]}" = "${gt[1]}" ] && [ "${gt[0]}" = "$ref" ]
                then
                    let wt_count++
                fi
            done
            if [ $gt_count -eq 0 ]
            then
                af="NA"
            else
                cmd="echo \"$al_count / ($gt_count * 2 ) \" | bc -l"
                aff=` eval "$cmd" `
                af=` printf "%6.4f" "$aff" `
            fi
            cmd="echo \"$gt_count / ($col_count ) \" | bc -l"
            gf=` eval "$cmd" `
            cmd="echo \"$al_count / ($col_count *2 ) \" | bc -l"
            pf=` eval "$cmd" `
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%6.4f\t%s\t%6.4f\n" "$rec_out" "$wt_count" "$het_count" "$hom_count" "$oth_count" "$na_count" "$gt_count" "$gf" "$af" "$pf"
        done
    done
}

col_prefix=$( echo $dataset_name | awk '{print toupper($0)}' )
# create header
header="#Chr"
header+="\tStart"
header+="\tEnd"
header+="\tRef"
header+="\tAlt"
header+="\t$col_prefix"_WT
header+="\t$col_prefix"_HET
header+="\t$col_prefix"_HOM
header+="\t$col_prefix"_OTH
header+="\t$col_prefix"_NA
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

#---------- idx stat file --------------
idx_stat_file="$out_file.idx"
idx_cmd="$PYCMM/bash/compileAnnnovarIndex.pl"
idx_cmd+=" $out_file"
idx_cmd+=" 1000"
idx_cmd+=" > $idx_stat_file"
eval_cmd "$idx_cmd"
#---------- idx stat file --------------

new_section_txt "F I N I S H <$script_name>"
