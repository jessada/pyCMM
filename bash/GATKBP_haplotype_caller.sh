#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-I {BAM file}       BAM file (required)
-R {file}           reference file (required)
-d {file}           dbSNP file
-L {file}           target interval list
-o {file}           output GVCF file (required)
-h                  this help
EOF
)

# parse option
while getopts ":I:R:d:L:o:h" OPTION; do
  case "$OPTION" in
    I)
      bam_file="$OPTARG"
      ;;
    R)
      ref="$OPTARG"
      ;;
    d)
      dbsnp="$OPTARG"
      ;;
    L)
      targets_interval_list="$OPTARG"
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

[ ! -z $bam_file ] || die "a BAM input file is required (-I)"
[ ! -z $ref ] || die "reference file is required (-R)"
[ ! -z $out_file ] || die "please indicate output file name (-o)"
[ -f "$bam_file" ] || die "$bam_file is not found"
[ -f "$ref" ] || die "$ref is not found"

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
info_msg "  Call SNPs and indels simultaneously via local re-assembly of haplotypes in an active region"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "BAM input file (-I)" "$bam_file"
display_param "reference file (-R)" "$ref"
if [ ! -z "$targets_interval_list" ]
then
    display_param "target interval list (-L)" "$targets_interval_list"
fi
if [ ! -z "$dbsnp" ]
then
    display_param "dbSNP (-d)" "$dbsnp"
fi
display_param "output GVCF file (-o)" "$out_file"

# ****************************************  executing  ****************************************
 
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T HaplotypeCaller"
cmd+=" -R $ref"
cmd+=" -I $bam_file"
cmd+=" --emitRefConfidence GVCF"
cmd+=" --variant_index_type LINEAR"
cmd+=" --variant_index_parameter 128000"
#cmd+=" -allowPotentiallyMisencodedQuals"
if [ ! -z "$dbsnp" ]
then
    cmd+=" --dbsnp $dbsnp"
fi
if [ ! -z "$targets_interval_list" ]
then
    cmd+=" -L $targets_interval_list"
fi
cmd+=" -o $out_file"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
