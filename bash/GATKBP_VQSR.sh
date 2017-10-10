#!/bin/bash 
set -e
set -o pipefail

source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}         input VCF file (required)
-o {file}         output VCF file with VQSR (required)
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
[ ! -z $out_file ] || die "please indicate output VQSR VCF file name (-o)"
[ -f "$input_vcf" ] || die "$input_vcf is not a valid file name."

working_dir=`mktemp -d`
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
info_msg "  Perform Variant Quality Score Recalibration of VCF file"
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
display_param "output VQSR VCF file (-o)" "$out_file"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************
 
snp_recal_file="$working_dir/sample.snp.recal"
snp_tranches_file="$working_dir/sample.snp.tranches"
indel_recal_file="$working_dir/sample.indel.recal"
indel_tranches_file="$working_dir/sample.indel.tranches"
out_filtered_snp_file="$working_dir/sample.recal.snp.filtered.vcf"
out_filtered_indel_file="$out_file"

resource_dir="$REFERENCE_DATA_DIR"
ref="$GRCH37_REF"
hapmap_vcf="$resource_dir/hapmap_3.3.b37.vcf"
KG_omni="$resource_dir/1000G_omni2.5.b37.vcf"
KG_phase1="$resource_dir/1000G_phase1.snps.high_confidence.b37.vcf"
db_snp="$resource_dir/dbsnp_138.b37.vcf"
Mills="$resource_dir/Mills_and_1000G_gold_standard.indels.b37.vcf"

cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T VariantRecalibrator"
cmd+=" -R $ref"
cmd+=" -input $input_vcf"
cmd+=" -recalFile $snp_recal_file"
cmd+=" -tranchesFile $snp_tranches_file"
cmd+=" -nt 4"
cmd+=" -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf"
cmd+=" -resource:omni,known=false,training=true,truth=true,prior=12.0 $KG_omni"
cmd+=" -resource:1000G,known=false,training=true,truth=false,prior=10.0 $KG_phase1"
cmd+=" -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $db_snp"
cmd+=" -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff"
#cmd+=" -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP"
cmd+=" -mode SNP"
eval_cmd "$cmd"

cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T VariantRecalibrator"
cmd+=" -R $ref"
cmd+=" -input $input_vcf"
cmd+=" -recalFile $indel_recal_file"
cmd+=" -tranchesFile $indel_tranches_file"
cmd+=" -nt 4"
cmd+=" --maxGaussians 4"
cmd+=" -resource:mills,known=false,training=true,truth=true,prior=12.0 $Mills"
cmd+=" -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $db_snp"
cmd+=" -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff "
#cmd+=" -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum"
cmd+=" -mode INDEL"
eval_cmd "$cmd"

cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T ApplyRecalibration"
cmd+=" -R $ref"
cmd+=" -input $input_vcf"
cmd+=" -recalFile $snp_recal_file"
cmd+=" -tranchesFile $snp_tranches_file"
cmd+=" -o $out_filtered_snp_file"
cmd+=" --ts_filter_level 99.5"
cmd+=" -mode SNP"
eval_cmd "$cmd"

cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T ApplyRecalibration"
cmd+=" -R $ref"
cmd+=" -input $out_filtered_snp_file"
cmd+=" -recalFile $indel_recal_file"
cmd+=" -tranchesFile $indel_tranches_file"
cmd+=" -o $out_filtered_indel_file"
cmd+=" --ts_filter_level 99.0"
cmd+=" -mode INDEL"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
