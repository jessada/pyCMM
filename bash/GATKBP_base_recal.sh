#!/bin/bash 
source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-I {file}     input bam file to be realigned(required)
-R {file}     reference file (required)
-k {file}     known indels (required)
-t {file}     output recalibration table (required)
-a {file}     output post-recalibration table (required)
-p {file}     pdf output recalibration plot (required)
-o {file}     output recalibrated bam file (required)
-h            this help
EOF
)

# parse option
while getopts ":I:R:k:t:a:p:o:h" OPTION; do
  case "$OPTION" in
    I)
      bam_file="$OPTARG"
      ;;
    R)
      ref="$OPTARG"
      ;;
    k)
      known_indels="$OPTARG"
      ;;
    t)
      recal_table="$OPTARG"
      ;;
    a)
      post_recal_table="$OPTARG"
      ;;
    p)
      pdf_plots="$OPTARG"
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

[ ! -z $bam_file ] || die "a bam file to be re-aligned is required (-I)"
[ ! -z $ref ] || die "reference file is required (-R)"
[ ! -z $known_indels ] || die "list of known indels in vcf format is required (-k)"
[ ! -z $recal_table ] || die "please specify the output file name of recalibration table (-t)"
[ ! -z $post_recal_table ] || die "please specify the output file name of post-recalibration table (-a)"
[ ! -z $pdf_plots ] || die "output file name of recalibration plots is required (-p)"
[ ! -z $out_file ] || die "bam output file name is required (-o)"
[ -f "$bam_file" ] || die "$bam_file is not found"
[ -f "$ref" ] || die "$ref is not found"

IFS=',' read -ra known_indels_files <<< "$known_indels"
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
info_msg "  Base Quality Score Recalibration"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "bam file (-I)" "$bam_file"
display_param "reference file (-R)" "$ref"
display_param "known indels (-k)"
known_indels_param=""
for (( indel_idx=0; indel_idx<$((${#known_indels_files[@]})); indel_idx++ ))
do
    display_param "  known indels file #$(( indel_idx+1 )) " "${known_indels_files[$indel_idx]}"
    known_indels_param+=" -knownSites ${known_indels_files[$indel_idx]}"
done
display_param "output recalibration table (-t)" "$recal_table"
display_param "output post-recalibration table (-a)" "$post_recal_table"
display_param "pdf recalibration plots (-p)" "$pdf_plots"
display_param "bam output file (-o)" "$out_file"
# ****************************************  executing  ****************************************
# 1. Analyze patterns of covariation in the sequence dataset
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T BaseRecalibrator"
cmd+=" -R $ref"
cmd+=" -I $bam_file"
cmd+=" $known_indels_param"
cmd+=" -o $recal_table"
eval_cmd "$cmd"

## 2. Do a second pass to analyze covariation remaining after recalibration
#cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
#cmd+=" -T BaseRecalibrator"
#cmd+=" -R $ref"
#cmd+=" -I $bam_file"
#cmd+=" $known_indels_param"
#cmd+=" -BQSR $recal_table"
#cmd+=" -o $post_recal_table"
#eval_cmd "$cmd"
#
## 3. Generate before/after plots
#cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
#cmd+=" -T AnalyzeCovariates"
#cmd+=" -R $ref"
#cmd+=" -before $recal_table"
#cmd+=" -after $post_recal_table"
#cmd+=" -plots $pdf_plots"
#eval_cmd "$cmd"

# 4. Apply the recalibration to your sequence data
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T PrintReads"
cmd+=" -R $ref"
cmd+=" -I $bam_file"
cmd+=" -BQSR $recal_table"
cmd+=" -o $out_file"
eval_cmd "$cmd"
new_section_txt "F I N I S H <$script_name>"

