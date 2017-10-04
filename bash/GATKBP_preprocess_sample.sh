#!/bin/bash 
set -e
set -o pipefail

module load bwa
module load picard
module load GATK

source $PYCMM/bash/cmm_functions.sh

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-I {fastq list}     list of fastq files (comma separated) (required)
-R {file}           reference file (required)
-k {file}           known indels (required)
-d {file}           dbSNP file (required)
-L {file}           target interval list
-g {sample group}   sample group (required) 
-n {sample id}      sample id (required)
-w {dir}            working directory (default: $TMPDIR)
-b {file}           output bam file (required)
-o {file}           output gvcf file (required)
-h                  this help
EOF
)

# parse option
while getopts ":I:R:k:d:L:g:n:w:b:o:h" OPTION; do
  case "$OPTION" in
    I)
      sample_list="$OPTARG"
      ;;
    R)
      ref="$OPTARG"
      ;;
    k)
      known_indels="$OPTARG"
      ;;
    d)
      dbsnp="$OPTARG"
      ;;
    L)
      targets_interval_list="$OPTARG"
      ;;
    g)
      sample_group="$OPTARG"
      ;;
    n)
      sample_id="$OPTARG"
      ;;
    w)
      working_dir="$OPTARG"
      ;;
    b)
      out_recal_bam="$OPTARG"
      ;;
    o)
      out_gvcf="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $sample_list ] || die "list of fastq files is required (-I)"
[ ! -z $sample_group ] || die "sample group is required (-g)"
[ ! -z $sample_id ] || die "sample name is required (-n)"
[ ! -z $out_recal_bam ] || die "output bam file is required (-b)"
[ ! -z $out_gvcf ] || die "output gvcf file is required (-o)"
[ ! -z $ref ] || die "reference file is required (-R)"
[ ! -z $known_indels ] || die "known indel file(s) are required (-k)"
[ ! -z $dbsnp ] || die "dbSNP file is required (-d)"
IFS=',' read -ra fastq_files <<< "$sample_list"
sample1="${fastq_files[0]}"
sample2="${fastq_files[1]}"
IFS=',' read -ra known_indels_files <<< "$known_indels"
for (( indel_idx=0; indel_idx<$((${#known_indels_files[@]})); indel_idx++ ))
do
    [ -f "${known_indels_files[$indel_idx]}" ] || die "${known_indels_files[$indel_idx]} is not found"
done
[ -f "$sample1" ] || die "$sample1 is not found"
[ -f "$ref" ] || die "$ref is not found"
[ -f "$dbsnp" ] || die "$dbsnp is not found"

time_stamp=$( date )
if [ -z $working_dir ]
then
    working_dir=`mktemp -d`
fi

tmp_bwa_mem_out="$working_dir/sample.sam"
tmp_sorted_bam="$working_dir/sample.bam"
tmp_dedup_bam="$working_dir/sample.dedup.bam"
tmp_metrics_file="$working_dir/sample.metics.txt"
tmp_realign_target_intevals="$working_dir/sample.target.intervals"
tmp_realigned_bam="$working_dir/sample.dedup.realigned.bam"
tmp_recal_table="$working_dir/sample.recal.table"

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  Preprocess sample based on GATK Best Practice workflow"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "script path" "$PYCMM"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "sample prefix (-I)" "$sample_prefix"
display_param "  sample 1" "$sample1"
display_param "  sample 2" "$sample2"
display_param "sample group (-g)" "$sample_group"
display_param "sample name (-n)" "$sample_id"
display_param "reference file (-R)" "$ref"
display_param "known indels (-k)"
known_indels_param=""
known_sites=""
for (( indel_idx=0; indel_idx<$((${#known_indels_files[@]})); indel_idx++ ))
do
    display_param "  known indels file #$(( indel_idx+1 )) " "${known_indels_files[$indel_idx]}"
    known_indels_param+=" -known ${known_indels_files[$indel_idx]}"
    known_sites+=" -knownSites ${known_indels_files[$indel_idx]}"
done
display_param "dbSNP file (-d)" "$dbsnp"
if [ ! -z "$targets_interval_list" ]
then
    display_param "target interval list (-L)" "$targets_interval_list"
fi
display_param "working directory (-w)" "$working_dir"
display_param "  bwa-mem output" "$tmp_bwa_mem_out"
display_param "  sorted bam file" "$tmp_sorted_bam"
display_param "  dedup bam file" "$tmp_dedup_bam"
display_param "  realign target intervals" "$tmp_realign_target_intevals"
display_param "  realigned bam file" "$tmp_realigned_bam"
display_param "recal bam file" "$out_recal_bam"
display_param "output gvcf file" "$out_gvcf"

# ****************************************  executing  ****************************************

new_section_txt "executing bwa-mem"
read_group="$sample_group.$sample_id"

cmd=" bwa mem -M -R '@RG\tID:$read_group\tSM:$sample_id\tPL:illumina\tLB:lib1\tPU:unit1'"
cmd+=" $ref"
cmd+=" $sample1"
if [ -f $sample2 ]
then
    cmd+=" $sample2"
fi
cmd+=" > $tmp_bwa_mem_out" 
eval_cmd "$cmd"

new_section_txt "executing SortSam"
cmd=" java -jar $CLASSPATH/picard.jar SortSam"
cmd+=" INPUT=$tmp_bwa_mem_out"
cmd+=" OUTPUT=$tmp_sorted_bam"
cmd+=" SORT_ORDER=coordinate"
eval_cmd "$cmd"

new_section_txt "executing MarkDuplicates"
cmd=" java -jar $CLASSPATH/picard.jar MarkDuplicates"
cmd+=" INPUT=$tmp_sorted_bam"
cmd+=" OUTPUT=$tmp_dedup_bam"
cmd+=" METRICS_FILE=$tmp_metrics_file"
eval_cmd "$cmd"

cmd=" java -jar $CLASSPATH/picard.jar BuildBamIndex"
cmd+=" INPUT=$tmp_dedup_bam"
eval_cmd "$cmd"

new_section_txt "Create a target list of intervals to be realigned"
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T RealignerTargetCreator"
cmd+=" -I $tmp_dedup_bam"
cmd+=" -R $ref"
cmd+=" $known_indels_param"
cmd+=" -o $tmp_realign_target_intevals"
eval_cmd "$cmd"

new_section_txt "Perform realignment of the target intervals"
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T IndelRealigner"
cmd+=" -R $ref"
cmd+=" -I $tmp_dedup_bam"
cmd+=" -targetIntervals $tmp_realign_target_intevals"
cmd+=" $known_indels_param"
cmd+=" -o $tmp_realigned_bam"
eval_cmd "$cmd"

new_section_txt "Analyze patterns of covariation in the sequence dataset"
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T BaseRecalibrator"
cmd+=" -R $ref"
cmd+=" -I $tmp_realigned_bam"
cmd+=" $known_sites"
cmd+=" -o $tmp_recal_table"
eval_cmd "$cmd"

new_section_txt "Apply the recalibration to your sequence data"
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T PrintReads"
cmd+=" -R $ref"
cmd+=" -I $tmp_realigned_bam"
cmd+=" -BQSR $tmp_recal_table"
cmd+=" -o $out_recal_bam"
eval_cmd "$cmd"

new_section_txt " executing HaplotypeCaller"
cmd=" java -jar $GATK_dir/GenomeAnalysisTK.jar"
cmd+=" -T HaplotypeCaller"
cmd+=" -R $ref"
cmd+=" -I $out_recal_bam"
cmd+=" --emitRefConfidence GVCF"
cmd+=" --variant_index_type LINEAR"
cmd+=" --variant_index_parameter 128000"
cmd+=" --disable_auto_index_creation_and_locking_when_reading_rods"
if [ ! -z "$dbsnp" ]
then
    cmd+=" --dbsnp $dbsnp"
fi
if [ ! -z "$targets_interval_list" ]
then
    cmd+=" -L $targets_interval_list"
fi
cmd+=" -o $out_gvcf"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"

