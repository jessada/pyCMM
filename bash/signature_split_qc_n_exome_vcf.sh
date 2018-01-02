#!/bin/bash 

#set -e
#set -u
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
-r {file}         list of chromosome regions (required)
-o {dir}          project output directory (required)
-h                this help
EOF
)

# parse option
while getopts ":i:r:o:h" OPTION; do
  case "$OPTION" in
    i)
      input_vcf="$OPTARG"
      ;;
    r)
      regions_list="$OPTARG"
      ;;
    o)
      project_out_dir="$OPTARG"
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
[ ! -z $regions_list ] || die "chromosome regions list is required (-r)"
[ ! -z $project_out_dir ] || die "please indicate project output file directory (-o)"
[ -f "$input_vcf" ] || die "$input_vcf is not a valid file name."

time_stamp=$( date )

cd $PYCMM
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

slurm_log_dir="$project_out_dir/slurm_log"
if [ ! -d "$slurm_log_dir" ]
then
    mkdir -p "$slurm_log_dir"
fi
qced_vcf_out_dir="$project_out_dir/qced_vcf"
if [ ! -d "$qced_vcf_out_dir" ]
then
    mkdir -p "$qced_vcf_out_dir"
fi
chunk_vcf_out_dir="$project_out_dir/chunk_vcf"
if [ ! -d "$chunk_vcf_out_dir" ]
then
    mkdir -p "$chunk_vcf_out_dir"
fi
list_file_out_dir="$project_out_dir/list_file"
if [ ! -d "$list_file_out_dir" ]
then
    mkdir -p "$list_file_out_dir"
fi

file_name=$(basename "$input_vcf")
file_prefix="${file_name%%.*}"

qc_n_exome_script="$PYCMM/bash/signature_qc_n_exome_vcf.sh"
concat_vcfs_script="$PYCMM/bash/concat_vcfs.sh"

out_vcf_gz="$qced_vcf_out_dir/$file_prefix.qc_n_exome.vcf.gz"

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
display_param "chromosome regions list (-r)" "$regions_list"
display_param "project output directory (-o)" "$project_out_dir"
display_param "  - slurm log directory" "$slurm_log_dir"
display_param "  - QCed vcf output directory" "$qced_vcf_out_dir"
display_param "  - chunk vcf output directory" "$chunk_vcf_out_dir"
display_param "  - list file output directory" "$list_file_out_dir"
display_param "qc_n_exome output file" "$out_vcf_gz"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************
 
new_section_txt "Split vcf file into small small chunks"

qced_vcf_list="$list_file_out_dir/qced_vcf.list"
> $qced_vcf_list
job_id_list="$list_file_out_dir/job_id.list"
> $job_id_list

while read region
do
    region_txt=`echo $region | sed 's/:/_/' | sed 's/-/_/'`
    chunk_vcf_file="$chunk_vcf_out_dir/$file_prefix.$region_txt.vcf"
    qced_vcf_file="$qced_vcf_out_dir/$file_prefix.$region_txt.qced.vcf"

    cmd="tabix -h"
    cmd+=" $input_vcf"
    cmd+=" $region"
    cmd+=" > $chunk_vcf_file"
    eval_cmd "$cmd"
    
    cmd="bgzip -f"
    cmd+=" $chunk_vcf_file"
    eval_cmd "$cmd"
    
    cmd="tabix"
    cmd+=" -p vcf"
    cmd+=" $chunk_vcf_file.gz"
    eval_cmd "$cmd"

    job_name="qc_exome.$file_prefix.$region_txt"
    std_out_slurm="$slurm_log_dir/$job_name.%j.out.log"
    std_err_slurm="$slurm_log_dir/$job_name.%j.err.log"

    submit_job_cmd="sbatch"
    submit_job_cmd+=" -A $SLURM_JOB_ACCOUNT"
    submit_job_cmd+=" -p core"
    submit_job_cmd+=" -n 8"
    submit_job_cmd+=" -t 4-00:00:00"
    submit_job_cmd+=" -J $job_name"
    submit_job_cmd+=" -o $std_out_slurm"
    submit_job_cmd+=" -e $std_err_slurm"
    submit_job_cmd+=" $qc_n_exome_script"
    submit_job_cmd+=" -i $chunk_vcf_file.gz"
    submit_job_cmd+=" -o $qced_vcf_file"
    job_info=$( eval_cmd "$submit_job_cmd" )
    job_id=${job_info##* }

    echo "$qced_vcf_file.gz" >> $qced_vcf_list
    echo "$job_id" >> $job_id_list
done < $regions_list


new_section_txt "Submit job to concat all the vcf files back"

job_name="concat_qced_vcf_$file_prefix"
std_out_slurm="$slurm_log_dir/$job_name.%j.out.log"
std_err_slurm="$slurm_log_dir/$job_name.%j.err.log"

submit_job_cmd="sbatch"
submit_job_cmd+=" -A $SLURM_JOB_ACCOUNT"
submit_job_cmd+=" -p core"
submit_job_cmd+=" -n 8"
submit_job_cmd+=" -t 4-00:00:00"
submit_job_cmd+=" -J $job_name"
submit_job_cmd+=" -o $std_out_slurm"
submit_job_cmd+=" -e $std_err_slurm"
submit_job_cmd+=" --dependency=afterok"
while read job_id
do
    submit_job_cmd+=":$job_id"
done < $job_id_list
submit_job_cmd+=" $concat_vcfs_script"
submit_job_cmd+=" -i $qced_vcf_list"
submit_job_cmd+=" -o $out_vcf_gz"
eval_cmd "$submit_job_cmd"

new_section_txt "F I N I S H <$script_name>"
