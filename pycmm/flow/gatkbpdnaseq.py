import sys
from os.path import join as join_path
from pycmm.utils import mylogger
from pycmm.utils.jobman import JobManager
from pycmm.settings import GATK_ALLOC_TIME

jobman = JobManager()

#def cancel_job(job_name):
#    job_status = jobman.get_job_status(job_name)
#    jobman.cancel_job(job_name)
#
def bwa_mem(job_name,
            project_code,
            slurm_log_file,
            sample_prefix,
            sample_group,
            ref,
            out_raw_aligned_reads,
            email=False,
            ):
    """ To generate a SAM file containing aligned reads """

    job_script = "$PYCMM/bash/GATKBP_bwa_mem.sh"
    job_params = " -S " + sample_prefix 
    job_params += " -g " + sample_group
    job_params += " -r " + ref
    job_params += " -o " + out_raw_aligned_reads
    jobman.submit_job(job_name,
                      project_code,
                      "core",
                      "1",
                      GATK_ALLOC_TIME,
                      slurm_log_file,
                      job_script,
                      job_params,
                      email=email,
                      )

def sort_sam(job_name,
               project_code,
               slurm_log_file,
               sam_file,
               out_sorted_reads,
               email=False,
               prerequisite=None,
               ):
    """ Use picard command to sort the SAM file and convert it to BAM """

    job_script = "$PYCMM/bash/GATKBP_sort_sam.sh"
    job_params = " -S " + sam_file 
    job_params += " -o " + out_sorted_reads
    jobman.submit_job(job_name,
                      project_code,
                      "core",
                      "2",
                      GATK_ALLOC_TIME,
                      slurm_log_file,
                      job_script,
                      job_params,
                      email=email,
                      prerequisite=prerequisite,
                      )

def mark_dup(job_name,
             project_code,
             slurm_log_file,
             bam_file,
             out_dedup_reads,
             out_metrics_file,
             email=False,
             prerequisite=None,
             ):
    """ Use picard command to mark duplicates """

    job_script = "$PYCMM/bash/GATKBP_mark_dup.sh"
    job_params = " -B " + bam_file 
    job_params += " -o " + out_dedup_reads
    job_params += " -m " + out_metrics_file
    jobman.submit_job(job_name,
                      project_code,
                      "core",
                      "4",
                      GATK_ALLOC_TIME,
                      slurm_log_file,
                      job_script,
                      job_params,
                      email=email,
                      prerequisite=prerequisite,
                      )

def haplotype_caller(job_name,
                     project_code,
                     slurm_log_file,
                     bam_file,
                     out_gvcf,
                     ref,
                     targets_interval_list=None,
                     email=False,
                     prerequisite=None,
                     ):
    """
    Call SNPs and indels simultaneously via local re-assembly of haplotypes
    in an active region
    """

    job_script = "$PYCMM/bash/GATKBP_haplotype_caller.sh"
    job_params = " -B " + bam_file 
    job_params += " -o " + out_gvcf
    job_params += " -r " + ref
    if targets_interval_list is not None:
        job_params += " -L " + targets_interval_list
    jobman.submit_job(job_name,
                      project_code,
                      "core",
                      "2",
                      GATK_ALLOC_TIME,
                      slurm_log_file,
                      job_script,
                      job_params,
                      email=email,
                      prerequisite=prerequisite,
                      )

def combine_gvcfs(job_name,
                  project_code,
                  slurm_log_file,
                  gvcf_list,
                  out_merged_gvcf,
                  ref,
                  working_dir,
                  email=False,
                  prerequisite=None,
                  ):
    """
    Combines any number of gVCF files that were produced by
    the Haplotype Caller into a single joint gVCF file
    """

    job_script = "$PYCMM/bash/GATKBP_combine_gvcfs.sh"
    gvcf_list_file = join_path(working_dir,
                               job_name+"_gvcf.list")
    f_gvcf = open(gvcf_list_file, "w")
    map(lambda x: f_gvcf.write(x+"\n"), gvcf_list)
    job_params = " -G " + gvcf_list_file 
    job_params += " -o " + out_merged_gvcf
    job_params += " -r " + ref
    jobman.submit_job(job_name,
                      project_code,
                      "core",
                      "4",
                      GATK_ALLOC_TIME,
                      slurm_log_file,
                      job_script,
                      job_params,
                      email=email,
                      prerequisite=prerequisite,
                      )

def preprocess_sample(project_code,
                      sample_name,
                      sample_prefix,
                      sample_group,
                      ref,
                      out_dir,
                      working_dir,
                      slurm_log_dir,
                      time_stamp,
                      out_merged_gvcf,
                      email=False,
                      ):
    """
    An aggregate flow to pre-process a DNASeq sample. The flow consists of
      - bwa-mem alignment
      - reads sorting
      - marking duplicates
      - haplotype caller
    The flow output a gvcf file for one sample. The return value is the last
    job name of the process.
    """

    #define shared variable
    log_file_suffix = "_"
    if type(time_stamp) is str:
        log_file_suffix += time_stamp
    else:
        log_file_suffix += time_stamp.strftime("%Y%m%d%H%M%S")
    log_file_suffix += ".log"
    raw_aligned_reads = join_path(working_dir,
                                  sample_name+"_raw_aligned.sam")
    sorted_reads = join_path(working_dir,
                             sample_name+"_sorted_reads.bam")
    dedup_reads = join_path(out_dir,
                            sample_name+"_dedup_reads.bam")
    metrics_file = join_path(working_dir,
                             sample_name+"_metrics.txt")

    # running bwa-mem alignment
    job_name_bwa_mem = sample_name + "_bwa_mem"
    slurm_log_file = join_path(slurm_log_dir,
                               job_name_bwa_mem+log_file_suffix)
    bwa_mem(job_name_bwa_mem,
            project_code,
            slurm_log_file,
            sample_prefix,
            sample_group,
            ref,
            raw_aligned_reads,
            email=email,
            )

    # running reads sorting
    job_name_sort_sam = sample_name + "_sort_sam"
    slurm_log_file = join_path(slurm_log_dir,
                               job_name_sort_sam+log_file_suffix)
    sort_sam(job_name_sort_sam,
             project_code,
             slurm_log_file,
             raw_aligned_reads,
             sorted_reads,
             email=email,
             prerequisite=[job_name_bwa_mem],
             )

    # running marking duplicates
    job_name_mark_dup = sample_name + "_mark_dup"
    slurm_log_file = join_path(slurm_log_dir,
                               job_name_mark_dup+log_file_suffix)
    mark_dup(job_name_mark_dup,
             project_code,
             slurm_log_file,
             sorted_reads,
             dedup_reads,
             metrics_file,
             email=email,
             prerequisite=[job_name_sort_sam],
             )

    # running haplotype caller
    job_name_hap_call = sample_name + "_hap_call"
    slurm_log_file = join_path(slurm_log_dir,
                               job_name_hap_call+log_file_suffix)
    haplotype_caller(job_name_hap_call,
                     project_code,
                     slurm_log_file,
                     dedup_reads,
                     out_merged_gvcf,
                     ref,
                     email=email,
                     prerequisite=[job_name_mark_dup],
                     )
    return job_name_mark_dup

def preprocess_dataset(project_code,
                       sample_list,
                       ref,
                       out_dir,
                       working_dir,
                       slurm_log_dir,
                       time_stamp,
                       out_gvcf,
                       email=False,
                       ):
    """
    An aggregate flow to pre-process a DNASeq sample. The flow consists of
      - bwa-mem alignment
      - reads sorting
      - marking duplicates
      - haplotype caller
    The flow output a gvcf file for one sample. The return value is the last
    job name of the process.
    """

    return None
