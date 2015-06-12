import sys
import re
import datetime
import pyaml
import yaml
from os import listdir
from os.path import join as join_path
from os.path import isdir
from os.path import isfile
from pycmm.template import pyCMMBase
from pycmm.utils import mylogger
from pycmm.utils.jobman import JobManager
from pycmm.utils.jobman import JOB_STATUS_COMPLETED
from pycmm.utils import exec_sh
from pycmm.settings import GATK_ALLOC_TIME

METADATA_PATTERN = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')

BWA_MEM_SCRIPT = "$PYCMM/bash/GATKBP_bwa_mem.sh"
CONCAT_FILES_SCRIPT = "$PYCMM/bash/concat_files.sh"
SORT_SAM_SCRIPT = "$PYCMM/bash/GATKBP_sort_sam.sh"
MARK_DUP_SCRIPT = "$PYCMM/bash/GATKBP_mark_dup.sh"
INDELS_TARGET_INTERVALS_SCRIPT = "$PYCMM/bash/GATKBP_create_intervals.sh"
INDELS_REALIGN_SCRIPT = "$PYCMM/bash/GATKBP_indels_realign.sh"
BASE_RECAL_SCRIPT = "$PYCMM/bash/GATKBP_base_recal.sh"
HAPLOTYPE_CALLER_SCRIPT = "$PYCMM/bash/GATKBP_haplotype_caller.sh"
COMBINE_GVCFS_SCRIPT = "$PYCMM/bash/GATKBP_combine_gvcfs.sh"
GENOTYPE_GVCFS_SCRIPT = "$PYCMM/bash/GATKBP_genotype_gvcfs.sh"

# jobs metadata section
JOBS_SETUP_DATASET_NAME_KEY = "DATASET_NAME"
JOBS_SETUP_PROJECT_CODE_KEY = "PROJECT_CODE"
JOBS_SETUP_KNOWN_INDELS_KEY = "KNOWN_INDELS"
JOBS_SETUP_DBSNP_KEY = "DBSNP"
JOBS_SETUP_REFFERENCE_KEY = "REFERENCE"
JOBS_SETUP_OUTPUT_DIR_KEY = "OUTPUT_DIR"
JOBS_SETUP_VARIANTS_CALLING_KEY = "VARIANTS_CALLING"
JOBS_SETUP_JOBS_REPORT_FILE_KEY = "JOBS_REPORT_FILE"
JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY = "TARGETS_INTERVAL_LIST"
JOBS_SETUP_DATASET_USAGE_MAIL_KEY = "USAGE_MAIL"

# jobs sample info section
JOBS_SETUP_SAMPLES_SECTION = "SAMPLES"
JOBS_SETUP_SAMPLE_NAME_KEY = "sample_name"
JOBS_SETUP_FASTQ_PAIRS_KEY = "fastq_files"
JOBS_SETUP_R1_KEY = "R1"
JOBS_SETUP_R2_KEY = "R2"
JOBS_SETUP_SAMPLE_GROUP_KEY = "sample_group"
JOBS_SETUP_PLATFORM_KEY = "platform"
JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY = "usage_mail"
JOBS_SETUP_PREPROCESS_SAMPLE_KEY = "preprocess_sample"

class SampleRecord(pyCMMBase):
    """ A structure to parse and keep sample information """

    def __init__(self,
                 sample_info,
                 output_dir,
                 ):
        self.__parse_sample_info(sample_info)
        self.__output_dir = output_dir

    def get_raw_repr(self):
        return {"sample name": self.sample_name,
                "fastq pairs": self.fastq_pairs,
                "sample group": self.sample_group,
                "platform": self.platform,
                "library": self.library,
                "unit": self.unit,
                "usage limit": self.usage_limit,
                "pre-process sample": self.preprocess_sample,
                }

    def __parse_sample_info(self, sample_info):
        self.__sample_name = sample_info[JOBS_SETUP_SAMPLE_NAME_KEY]
        self.__platform = sample_info[JOBS_SETUP_PLATFORM_KEY]
        self.__sample_group = sample_info[JOBS_SETUP_SAMPLE_GROUP_KEY]
        self.__fastq_pairs = sample_info[JOBS_SETUP_FASTQ_PAIRS_KEY]
        self.__usage_mail = sample_info[JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY]
        self.__preprocess_sample = sample_info[JOBS_SETUP_PREPROCESS_SAMPLE_KEY]

    @property
    def sample_name(self):
        return self.__sample_name

    @property
    def sample_group(self):
        return self.__sample_group

    @property
    def platform(self):
        return self.__platform

    @property
    def library(self):
        return 'lib1'

    @property
    def unit(self):
        return 'unit1'

    @property
    def usage_mail(self):
        return self.__usage_mail

    @property
    def preprocess_sample(self):
        return self.__preprocess_sample

    @property
    def fastq_pairs(self):
        return self.__fastq_pairs

    @property
    def raw_aligned_reads_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_raw_aligned_reads.sam")

    @property
    def sorted_reads_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_sorted_reads.bam")

    @property
    def dedup_reads_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_dedup_reads.bam")

    @property
    def metrics_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_metrics.txt")

    @property
    def indels_target_intervals_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_indels_target_intervals.list")

    @property
    def realigned_reads_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_realigned_reads.bam")

    @property
    def recal_table_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_recal_data.table")

    @property
    def post_recal_table_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_post_recal_data.table")

    @property
    def recalibration_plots(self):
        return join_path(join_path(self.__output_dir,
                                   'pdf'),
                         self.sample_name+"_recalibration_plots.pdf")

    @property
    def recal_reads_file(self):
        return join_path(join_path(self.__output_dir,
                                   'bam'),
                         self.sample_name+"_recal_reads.bam")

    @property
    def raw_gvcf_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_raw.g.vcf")


class GATKBPPipeline(JobManager):
    """ A class to control GATK best practice pipeline """

    def __init__(self,
                 jobs_setup_file,
                 ):
        mylogger.getLogger(__name__)
        self.__load_jobs_info(jobs_setup_file)
        JobManager.__init__(self,
                            jobs_report_file=self.jobs_report_file)
        self.__vcf_out_dir = join_path(self.output_dir,
                                       "vcf")
        self.__bam_out_dir = join_path(self.output_dir,
                                       "bam")
        self.__pdf_out_dir = join_path(self.output_dir,
                                       "pdf")
        self.__samples_working_dir = join_path(self.output_dir,
                                               "tmp")
        self.__slurm_log_dir = join_path(self.output_dir,
                                         "slurm_log")
        self.__time_stamp = datetime.datetime.now()
        self.__create_directories()

    def get_raw_repr(self):
        return {"dataset name": self.dataset_name,
                "project code": self.project_code,
                "known indels": self.known_indels,
                "reference": self.reference,
                "output directory": self.output_dir,
                "varaints calling": self.variants_calling,
                "jobs report file": self.jobs_report_file,
                "targets interval list": self.targets_interval_list,
                "usage_mail": self.usage_mail,
                }

    @property
    def dataset_name(self):
        return self.__meta_data[JOBS_SETUP_DATASET_NAME_KEY]

    @property
    def project_code(self):
        return self.__meta_data[JOBS_SETUP_PROJECT_CODE_KEY]

    @property
    def known_indels(self):
        if JOBS_SETUP_KNOWN_INDELS_KEY in self.__meta_data:
            return self.__meta_data[JOBS_SETUP_KNOWN_INDELS_KEY]
        else:
            return None

    @property
    def dbsnp(self):
        if JOBS_SETUP_DBSNP_KEY in self.__meta_data:
            return self.__meta_data[JOBS_SETUP_DBSNP_KEY]
        else:
            return None

    @property
    def reference(self):
        return self.__meta_data[JOBS_SETUP_REFFERENCE_KEY]

    @property
    def output_dir(self):
        return self.__meta_data[JOBS_SETUP_OUTPUT_DIR_KEY]

    @property
    def variants_calling(self):
        return self.__meta_data[JOBS_SETUP_VARIANTS_CALLING_KEY]

    @property
    def jobs_report_file(self):
        return self.__meta_data[JOBS_SETUP_JOBS_REPORT_FILE_KEY]

    @property
    def targets_interval_list(self):
        if JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY in self.__meta_data:
            return self.__meta_data[JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY]
        else:
            return None

    @property
    def usage_mail(self):
        return self.__meta_data[JOBS_SETUP_DATASET_USAGE_MAIL_KEY]

    @property
    def samples(self):
        return self.__samples

    @property
    def vcf_out_dir(self):
        return self.__vcf_out_dir

    @property
    def bam_out_dir(self):
        return self.__bam_out_dir

    @property
    def pdf_out_dir(self):
        return self.__pdf_out_dir

    @property
    def slurm_log_dir(self):
        return self.__slurm_log_dir

    @property
    def samples_working_dir(self):
        return self.__samples_working_dir

    @property
    def time_stamp(self):
        return self.__time_stamp

    @property
    def combined_gvcfs_file(self):
        return join_path(self.samples_working_dir,
                         self.dataset_name+"_cohort.g.vcf")

    @property
    def genotyped_gvcfs_file(self):
        return join_path(self.vcf_out_dir,
                         self.dataset_name+"_genotyped.vcf")

    def __load_jobs_info(self, jobs_setup_file):
        self.__meta_data = {}
        stream = file(jobs_setup_file, "r")
        document = yaml.safe_load(stream)
        # load metadata
        for item_key in document:
            if item_key != JOBS_SETUP_SAMPLES_SECTION:
                self.__meta_data[item_key] = document[item_key]
        # load samples information
        document_samples = document[JOBS_SETUP_SAMPLES_SECTION]
        self.__samples = {}
        for sample in document_samples:
            sample_rec = SampleRecord(sample, self.output_dir)
            self.__samples[sample_rec.sample_name] = sample_rec

    def __parse_metadata(self, meta_string):
        match = METADATA_PATTERN.match(meta_string)
        return match.group('key'), match.group('val')

    def __create_directories(self):
        self.create_dir(self.output_dir)
        self.create_dir(self.vcf_out_dir)
        self.create_dir(self.bam_out_dir)
        self.create_dir(self.slurm_log_dir)
        self.create_dir(self.samples_working_dir)
        for sample_name in self.samples:
            sample_dir = join_path(self.samples_working_dir,
                                   sample_name)
            self.create_dir(sample_dir)

    def __garbage_collecting(self):
        for job_name in self.job_dict:
            job_rec = self.job_dict[job_name]
            if job_name.endswith("_base_recal") and (job_rec.job_status == JOB_STATUS_COMPLETED):
                sample_name = job_name.strip("_base_recal")
                sample_rec = self.__samples[sample_name]
                mylogger.info("deleting " + sample_rec.raw_aligned_reads_file)
                self.delete_file(sample_rec.raw_aligned_reads_file)
                mylogger.info("deleting " + sample_rec.sorted_reads_file)
                self.delete_file(sample_rec.sorted_reads_file)
                mylogger.info("deleting " + sample_rec.dedup_reads_file)
                self.delete_file(sample_rec.dedup_reads_file)
                mylogger.info("deleting " + sample_rec.realigned_reads_file)
                self.delete_file(sample_rec.realigned_reads_file)

    def monitor_action(self):
        JobManager.monitor_action(self)
        self.__garbage_collecting()

    def bwa_mem(self,
                sample_name,
                ):
        """ To generate a SAM file containing aligned reads """
    
        sample_rec = self.__samples[sample_name]
        slurm_log_file = join_path(self.slurm_log_dir,
                                   sample_name)
        slurm_log_file += "_bwa_mem_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = sample_name + "_bwa_mem"
        job_script = BWA_MEM_SCRIPT
        job_params = " -I " + sample_rec.fastq_pairs[0]['R1'] 
        job_params += "," + sample_rec.fastq_pairs[0]['R2']
        job_params += " -g " + sample_rec.sample_group
        job_params += " -R " + self.reference
        job_params += " -o " + sample_rec.raw_aligned_reads_file
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "1",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=sample_rec.usage_mail,
                        )
        return job_name

    def sort_sam(self,
                 sample_name,
                 sam_file=None,
                 prereq=None,
                 ):
        """ Use picard command to sort the SAM file and convert it to BAM """
    
        sample_rec = self.__samples[sample_name]
        slurm_log_file = join_path(self.slurm_log_dir,
                                   sample_name)
        slurm_log_file += "_sort_sam_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = sample_name + "_sort_sam"
        job_script = SORT_SAM_SCRIPT
        if sam_file is None:
            job_params = " -I " + sample_rec.raw_aligned_reads_file 
        else:
            job_params = " -I " + sam_file 
        job_params += " -o " + sample_rec.sorted_reads_file
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "2",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=sample_rec.usage_mail,
                        prereq=prereq,
                        )
        return job_name

    def mark_dup(self,
                 sample_name,
                 bam_file=None,
                 prereq=None,
                 ):
        """ Use picard command to mark duplicates """
    
        sample_rec = self.__samples[sample_name]
        slurm_log_file = join_path(self.slurm_log_dir,
                                   sample_name)
        slurm_log_file += "_mark_dup_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = sample_name + "_mark_dup"
        job_script = MARK_DUP_SCRIPT
        if bam_file is None:
            job_params = " -I " + sample_rec.sorted_reads_file
        else:
            job_params = " -I " + bam_file
        job_params += " -o " + sample_rec.dedup_reads_file
        job_params += " -m " + sample_rec.metrics_file
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "4",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=sample_rec.usage_mail,
                        prereq=prereq,
                        )
        return job_name

    def create_intervals(self,
                         sample_name,
                         bam_file=None,
                         prereq=None,
                         ):
        """ Create a target list of intervals to be realigned """
        sample_rec = self.__samples[sample_name]
        slurm_log_file = join_path(self.slurm_log_dir,
                                   sample_name)
        slurm_log_file += "_create_intervals_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = sample_name + "_create_intervals"
        job_script = INDELS_TARGET_INTERVALS_SCRIPT
        if bam_file is None:
            job_params = " -I " + sample_rec.dedup_reads_file
        else:
            job_params = " -I " + bam_file
        job_params += " -k " + ",".join(self.known_indels)
        job_params += " -R " + self.reference
        job_params += " -o " + sample_rec.indels_target_intervals_file
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "2",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=sample_rec.usage_mail,
                        prereq=prereq,
                        )
        return job_name

    def indels_realign(self,
                       sample_name,
                       bam_file=None,
                       indels_target_intervals_file=None,
                       prereq=None,
                       ):
        """ Local Realignment around Indels """
        sample_rec = self.__samples[sample_name]
        slurm_log_file = join_path(self.slurm_log_dir,
                                   sample_name)
        slurm_log_file += "_indels_realign_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = sample_name + "_indels_realign"
        job_script = INDELS_REALIGN_SCRIPT
        if bam_file is None:
            job_params = " -I " + sample_rec.dedup_reads_file
        else:
            job_params = " -I " + bam_file
        job_params += " -k " + ",".join(self.known_indels)
        job_params += " -R " + self.reference
        if indels_target_intervals_file is None:
            job_params += " -t " + sample_rec.indels_target_intervals_file
        else:
            job_params += " -t " + indels_target_intervals_file
        job_params += " -o " + sample_rec.realigned_reads_file
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "2",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=sample_rec.usage_mail,
                        prereq=prereq,
                        )
        return job_name

    def base_recal(self,
                   sample_name,
                   bam_file=None,
                   prereq=None,
                   ):
        """ Base Quality Score Recalibration """
        sample_rec = self.__samples[sample_name]
        slurm_log_file = join_path(self.slurm_log_dir,
                                   sample_name)
        slurm_log_file += "_base_recal_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = sample_name + "_base_recal"
        job_script = BASE_RECAL_SCRIPT
        if bam_file is None:
            job_params = " -I " + sample_rec.realigned_reads_file
        else:
            job_params = " -I " + bam_file
        job_params += " -R " + self.reference
        job_params += " -k " + ",".join(self.known_indels)
        job_params += "," + self.dbsnp
        job_params += " -t " + sample_rec.recal_table_file
        job_params += " -a " + sample_rec.post_recal_table_file
        job_params += " -p " + sample_rec.recalibration_plots
        job_params += " -o " + sample_rec.recal_reads_file
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "3",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=sample_rec.usage_mail,
                        prereq=prereq,
                        )
        return job_name

    def hap_cal(self,
                sample_name,
                bam_file=None,
                prereq=None,
                ):
        """
        Call SNPs and indels simultaneously via local re-assembly of haplotypes
        in an active region
        """
    
        sample_rec = self.__samples[sample_name]
        slurm_log_file = join_path(self.slurm_log_dir,
                                   sample_name)
        slurm_log_file += "_hap_cal_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = sample_name + "_hap_cal"
        job_script = HAPLOTYPE_CALLER_SCRIPT
        if bam_file is None:
            job_params = " -I " + sample_rec.recal_reads_file
        else:
            job_params = " -I " + bam_file
        job_params += " -R " + self.reference
        if self.targets_interval_list is not None:
            job_params += " -L " + self.targets_interval_list
        if self.dbsnp is not None:
            job_params += " -d " + self.dbsnp
        job_params += " -o " + sample_rec.raw_gvcf_file
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "2",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=sample_rec.usage_mail,
                        prereq=prereq,
                        )
        return job_name

    def combine_gvcfs(self,
                      gvcf_list=None,
                      prereq=None,
                      ):
        """
        Combines any number of gVCF files that were produced by
        the Haplotype Caller into a single joint gVCF file
        """
    
        slurm_log_file = join_path(self.slurm_log_dir,
                                   self.dataset_name)
        slurm_log_file += "_comb_gvcfs_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = self.dataset_name + "_comb_gvcfs"
        job_script = COMBINE_GVCFS_SCRIPT
        gvcf_list_file = join_path(self.samples_working_dir,
                                   job_name+"_gvcf.list")
        f_gvcf = open(gvcf_list_file, "w")
        if gvcf_list is None:
            for sample_name in self.samples:
                sample_rec = self.samples[sample_name]
                f_gvcf.write(sample_rec.raw_gvcf_file + "\n")
        else:
            map(lambda x: f_gvcf.write(x+"\n"), gvcf_list)
        f_gvcf.close()
        job_params = " -G " + gvcf_list_file 
        job_params += " -o " + self.combined_gvcfs_file
        job_params += " -r " + self.reference
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "2",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=self.usage_mail,
                        prereq=prereq,
                        )
        return job_name

    def genotype_gvcfs(self,
                      gvcf_list=None,
                      prereq=None,
                      ):
        """
        Combines any number of gVCF files that were produced by
        the Haplotype Caller into a single joint gVCF file
        """
    
        slurm_log_file = join_path(self.slurm_log_dir,
                                   self.dataset_name)
        slurm_log_file += "_genotype_gvcfs_"
        slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        slurm_log_file += ".log"
        job_name = self.dataset_name + "_gentype_gvcfs"
        job_script = GENOTYPE_GVCFS_SCRIPT
        gvcf_list_file = join_path(self.samples_working_dir,
                                   job_name+"_gvcf.list")
        f_gvcf = open(gvcf_list_file, "w")
        if gvcf_list is None:
            f_gvcf.write(self.combined_gvcfs_file + "\n")
        else:
            map(lambda x: f_gvcf.write(x+"\n"), gvcf_list)
        f_gvcf.close()

        job_params = " -G " + gvcf_list_file 
        job_params += " -o " + self.genotyped_gvcfs_file
        job_params += " -r " + self.reference
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "2",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=self.usage_mail,
                        prereq=prereq,
                        )
        return job_name

    def preprocess_sample(self,
                          sample_name,
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
    
        mylogger.debug("preprocess sample " + sample_name)
        job_name_bwa_mem = self.bwa_mem(sample_name)
        job_name_sort_sam = self.sort_sam(sample_name,
                                          prereq=[job_name_bwa_mem])
        job_name_mark_dup = self.mark_dup(sample_name,
                                          prereq=[job_name_sort_sam])
        job_name_create_intervals = self.create_intervals(sample_name,
                                                          prereq=[job_name_mark_dup])
        job_name_indels_realign = self.indels_realign(sample_name,
                                                      prereq=[job_name_create_intervals])
        job_name_base_recal = self.base_recal(sample_name,
                                              prereq=[job_name_indels_realign])
        job_name_hap_cal = self.hap_cal(sample_name,
                                        prereq=[job_name_base_recal])
        return job_name_hap_cal

    def preprocess_dataset(self):
        """
        An aggregate flow to pre-process a DNASeq sample. The flow consists of
          - bwa-mem alignment
          - reads sorting
          - marking duplicates
          - haplotype caller
        The flow output a gvcf file for one sample. The return value is the last
        job name of the process.
        """

        print "abc"
        prereq = []
        for sample_name in self.samples:
            if self.__samples[sample_name].preprocess_sample:
                prereq.append(self.preprocess_sample(sample_name))
        job_name_combine_gvcfs = self.combine_gvcfs(prereq=prereq) 
        job_name_genotype_gvcfs = self.genotype_gvcfs(prereq=[job_name_combine_gvcfs]) 
        return job_name_genotype_gvcfs

def create_jobs_setup_file(dataset_name,
                           sample_group,
                           project_code,
                           reference_file,
                           project_out_dir,
                           samples_root_dir,
                           jobs_report_file=None,
                           platform='Illumina',
                           known_indels_file=None,
                           dbsnp_file=None,
                           variants_calling="YES",
                           targets_interval_list=None,
                           dataset_usage_mail="NO",
                           sample_usage_mail={},
                           out_jobs_setup_file=None,
                           ):
    mylogger.getLogger(__name__)
    if jobs_report_file is None:
        jobs_report_file = join_path(project_out_dir,
                                     dataset_name+"_rpt.txt")
    if out_jobs_setup_file is None:
        out_jobs_setup_file = join_path(samples_root_dir,
                                       dataset_name+"_job_setup.txt")
    stream = file(out_jobs_setup_file, 'w')
    job_setup_document = {}
    job_setup_document[JOBS_SETUP_DATASET_NAME_KEY] = dataset_name
    job_setup_document[JOBS_SETUP_PROJECT_CODE_KEY] = project_code
    if known_indels_file is not None:
        job_setup_document[JOBS_SETUP_KNOWN_INDELS_KEY] = known_indels_file
    if dbsnp_file is not None:
        job_setup_document[JOBS_SETUP_DBSNP_KEY] = dbsnp_file
    job_setup_document[JOBS_SETUP_REFFERENCE_KEY] = reference_file
    job_setup_document[JOBS_SETUP_OUTPUT_DIR_KEY] = project_out_dir
    job_setup_document[JOBS_SETUP_VARIANTS_CALLING_KEY] = variants_calling
    if targets_interval_list is not None:
        job_setup_document[JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY] = targets_interval_list
    job_setup_document[JOBS_SETUP_JOBS_REPORT_FILE_KEY] = jobs_report_file
    job_setup_document[JOBS_SETUP_DATASET_USAGE_MAIL_KEY] = dataset_usage_mail

    samples = []
    # look for all sub-directory of samples_root_dir
    for item in listdir(samples_root_dir):
        item_path = join_path(samples_root_dir, item)
        if isdir(item_path):
            # look for fastq.gz files to create sample information
            sample = None
            for subdir_item in listdir(item_path):
                if subdir_item.endswith('.fastq.gz'):
                    sample = {}
                    sample[JOBS_SETUP_SAMPLE_NAME_KEY] = item
                    sample[JOBS_SETUP_PLATFORM_KEY] = platform
                    sample[JOBS_SETUP_SAMPLE_GROUP_KEY] = sample_group
                    if item in sample_usage_mail:
                        sample[JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY] = 'YES'
                    else:
                        sample[JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY] = 'NO'
                    sample[JOBS_SETUP_PREPROCESS_SAMPLE_KEY] = 'YES'
                    break
            # look for fastq.gz files to fastq files info for R1,R2 pairs or R1 alone
            fastq_files = []
            for subdir_item in listdir(item_path):
                if not subdir_item.endswith('.fastq.gz'):
                    continue
                if 'R1' not in subdir_item:
                    continue
                pair = {}
                R1_file = join_path(item_path, subdir_item)
                pair[JOBS_SETUP_R1_KEY] = R1_file
                R2_file = R1_file.replace('R1', 'R2')
                if isfile(R2_file):
                    pair[JOBS_SETUP_R2_KEY] = R2_file
                fastq_files.append(pair)
            # look for fastq.gz files to fastq files info for single file (without 'R1')
            for subdir_item in listdir(item_path):
                if not subdir_item.endswith('.fastq.gz'):
                    continue
                if 'R1' in subdir_item:
                    continue
                pair = {}
                file_path = join_path(item_path, subdir_item)
                if 'R2' in file_path and isfile(file_path.replace('R2', 'R1')):
                    continue
                pair[JOBS_SETUP_R1_KEY] = file_path
                fastq_files.append(pair)
            if sample is not None:
                sample[JOBS_SETUP_FASTQ_PAIRS_KEY] = fastq_files
                samples.append(sample)
    job_setup_document[JOBS_SETUP_SAMPLES_SECTION] = samples
    pyaml.dump(job_setup_document, stream)
