import sys
import re
import datetime
from os.path import join as join_path
from pycmm.template import pyCMMBase
from pycmm.utils import mylogger
from pycmm.utils.jobman import JobManager
from pycmm.settings import GATK_ALLOC_TIME

METADATA_PATTERN = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')

BWA_MEM_SCRIPT = "$PYCMM/bash/GATKBP_bwa_mem.sh"
SORT_SAM_SCRIPT = "$PYCMM/bash/GATKBP_sort_sam.sh"
MARK_DUP_SCRIPT = "$PYCMM/bash/GATKBP_mark_dup.sh"
HAPLOTYPE_CALLER_SCRIPT = "$PYCMM/bash/GATKBP_haplotype_caller.sh"
COMBINE_GVCFS_SCRIPT = "$PYCMM/bash/GATKBP_combine_gvcfs.sh"
GENOTYPE_GVCFS_SCRIPT = "$PYCMM/bash/GATKBP_genotype_gvcfs.sh"

JOBS_SETUP_DATASET_NAME_KEY = "DATASET_NAME"
JOBS_SETUP_PROJECT_CODE_KEY = "PROJECT_CODE"
JOBS_SETUP_REFFERENCE_KEY = "REFERENCE"
JOBS_SETUP_OUTPUT_DIR_KEY = "OUTPUT_DIR"
JOBS_SETUP_VARIANTS_CALLING_KEY = "VARIANTS_CALLING"
JOBS_SETUP_JOBS_REPORT_FILE_KEY = "JOBS_REPORT_FILE"
JOBS_SETUP_USAGE_MAIL_KEY = "USAGE_MAIL"

SAMPLE_RECORD_JOBS_NAME_IDX = 0
SAMPLE_RECORD_FASTQ1_IDX = 1
SAMPLE_RECORD_FASTQ2_IDX = 2
SAMPLE_RECORD_SAMPLE_GROUP_IDX = 3
SAMPLE_RECORD_PLATFORM_IDX = 4
SAMPLE_RECORD_LIBRARY_IDX = 5
SAMPLE_RECORD_UNIT_IDX = 6
SAMPLE_RECORD_USAGE_MAIL_IDX = 7
SAMPLE_RECORD_PREPROCESS_SAMPLE_IDX = 8

class SampleRecord(pyCMMBase):
    """ A structure to parse and keep sample information """

    def __init__(self,
                 sample_line,
                 output_dir,
                 ):
        self.__data = sample_line.strip().split('\t')
        self.__output_dir = output_dir

    def get_raw_repr(self):
        return {"sample name": self.sample_name,
                "fastq1 file": self.fastq1_file,
                "fastq2 file": self.fastq2_file,
                "sample group": self.sample_group,
                "platform": self.platform,
                "library": self.library,
                "unit": self.unit,
                "usage limit": self.usage_limit,
                "pre-process sample": self.preprocess_sample,
                }

    @property
    def sample_name(self):
        return self.__data[SAMPLE_RECORD_JOBS_NAME_IDX]

    @property
    def fastq1_file(self):
        return self.__data[SAMPLE_RECORD_FASTQ1_IDX]

    @property
    def fastq2_file(self):
        return self.__data[SAMPLE_RECORD_FASTQ2_IDX]

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
        return join_path(join_path(self.__output_dir,
                                   'bam'),
                         self.sample_name+"_dedup_reads.bam")

    @property
    def metrics_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_metrics.txt")

    @property
    def raw_gvcf_file(self):
        return join_path(join_path(join_path(self.__output_dir,
                                             'tmp'),
                                   self.sample_name),
                         self.sample_name+"_raw.g.vcf")

    @property
    def sample_group(self):
        return self.__data[SAMPLE_RECORD_SAMPLE_GROUP_IDX]

    @property
    def platform(self):
        return self.__data[SAMPLE_RECORD_PLATFORM_IDX]

    @property
    def library(self):
        return self.__data[SAMPLE_RECORD_LIBRARY_IDX]

    @property
    def unit(self):
        return self.__data[SAMPLE_RECORD_UNIT_IDX]

    @property
    def usage_mail(self):
        return self.__data[SAMPLE_RECORD_USAGE_MAIL_IDX] == "YES"

    @property
    def preprocess_sample(self):
        return self.__data[SAMPLE_RECORD_PREPROCESS_SAMPLE_IDX] == "YES"


class GATKBPPipeline(JobManager):
    """ A class to control GATK best practice pipeline """

    def __init__(self,
                 jobs_setup_file,
                 ):
        self.__load_jobs_info(jobs_setup_file)
        self.__meta_data[JOBS_SETUP_VARIANTS_CALLING_KEY] == "YES"
        self.__meta_data[JOBS_SETUP_USAGE_MAIL_KEY] == "YES"
        JobManager.__init__(self,
                            job_report_file=self.jobs_report_file)
        self.__vcf_out_dir = join_path(self.output_dir,
                                       "vcf")
        self.__bam_out_dir = join_path(self.output_dir,
                                       "bam")
        self.__samples_working_dir = join_path(self.output_dir,
                                               "tmp")
        self.__slurm_log_dir = join_path(self.output_dir,
                                         "slurm_log")
        self.__time_stamp = datetime.datetime.now()
        self.__create_directories()

    def get_raw_repr(self):
        return {"dataset name": self.dataset_name,
                "project code": self.project_code,
                "reference": self.reference,
                "output directory": self.output_dir,
                "varaints calling": self.variants_calling,
                "jobs report file": self.jobs_report_file,
                "usage_mail": self.usage_mail,
                }

    @property
    def dataset_name(self):
        return self.__meta_data[JOBS_SETUP_DATASET_NAME_KEY]

    @property
    def project_code(self):
        return self.__meta_data[JOBS_SETUP_PROJECT_CODE_KEY]

    @property
    def reference(self):
        return self.__meta_data[JOBS_SETUP_REFFERENCE_KEY]

    @property
    def output_dir(self):
        return self.__meta_data[JOBS_SETUP_OUTPUT_DIR_KEY]

    @property
    def variants_calling(self):
        return self.__meta_data[JOBS_SETUP_VARIANTS_CALLING_KEY] == "YES"

    @property
    def jobs_report_file(self):
        return self.__meta_data[JOBS_SETUP_JOBS_REPORT_FILE_KEY]

    @property
    def usage_mail(self):
        return self.__meta_data[JOBS_SETUP_USAGE_MAIL_KEY] == "YES"

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

    def __load_jobs_info(self, job_setup_file):
        mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
        f_jobs_info = open(job_setup_file, 'r')
        reader = (line.strip() for line in f_jobs_info)
        line = reader.next()
        self.__meta_data = {}
        while line.startswith('##'):
            # parse meta data
            key, val = self.__parse_metadata(line)
            self.__meta_data[key] = val
            line = reader.next()
        self.__samples = {}
        for line in reader:
            # parse sample information
            sample_rec = SampleRecord(line, self.output_dir)
            sample_name = sample_rec.sample_name
            if sample_name in self.__samples:
                raise Exception("Sample " + sample_name + " is duplicated in jobs setup file")
            self.__samples[sample_name] = sample_rec
        f_jobs_info.close()

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

    def __get_raw_aligned_file(self,
                               sample_name,
                               ):

        bam_out_dir = join_path(self.samples_working_dir,
                                sample_name)
        return join_path(bam_out_dir,
                         sample_name+"_raw_aligned.sam")

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
        job_params = " -S " + sample_rec.fastq1_file 
        job_params += "," + sample_rec.fastq2_file
        job_params += " -g " + sample_rec.sample_group
        job_params += " -r " + self.reference
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
                 prerequisite=None,
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
            job_params = " -S " + sample_rec.raw_aligned_reads_file 
        else:
            job_params = " -S " + sam_file 
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
                        prerequisite=prerequisite,
                        )
        return job_name

    def mark_dup(self,
                 sample_name,
                 bam_file=None,
                 prerequisite=None,
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
            job_params = " -B " + sample_rec.sorted_reads_file
        else:
            job_params = " -B " + bam_file
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
                        prerequisite=prerequisite,
                        )
        return job_name

    def hap_cal(self,
                sample_name,
                bam_file=None,
                targets_interval_list=None,
                prerequisite=None,
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
            job_params = " -B " + sample_rec.dedup_reads_file
        else:
            job_params = " -B " + bam_file
        job_params += " -o " + sample_rec.raw_gvcf_file
        job_params += " -r " + self.reference
        if targets_interval_list is not None:
            job_params += " -L " + targets_interval_list
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "2",
                        GATK_ALLOC_TIME,
                        slurm_log_file,
                        job_script,
                        job_params,
                        email=sample_rec.usage_mail,
                        prerequisite=prerequisite,
                        )
        return job_name

    def combine_gvcfs(self,
                      gvcf_list=None,
                      prerequisite=None,
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
                        prerequisite=prerequisite,
                        )
        return job_name

    def genotype_gvcfs(self,
                      gvcf_list=None,
                      prerequisite=None,
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
                        prerequisite=prerequisite,
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
    
        job_name_bwa_mem = self.bwa_mem(sample_name)
        job_name_sort_sam = self.sort_sam(sample_name,
                                          prerequisite=[job_name_bwa_mem])
        job_name_mark_dup = self.mark_dup(sample_name,
                                          prerequisite=[job_name_sort_sam])
        job_name_hap_cal = self.hap_cal(sample_name,
                                        prerequisite=[job_name_mark_dup])
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

        prerequisite = []
        for sample_name in self.samples:
            prerequisite.append(self.preprocess_sample(sample_name))
        job_name_combine_gvcfs = self.combine_gvcfs(prerequisite=prerequisite) 
        job_name_genotype_gvcfs = self.genotype_gvcfs(prerequisite=[job_name_combine_gvcfs]) 
        return job_name_genotype_gvcfs
