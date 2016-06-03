import pyaml
from os import listdir
from os.path import join as join_path
from os.path import isdir
from os.path import isfile
from collections import OrderedDict
from pycmm.cmmlib import CMMParams
from pycmm.cmmlib.familylib import JOBS_SETUP_SAMPLES_INFOS_KEY
from pycmm.cmmlib.familylib import JOBS_SETUP_SAMPLE_ID_KEY
from pycmm.cmmlib.familylib import JOBS_SETUP_FAMILY_ID_KEY
from pycmm.cmmlib.familylib import JOBS_SETUP_MEMBERS_LIST_KEY
from pycmm.cmmlib.familylib import NO_FAMILY
from pycmm.cmmlib.familylib import Sample
from pycmm.cmmlib.familylib import Family
from pycmm.cmmlib.fastqlib import Fastq
from pycmm.flow import CMMPipeline
from pycmm.flow import init_jobs_setup_file
#from pycmm.utils.jobman import JOB_STATUS_COMPLETED
from pycmm.utils.ver import VersionManager

#BWA_MEM_SCRIPT = "$PYCMM/bash/GATKBP_bwa_mem.sh"
#SORT_SAM_SCRIPT = "$PYCMM/bash/GATKBP_sort_sam.sh"
#MARK_DUP_SCRIPT = "$PYCMM/bash/GATKBP_mark_dup.sh"
#INDELS_TARGET_INTERVALS_SCRIPT = "$PYCMM/bash/GATKBP_create_intervals.sh"
#INDELS_REALIGN_SCRIPT = "$PYCMM/bash/GATKBP_indels_realign.sh"
#BASE_RECAL_SCRIPT = "$PYCMM/bash/GATKBP_base_recal.sh"
#HAPLOTYPE_CALLER_SCRIPT = "$PYCMM/bash/GATKBP_haplotype_caller.sh"
PREPROCESS_SAMPLE_SCRIPT = "$PYCMM/bash/GATKBP_preprocess_sample.sh"
SPLIT_GVCFS_SCRIPT = "$PYCMM/bash/split_gvcfs.sh"
COMBINE_GVCFS_SCRIPT = "$PYCMM/bash/GATKBP_combine_gvcfs.sh"
GENOTYPE_GVCFS_SCRIPT = "$PYCMM/bash/GATKBP_genotype_gvcfs.sh"

# *************** gatk parameters section ***************
JOBS_SETUP_GATK_PARAMS_SECTION = "GATK_PARAMS"
JOBS_SETUP_KNOWN_INDELS_KEY = "KNOWN_INDELS"
JOBS_SETUP_DBSNP_KEY = "DBSNP"
JOBS_SETUP_REFFERENCE_KEY = "REFERENCE"
JOBS_SETUP_VARIANTS_CALLING_KEY = "VARIANTS_CALLING"
JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY = "TARGETS_INTERVAL_LIST"
JOBS_SETUP_SPLIT_REGIONS_FILE_KEY = "SPLIT_REGIONS_FILE"
JOBS_SETUP_OPTIONS_KEY = "OPTIONS"
JOBS_SETUP_USAGE_MAIL = "USAGE_MAIL"

# *************** samples information section ***************
JOBS_SETUP_PLATFORM_KEY = "platform"
JOBS_SETUP_SAMPLE_GROUP_KEY = "sample_group"
JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY = "usage_mail"
JOBS_SETUP_PREPROCESS_SAMPLE_KEY = "preprocess_sample"
JOBS_SETUP_R1_KEY = "R1"
JOBS_SETUP_R2_KEY = "R2"
JOBS_SETUP_FASTQ_PAIRS_KEY = "fastq_files"
JOBS_SETUP_FASTQ_ENCODING_KEY = "encoding"


class GATKSample(Sample):
    """  To parse and structure GATK sample """

    def __init__(self, sample_info, **kwargs):
        kwargs["member_info"] = sample_info
        super(GATKSample, self).__init__(**kwargs)

    def get_raw_repr(self, **kwargs):
        raw_repr = super(GATKSample, self).get_raw_repr(**kwargs)
        raw_repr["library"] = self.library
        raw_repr["unit"] = self.unit
        raw_repr["sample group"] = self.sample_group
        raw_repr["fastq pairs"] = self.fastq_pairs
        raw_repr["fastq encoding"] = self.encoding
        raw_repr["pre-process sample"] = self.preprocess_sample
        raw_repr["platform"] = self.platform
        raw_repr["usage mail"] = self.usage_mail
        raw_repr["bam directory"] = self.bam_out_dir
        raw_repr["gvcf directory"] = self.gvcf_out_dir
        return raw_repr

    @property
    def library(self):
        return 'lib1'

    @property
    def unit(self):
        return 'unit1'

    @property
    def fastq_pairs(self):
        return self._get_job_config(JOBS_SETUP_FASTQ_PAIRS_KEY,
                                    required=True)

    @property
    def preprocess_sample(self):
        return self._get_job_config(JOBS_SETUP_PREPROCESS_SAMPLE_KEY,
                                    required=True)

    @property
    def platform(self):
        return self._get_job_config(JOBS_SETUP_PLATFORM_KEY,
                                    required=True)

    @property
    def sample_group(self):
        return self._get_job_config(JOBS_SETUP_SAMPLE_GROUP_KEY,
                                    required=True)

    @property
    def encoding(self):
        return self._get_job_config(JOBS_SETUP_FASTQ_ENCODING_KEY,
                                    required=True)

    @property
    def usage_mail(self):
        return self._get_job_config(JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY,
                                    required=True)

#    @property
#    def working_dir(self):
#        return self.__working_dir
#
#    @working_dir.setter
#    def working_dir(self, value):
#        self.__working_dir = value
#
#    @property
#    def pdf_out_dir(self):
#        return self.__pdf_out_dir
#
#    @pdf_out_dir.setter
#    def pdf_out_dir(self, value):
#        self.__pdf_out_dir = value
#
    @property
    def bam_out_dir(self):
        return self.__bam_out_dir

    @bam_out_dir.setter
    def bam_out_dir(self, value):
        self.__bam_out_dir = value

    @property
    def gvcf_out_dir(self):
        return self.__gvcf_out_dir

    @gvcf_out_dir.setter
    def gvcf_out_dir(self, value):
        self.__gvcf_out_dir = value

#    @property
#    def raw_aligned_reads_file(self):
#        return join_path(self.working_dir,
#                         self.sample_id+"_raw_aligned_reads.sam")
#
#    @property
#    def sorted_reads_file(self):
#        return join_path(self.working_dir,
#                         self.sample_id+"_sorted_reads.bam")
#
#    @property
#    def dedup_reads_file(self):
#        return join_path(self.working_dir,
#                         self.sample_id+"_dedup_reads.bam")
#
#    @property
#    def metrics_file(self):
#        return join_path(self.working_dir,
#                         self.sample_id+"_metrics.txt")
#
#    @property
#    def indels_target_intervals_file(self):
#        return join_path(self.working_dir,
#                         self.sample_id+"_indels_target_intervals.list")
#
#    @property
#    def realigned_reads_file(self):
#        return join_path(self.working_dir,
#                         self.sample_id+"_realigned_reads.bam")
#
#    @property
#    def recal_table_file(self):
#        return join_path(self.working_dir,
#                         self.sample_id+"_recal_data.table")
#
#    @property
#    def post_recal_table_file(self):
#        return join_path(self.working_dir,
#                         self.sample_id+"_post_recal_data.table")
#
    @property
    def recal_bam_file(self):
        return join_path(self.bam_out_dir,
                         self.sample_id+".dedup.realigned.recal.bam")

#    @property
#    def recalibration_plots(self):
#        return join_path(self.pdf_out_dir,
#                         self.sample_id+"_recalibration_plots.pdf")
#
    @property
    def gvcf_file(self):
        return join_path(self.gvcf_out_dir,
                         self.sample_id+".g.vcf.gz")

class GATKFamily(Family):
    """  To parse and structure GATK dummy family """

    def __init__(self, **kwargs):
        self.__members = None
        super(GATKFamily, self).__init__(**kwargs)

    @property
    def members(self):
        if self.__members is None:
            self.__members = map(lambda x: GATKSample(sample_info=x),
                                 self._get_job_config(JOBS_SETUP_MEMBERS_LIST_KEY,
                                                      required=True)
                                 )
        return self.__members

class GATKParams(CMMParams):
    """  To handle and parse GATK parameters  """

    def __init__(self, **kwargs):
        super(GATKParams, self).__init__(**kwargs)
        self.__init_properties()

    def get_raw_repr(self, **kwargs):
        raw_repr = super(GATKParams, self).get_raw_repr(**kwargs)
        raw_repr["known indels"] = self.known_indels
        raw_repr["reference"] = self.reference
        raw_repr["dbsnp file"] = self.dbsnp
        raw_repr["variants calling"] = self.variants_calling
        raw_repr["targets interval list"] = self.targets_interval_list
        raw_repr["usage mail"] = self.dataset_usage_mail
        raw_repr["samples"] = self.samples
        return raw_repr

    def __init_properties(self):
        regions_list = self.split_regions_file
        if regions_list is not None:
            regions_list = map(lambda x: x.strip().replace(":", "_").replace("-", "_"),
                               open(regions_list, 'rt'))
        self.__split_regions_txt_list = regions_list

    @property
    def known_indels(self):
        return self._get_job_config(JOBS_SETUP_KNOWN_INDELS_KEY)

    @property
    def dbsnp(self):
        return self._get_job_config(JOBS_SETUP_DBSNP_KEY)

    @property
    def reference(self):
        return self._get_job_config(JOBS_SETUP_REFFERENCE_KEY, required=True)

    @property
    def variants_calling(self):
        return self._get_job_config(JOBS_SETUP_VARIANTS_CALLING_KEY)

    @property
    def targets_interval_list(self):
        return self._get_job_config(JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY)

    @property
    def split_regions_file(self):
        return self._get_job_config(JOBS_SETUP_SPLIT_REGIONS_FILE_KEY)

    @property
    def split_regions_txt_list(self):
        return self.__split_regions_txt_list

    @property
    def options(self):
        return self._get_job_config(JOBS_SETUP_OPTIONS_KEY, default_val=[])
    
    @property
    def dataset_usage_mail(self):
        return JOBS_SETUP_USAGE_MAIL in self.options

class GATKBPPipeline(CMMPipeline):
    """ A class to control GATK best practice pipeline """

    def __init__(self, **kwargs):
        kwargs['family_template'] = GATKFamily
        super(GATKBPPipeline, self).__init__(**kwargs)
        self.__init_properties()

    def get_raw_repr(self, **kwargs):
        raw_repr = super(GATKBPPipeline, self).get_raw_repr(**kwargs)
        raw_repr["dataset name"] = self.dataset_name
        raw_repr["GATK parameters"] = self.gatk_params
        raw_repr["samples"] = self.samples
        return raw_repr

    def __init_properties(self):
        for sample_id in self.samples:
            sample = self.samples[sample_id]
            sample.bam_out_dir = self.bam_out_dir
            sample.gvcf_out_dir = self.gvcf_out_dir

    @property
    def gatk_params(self):
        return GATKParams(entries=self._get_job_config(JOBS_SETUP_GATK_PARAMS_SECTION,
                                                       required=True))

    @property
    def samples(self):
        return self.samples_dict

    @property
    def bam_out_dir(self):
        return self._get_sub_project_dir("bam")

#    @property
#    def pdf_out_dir(self):
#        return self._get_sub_project_dir("pdf")
#
    @property
    def gvcf_out_dir(self):
        return self._get_sub_project_dir("gvcf")

    @property
    def split_gvcfs_out_dir(self):
        return join_path(self.working_dir,
                         "split_gvcfs")

    @property
    def vcf_out_dir(self):
        return self._get_sub_project_dir("vcf")

    def get_combined_gvcfs_file_name(self, region_txt=None):
        file_name = join_path(self.working_dir,
                              self.dataset_name)
        if region_txt is not None:
            file_name += "_" + region_txt
        file_name += "_cohort.g.vcf.gz"
        return file_name

    def get_genotyped_gvcfs_file_name(self, region_txt=None):
        file_name = join_path(self.vcf_out_dir,
                              self.dataset_name)
        if region_txt is not None:
            file_name += "_" + region_txt
        file_name += "_genotyped.vcf.gz"
        return file_name

    def get_third_party_software_version(self):
        vm = VersionManager()
        versions = OrderedDict()
        versions['GATK'] = vm.gatk_version
        return versions

#    def bwa_mem(self,
#                sample_id,
#                ):
#        """ To generate a SAM file containing aligned reads """
#
#        sample_rec = self.samples[sample_id]
#        job_script = BWA_MEM_SCRIPT
#        job_params = " -I " + sample_rec.fastq_pairs[0]['R1'] 
#        job_params += "," + sample_rec.fastq_pairs[0]['R2']
#        job_params += " -g " + sample_rec.sample_group
#        job_params += " -R " + self.gatk_params.reference
#        job_params += " -n " + sample_id
#        job_params += " -o " + sample_rec.raw_aligned_reads_file
#        if self.project_code is None:
#            self.exec_sh(job_script + job_params)
#            return sample_rec.raw_aligned_reads_file
#        else:
#            job_name = sample_id + "_bwa_mem"
#            self._submit_slurm_job(job_name,
#                                   "1",
#                                   job_script,
#                                   job_params,
#                                   email=sample_rec.usage_mail,
#                                   )
#            return job_name
#
#    def sort_sam(self,
#                 sample_id,
#                 sam_file=None,
#                 prereq=None,
#                 ):
#        """ Use picard command to sort the SAM file and convert it to BAM """
#
#        sample_rec = self.samples[sample_id]
#        job_script = SORT_SAM_SCRIPT
#        if sam_file is None:
#            job_params = " -I " + sample_rec.raw_aligned_reads_file 
#        else:
#            job_params = " -I " + sam_file 
#        job_params += " -o " + sample_rec.sorted_reads_file
#        if self.project_code is None:
#            self.exec_sh(job_script + job_params)
#            return sample_rec.sorted_reads_file
#        else:
#            job_name = sample_id + "_sort_sam"
#            self._submit_slurm_job(job_name,
#                                   "2",
#                                   job_script,
#                                   job_params,
#                                   email=sample_rec.usage_mail,
#                                   prereq=prereq,
#                                   )
#            return job_name
#
#    def mark_dup(self,
#                 sample_id,
#                 bam_file=None,
#                 prereq=None,
#                 ):
#        """ Use picard command to mark duplicates """
#
#        sample_rec = self.samples[sample_id]
#        job_script = MARK_DUP_SCRIPT
#        if bam_file is None:
#            job_params = " -I " + sample_rec.sorted_reads_file
#        else:
#            job_params = " -I " + bam_file
#        job_params += " -o " + sample_rec.dedup_reads_file
#        job_params += " -m " + sample_rec.metrics_file
#        if self.project_code is None:
#            self.exec_sh(job_script + job_params)
#            return sample_rec.dedup_reads_file
#        else:
#            job_name = sample_id + "_mark_dup"
#            self._submit_slurm_job(job_name,
#                                   "4",
#                                   job_script,
#                                   job_params,
#                                   email=sample_rec.usage_mail,
#                                   prereq=prereq,
#                                   )
#            return job_name
#
#    def create_intervals(self,
#                         sample_id,
#                         bam_file=None,
#                         prereq=None,
#                         ):
#        """ Create a target list of intervals to be realigned """
#        sample_rec = self.samples[sample_id]
#        job_script = INDELS_TARGET_INTERVALS_SCRIPT
#        if bam_file is None:
#            job_params = " -I " + sample_rec.dedup_reads_file
#        else:
#            job_params = " -I " + bam_file
#        if ((self.gatk_params.known_indels is not None) and
#            (len(self.gatk_params.known_indels) > 0)
#            ):
#            job_params += " -k " + ",".join(self.gatk_params.known_indels)
#        job_params += " -R " + self.gatk_params.reference
#        job_params += " -o " + sample_rec.indels_target_intervals_file
#        if self.project_code is None:
#            self.exec_sh(job_script + job_params)
#            return sample_rec.indels_target_intervals_file
#        else:
#            job_name = sample_id + "_create_intervals"
#            self._submit_slurm_job(job_name,
#                                   "2",
#                                   job_script,
#                                   job_params,
#                                   email=sample_rec.usage_mail,
#                                   prereq=prereq,
#                                   )
#            return job_name
#
#    def indels_realign(self,
#                       sample_id,
#                       bam_file=None,
#                       indels_target_intervals_file=None,
#                       prereq=None,
#                       ):
#        """ Local Realignment around Indels """
#        sample_rec = self.samples[sample_id]
#        job_script = INDELS_REALIGN_SCRIPT
#        if bam_file is None:
#            job_params = " -I " + sample_rec.dedup_reads_file
#        else:
#            job_params = " -I " + bam_file
#        if ((self.gatk_params.known_indels is not None) and
#            (len(self.gatk_params.known_indels) > 0)
#            ):
#            job_params += " -k " + ",".join(self.gatk_params.known_indels)
#        job_params += " -R " + self.gatk_params.reference
#        if indels_target_intervals_file is None:
#            job_params += " -t " + sample_rec.indels_target_intervals_file
#        else:
#            job_params += " -t " + indels_target_intervals_file
#        job_params += " -o " + sample_rec.realigned_reads_file
#        if self.project_code is None:
#            self.exec_sh(job_script + job_params)
#            return sample_rec.realigned_reads_file
#        else:
#            job_name = sample_id + "_indels_realign"
#            self._submit_slurm_job(job_name,
#                                   "2",
#                                   job_script,
#                                   job_params,
#                                   email=sample_rec.usage_mail,
#                                   prereq=prereq,
#                                   )
#            return job_name
#
#    def base_recal(self,
#                   sample_id,
#                   bam_file=None,
#                   prereq=None,
#                   ):
#        """ Base Quality Score Recalibration """
#        sample_rec = self.samples[sample_id]
#        job_script = BASE_RECAL_SCRIPT
#        if bam_file is None:
#            job_params = " -I " + sample_rec.realigned_reads_file
#        else:
#            job_params = " -I " + bam_file
#        job_params += " -R " + self.gatk_params.reference
#        if ((self.gatk_params.known_indels is not None) and
#            (len(self.gatk_params.known_indels) > 0)
#            ):
#            job_params += " -k " + ",".join(self.gatk_params.known_indels)
#        job_params += "," + self.gatk_params.dbsnp
#        job_params += " -t " + sample_rec.recal_table_file
#        job_params += " -a " + sample_rec.post_recal_table_file
#        job_params += " -p " + sample_rec.recalibration_plots
#        job_params += " -o " + sample_rec.recal_bam_file
#        if self.project_code is None:
#            self.exec_sh(job_script + job_params)
#            return sample_rec.recal_bam_file
#        else:
#            job_name = sample_id + "_base_recal"
#            self._submit_slurm_job(job_name,
#                                   "3",
#                                   job_script,
#                                   job_params,
#                                   email=sample_rec.usage_mail,
#                                   prereq=prereq,
#                                   )
#            return job_name
#
#    def hap_cal(self,
#                sample_id,
#                bam_file=None,
#                prereq=None,
#                ):
#        """
#        Call SNPs and indels simultaneously via local re-assembly of haplotypes
#        in an active region
#        """
#
#        sample_rec = self.samples[sample_id]
#        job_script = HAPLOTYPE_CALLER_SCRIPT
#        if bam_file is None:
#            job_params = " -I " + sample_rec.recal_bam_file
#        else:
#            job_params = " -I " + bam_file
#        job_params += " -R " + self.gatk_params.reference
#        if self.gatk_params.targets_interval_list is not None:
#            job_params += " -L " + self.gatk_params.targets_interval_list
#        if self.gatk_params.dbsnp is not None:
#            job_params += " -d " + self.gatk_params.dbsnp
#        job_params += " -o " + sample_rec.gvcf_file
#        if self.project_code is None:
#            self.exec_sh(job_script + job_params)
#            return sample_rec.gvcf_file
#        else:
#            job_name = sample_id + "_hap_cal"
#            self._submit_slurm_job(job_name,
#                                   "2",
#                                   job_script,
#                                   job_params,
#                                   email=sample_rec.usage_mail,
#                                   prereq=prereq,
#                                   )
#            return job_name

    def preprocess_sample(self,
                          sample_id,
                          ):
        """
        An aggregate flow to pre-process a DNASeq sample. The flow consists of
          - bwa-mem alignment
          - reads sorting
          - marking duplicates
          - indel realignment
          - base recalibration
          - haplotype caller
        The flow output a gvcf file for one sample. The return value is the last
        job name of the process.
        """

        self.info("preprocess sample " + sample_id)
        sample_rec = self.samples[sample_id]
        job_script = PREPROCESS_SAMPLE_SCRIPT
        job_params = " -I " + sample_rec.fastq_pairs[0]['R1'] 
        job_params += "," + sample_rec.fastq_pairs[0]['R2']
        job_params += " -g " + sample_rec.sample_group
        job_params += " -R " + self.gatk_params.reference
        job_params += " -n " + sample_id
        job_params += " -k " + ",".join(self.gatk_params.known_indels)
        job_params += " -d " + self.gatk_params.dbsnp
        if self.gatk_params.targets_interval_list is not None:
            job_params += " -L " + self.gatk_params.targets_interval_list
        job_params += " -b " + sample_rec.recal_bam_file
        job_params += " -o " + sample_rec.gvcf_file
        if self.project_code is None:
            self.exec_sh(job_script + job_params)
            return sample_rec.raw_aligned_reads_file
        else:
            job_name = sample_id + "_gvcf"
            self._submit_slurm_job(job_name,
                                   "4",
                                   job_script,
                                   job_params,
                                   email=sample_rec.usage_mail,
                                   )
            return job_name
#        job_name_bwa_mem = self.bwa_mem(sample_id)
#        job_name_sort_sam = self.sort_sam(sample_id,
#                                          prereq=[job_name_bwa_mem])
#        job_name_mark_dup = self.mark_dup(sample_id,
#                                          prereq=[job_name_sort_sam])
#        job_name_create_intervals = self.create_intervals(sample_id,
#                                                          prereq=[job_name_mark_dup])
#        job_name_indels_realign = self.indels_realign(sample_id,
#                                                      prereq=[job_name_create_intervals])
#        job_name_base_recal = self.base_recal(sample_id,
#                                              prereq=[job_name_indels_realign])
#        job_name_hap_cal = self.hap_cal(sample_id,
#                                        prereq=[job_name_base_recal])
        return job_name_hap_cal

    def __garbage_collecting(self):
        pass
#        for job_name in self.job_dict:
#            job_rec = self.job_dict[job_name]
#            if ((job_name.endswith("_base_recal")) and
#                (job_rec.job_status == JOB_STATUS_COMPLETED)
#                ):
#                sample_id = job_name.strip("_base_recal")
#                sample_rec = self.samples[sample_id]
#                self.info("deleting " + sample_rec.raw_aligned_reads_file)
#                self.delete_file(sample_rec.raw_aligned_reads_file)
#                self.info("deleting " + sample_rec.sorted_reads_file)
#                self.delete_file(sample_rec.sorted_reads_file)
#                self.info("deleting " + sample_rec.dedup_reads_file)
#                self.delete_file(sample_rec.dedup_reads_file)
#                self.info("deleting " + sample_rec.realigned_reads_file)
#                self.delete_file(sample_rec.realigned_reads_file)

    def monitor_init(self, **kwargs):
        super(GATKBPPipeline, self).monitor_init(**kwargs)
        self.preprocess_dataset()

    def monitor_action(self, **kwargs):
        super(GATKBPPipeline, self).monitor_action(**kwargs)
        self.__garbage_collecting()

    def split_gvcfs(self,
                    gvcf_list=None,
                    prereq=None,
                    ):
        """
        Combines any number of gVCF files that were produced by
        the Haplotype Caller into a single joint gVCF file
        """

        job_name = self.dataset_name + "_split_gvcfs"
        job_script = SPLIT_GVCFS_SCRIPT
        job_params = " -i " + self.gvcf_out_dir
        job_params += " -r " + self.gatk_params.split_regions_file
        job_params += " -o " + self.split_gvcfs_out_dir
        if self.project_code is None:
            self.exec_sh(job_script + job_params)
            return None
        else:
            self._submit_slurm_job(job_name,
                                   "1",
                                   job_script,
                                   job_params,
                                   email=self.gatk_params.dataset_usage_mail,
                                   prereq=prereq,
                                   )
            return job_name

    def combine_gvcfs(self,
                      gvcf_list=None,
                      prereq=None,
                      region_txt=None,
                      ):
        """
        Combines any number of gVCF files that were produced by
        the Haplotype Caller into a single joint gVCF file
        """

        job_name = self.dataset_name + "_comb_gvcfs_" + region_txt
        job_script = COMBINE_GVCFS_SCRIPT
        gvcf_list_file = join_path(self.working_dir,
                                   job_name+"_gvcf.list")
        f_gvcf = open(gvcf_list_file, "w")
        if (gvcf_list is None) and (region_txt is None):
            for sample_id in self.samples:
                sample_rec = self.samples[sample_id]
                f_gvcf.write(sample_rec.gvcf_file + "\n")
            output_file = self.get_combined_gvcfs_file_name()
        elif gvcf_list is None:
            for sample_id in self.samples:
                split_gz = join_path(join_path(self.split_gvcfs_out_dir,
                                               region_txt),
                                     sample_id+"."+region_txt+".g.vcf.gz")
                f_gvcf.write(split_gz + "\n")
            output_file = self.get_combined_gvcfs_file_name(region_txt=region_txt)
        else:
            map(lambda x: f_gvcf.write(x+"\n"), gvcf_list)
            output_file = self.get_combined_gvcfs_file_name()
        f_gvcf.close()
        job_params = " -G " + gvcf_list_file 
        job_params += " -o " + output_file
        job_params += " -r " + self.gatk_params.reference
        if self.project_code is None:
            self.exec_sh(job_script + job_params)
            return output_file
        else:
            self._submit_slurm_job(job_name,
                                   "6",
                                   job_script,
                                   job_params,
                                   email=self.gatk_params.dataset_usage_mail,
                                   prereq=prereq,
                                   )
            return job_name

    def genotype_gvcfs(self,
                       gvcf_list=None,
                       prereq=None,
                       region_txt=None,
                      ):
        """
        Combines any number of gVCF files that were produced by
        the Haplotype Caller into a single joint gVCF file
        """

        job_name = self.dataset_name + "_genotype_gvcfs" + region_txt
        job_script = GENOTYPE_GVCFS_SCRIPT
        gvcf_list_file = join_path(self.working_dir,
                                   job_name+"_gvcf.list")
        f_gvcf = open(gvcf_list_file, "w")
        if (gvcf_list is None) and (region_txt is None):
            f_gvcf.write(self.get_combined_gvcfs_file_name() + "\n")
            output_file = self.get_genotyped_gvcfs_file_name()
        elif gvcf_list is None:
            f_gvcf.write(self.get_combined_gvcfs_file_name(region_txt=region_txt) + "\n")
            output_file = self.get_genotyped_gvcfs_file_name(region_txt=region_txt)
        else:
            map(lambda x: f_gvcf.write(x+"\n"), gvcf_list)
            output_file = self.get_genotyped_gvcfs_file_name()
        f_gvcf.close()

        job_params = " -G " + gvcf_list_file
        job_params += " -o " + output_file
        job_params += " -r " + self.gatk_params.reference
        if self.project_code is None:
            self.exec_sh(job_script + job_params)
            return output_file
        else:
            self._submit_slurm_job(job_name,
                                   "4",
                                   job_script,
                                   job_params,
                                   email=self.gatk_params.dataset_usage_mail,
                                   prereq=prereq,
                                   )
            return job_name

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

        prereq = []
        for sample_id in self.samples:
            if self.samples[sample_id].preprocess_sample:
                prereq.append(self.preprocess_sample(sample_id))
        if ((self.gatk_params.variants_calling) and
            (self.gatk_params.split_regions_txt_list is not None)
            ):
            job_name_split_gvcfs = self.split_gvcfs(prereq=prereq)
            for region_txt in self.gatk_params.split_regions_txt_list:
                job_name_combine_gvcfs = self.combine_gvcfs(region_txt=region_txt,
                                                            prereq=[job_name_split_gvcfs])
                job_name_genotype_gvcfs = self.genotype_gvcfs(region_txt=region_txt,
                                                              prereq=[job_name_combine_gvcfs]) 
        elif self.gatk_params.variants_calling:
            job_name_combine_gvcfs = self.combine_gvcfs(prereq=prereq) 
            job_name_genotype_gvcfs = self.genotype_gvcfs(prereq=[job_name_combine_gvcfs]) 
            return job_name_genotype_gvcfs

def create_jobs_setup_file(project_name,
                           project_out_dir,
                           reference_file,
                           samples_root_dir,
                           platform='Illumina',
                           known_indels_file=None,
                           dbsnp_file=None,
                           variants_calling=False,
                           targets_interval_list=None,
                           split_regions_file=None,
                           dataset_usage_mail=False,
                           sample_usage_mail={},
                           preprocess_sample=True,
                           project_code=None,
                           job_alloc_time=None,
                           jobs_report_file=None,
                           out_jobs_setup_file=None,
                           ):
    job_setup_document, stream = init_jobs_setup_file(project_name=project_name,
                                                      project_out_dir=project_out_dir,
                                                      project_code=project_code,
                                                      job_alloc_time=job_alloc_time,
                                                      sample_info=None,
                                                      jobs_report_file=jobs_report_file,
                                                      out_jobs_setup_file=out_jobs_setup_file,
                                                      )
    gatk_params = {}
    if known_indels_file is not None:
        gatk_params[JOBS_SETUP_KNOWN_INDELS_KEY] = known_indels_file
    if dbsnp_file is not None:
        gatk_params[JOBS_SETUP_DBSNP_KEY] = dbsnp_file
    gatk_params[JOBS_SETUP_REFFERENCE_KEY] = reference_file
    gatk_params[JOBS_SETUP_VARIANTS_CALLING_KEY] = variants_calling
    if targets_interval_list is not None:
        gatk_params[JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY] = targets_interval_list
    if split_regions_file is not None:
        gatk_params[JOBS_SETUP_SPLIT_REGIONS_FILE_KEY] = split_regions_file
    options = []
    if dataset_usage_mail:
        options.append(JOBS_SETUP_USAGE_MAIL)
    if len(options) > 0:
        gatk_params[JOBS_SETUP_OPTIONS_KEY] = options
    job_setup_document[JOBS_SETUP_GATK_PARAMS_SECTION] = gatk_params

    samples = []
    # look for all sub-directories of samples_root_dir for sample groups
    for item in listdir(samples_root_dir):
        item_path = join_path(samples_root_dir, item)
        if isdir(item_path):
            # if it is a directory, lets assume that it's a sample group
            sample_group = item
            sample_group_path = item_path
            for subitem in listdir(sample_group_path):
                subitem_path = join_path(sample_group_path, subitem)
                if isdir(subitem_path):
                    # if it is a directory, lets assume that it's a sample directory
                    sample_id = subitem
                    sample_path = subitem_path
                    # look for fastq.gz files to create sample information
                    sample = None
                    for sample_item in listdir(sample_path):
                        if sample_item.endswith('.fastq.gz'):
                            sample = {}
                            sample[JOBS_SETUP_SAMPLE_ID_KEY] = '"' + sample_id + '"'
                            sample[JOBS_SETUP_PLATFORM_KEY] = platform
                            sample[JOBS_SETUP_SAMPLE_GROUP_KEY] = sample_group
                            if sample_id in sample_usage_mail:
                                sample[JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY] = 'YES'
                            else:
                                sample[JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY] = 'NO'
                            sample[JOBS_SETUP_PREPROCESS_SAMPLE_KEY] = preprocess_sample
                            break
                    # look for fastq.gz files to fastq files info for R1,R2 pairs or R1 alone
                    fastq_files = []
                    for sample_item in listdir(sample_path):
                        if not sample_item.endswith('.fastq.gz'):
                            continue
                        if 'R1' not in sample_item:
                            continue
                        pair = {}
                        R1_file = join_path(sample_path, sample_item)
                        pair[JOBS_SETUP_R1_KEY] = R1_file
                        R2_file = R1_file.replace('R1', 'R2')
                        if isfile(R2_file):
                            pair[JOBS_SETUP_R2_KEY] = R2_file
                        fastq_files.append(pair)
                    # look for fastq.gz files to fastq files info for single file (without 'R1')
                    for sample_item in listdir(sample_path):
                        if not sample_item.endswith('.fastq.gz'):
                            continue
                        if 'R1' in sample_item:
                            continue
                        pair = {}
                        file_path = join_path(sample_path, sample_item)
                        if 'R2' in file_path and isfile(file_path.replace('R2', 'R1')):
                            continue
                        pair[JOBS_SETUP_R1_KEY] = file_path
                        fastq_files.append(pair)
                    if sample is not None:
                        sample[JOBS_SETUP_FASTQ_PAIRS_KEY] = fastq_files
                        # check fastq encoding version
                        fastq = Fastq(fastq_files[0][JOBS_SETUP_R1_KEY])
                        sample[JOBS_SETUP_FASTQ_ENCODING_KEY] = fastq.encoding
                        samples.append(sample)
    # dummy family (no family)
    families = []
    family_info = {}
    family_info[JOBS_SETUP_FAMILY_ID_KEY] = '"' + NO_FAMILY + '"'
    family_info[JOBS_SETUP_MEMBERS_LIST_KEY] = samples
    families.append(family_info)
    job_setup_document[JOBS_SETUP_SAMPLES_INFOS_KEY] = families

    pyaml.dump(job_setup_document, stream)
