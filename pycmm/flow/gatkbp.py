import pyaml
from os import listdir
from os.path import join as join_path
from os.path import isdir
from os.path import isfile
from collections import OrderedDict
from pycmm.utils.ver import VersionManager
from pycmm.cmmlib import CMMParams
from pycmm.cmmlib.samplelib import JOBS_SETUP_SAMPLES_INFOS_KEY
from pycmm.cmmlib.samplelib import JOBS_SETUP_SAMPLE_ID_KEY
from pycmm.cmmlib.samplelib import JOBS_SETUP_FAMILY_ID_KEY
from pycmm.cmmlib.samplelib import JOBS_SETUP_MEMBERS_LIST_KEY
from pycmm.cmmlib.samplelib import NO_FAMILY
from pycmm.cmmlib.samplelib import Sample
from pycmm.cmmlib.samplelib import Family
from pycmm.cmmlib.samplelib import SamplesInfo
from pycmm.cmmlib.fastqlib import Fastq
from pycmm.cmmlib.fastqlib import is_r1
from pycmm.cmmlib.fastqlib import is_r2
from pycmm.cmmlib.fastqlib import r1tor2
from pycmm.flow import CMMPipeline
from pycmm.flow import get_func_arg
from pycmm.flow import init_jobs_setup_file

PREPROCESS_SAMPLE_SCRIPT = "$PYCMM/bash/GATKBP_preprocess_sample.sh"
SPLIT_GVCFS_SCRIPT = "$PYCMM/bash/split_gvcfs.sh"
COMBINE_GVCFS_SCRIPT = "$PYCMM/bash/GATKBP_combine_gvcfs.sh"
GENOTYPE_GVCFS_SCRIPT = "$PYCMM/bash/GATKBP_genotype_gvcfs.sh"
CONCAT_VCFS_SCRIPT = "$PYCMM/bash/concat_vcfs.sh"
VQSR_SCRIPT = "$PYCMM/bash/GATKBP_VQSR.sh"

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

DFLT_PLATFORM = "Illumina"


class GATKSample(Sample):
    """  To parse and structure GATK sample """

    def __init__(self, sample_info, *args, **kwargs):
        kwargs["member_info"] = sample_info
        super(GATKSample, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(GATKSample, self).get_raw_obj_str(*args, **kwargs)
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

    @property
    def recal_bam_file(self):
        return join_path(self.bam_out_dir,
                         self.sample_id+".dedup.realigned.recal.bam")

    @property
    def gvcf_file(self):
        return join_path(self.gvcf_out_dir,
                         self.sample_id+".g.vcf.gz")

class GATKBPSamplesInfo(SamplesInfo):
    """  To parse and structure all samples information for mutation report """

    def __init__(self,
                 samples_info,
                 *args,
                 **kwargs
                 ):
        super(GATKBPSamplesInfo, self).__init__(samples_info,
                                                GATKBPFamily,
                                                *args,
                                                **kwargs)
#        self.__parse_families()
#
#    def __parse_families(self):
#        if self.datasets is None:
#            self.__families = None
#            return
#        self.__families = OrderedDict()
#        for dataset in self.datasets.values():
#            for fam_id in dataset.families:
#                family = dataset.families[fam_id]
#                self.__families[fam_id] = family
#
#    @property
#    def datasets(self):
#        return self.items
#
#    @property
#    def families(self):
#        return self.__families

class GATKBPFamily(Family):
    """  To parse and structure GATK dummy family """

    def __init__(self, *args, **kwargs):
        self.__members = None
        super(GATKBPFamily, self).__init__(*args, **kwargs)

    @property
    def _id(self):
        return self.fam_id

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

    def __init__(self, *args, **kwargs):
        super(GATKParams, self).__init__(*args, **kwargs)
        self.__init_properties()

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(GATKParams, self).get_raw_obj_str(*args, **kwargs)
        raw_repr["known indels"] = self.known_indels
        raw_repr["reference"] = self.reference
        raw_repr["dbsnp file"] = self.dbsnp
        raw_repr["variants calling"] = self.variants_calling
        raw_repr["targets interval list"] = self.targets_interval_list
        raw_repr["usage mail"] = self.dataset_usage_mail
        raw_repr["split region file"] = self.split_regions_file
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

    def __init__(self, *args, **kwargs):
        kwargs['info_template'] = GATKBPSamplesInfo
        super(GATKBPPipeline, self).__init__(*args, **kwargs)
        self.__init_properties()

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(GATKBPPipeline, self).get_raw_obj_str(*args, **kwargs)
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
        file_name = join_path(self.working_dir,
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

    @property
    def vqsr_vcf(self):
        return join_path(self.vcf_out_dir,
                         self.dataset_name+".VQSR.vcf")

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
        return job_name_hap_cal

    def __garbage_collecting(self):
        pass

    def monitor_init(self, *args, **kwargs):
        super(GATKBPPipeline, self).monitor_init(*args, **kwargs)
        self.process_dataset()

    def monitor_action(self, *args, **kwargs):
        super(GATKBPPipeline, self).monitor_action(*args, **kwargs)
        self.__garbage_collecting()

    def split_gvcfs(self,
                    gvcf_list=None,
                    prereq=None,
                    ):
        """
        Split gvcf files by region
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

        if region_txt is None:
            job_name = self.dataset_name + "_comb_gvcfs"
        else:
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
        Call samples actual genotypes after combined
        """

        if region_txt is None:
            job_name = self.dataset_name + "_gnt_gvcfs"
        else:
            job_name = self.dataset_name + "_gnt_gvcfs_" + region_txt
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

    def concat_vcfs(self,
                    prereq=None,
                    ):
        """
        concat split vcf
        """

        job_name = self.dataset_name + "_cct_vcfs"
        job_script = CONCAT_VCFS_SCRIPT
        vcf_list_file = join_path(self.working_dir,
                                  job_name+"_vcf.list")
        f_vcf = open(vcf_list_file, "w")
        for region_txt in self.gatk_params.split_regions_txt_list:
            f_vcf.write(self.get_genotyped_gvcfs_file_name(region_txt=region_txt) + "\n")
        f_vcf.close()
        output_file = self.get_genotyped_gvcfs_file_name()

        job_params = " -i " + vcf_list_file
        job_params += " -o " + output_file
        if self.project_code is None:
            self.exec_sh(job_script + job_params)
            return output_file
        else:
            self._submit_slurm_job(job_name,
                                   "1",
                                   job_script,
                                   job_params,
                                   email=self.gatk_params.dataset_usage_mail,
                                   prereq=prereq,
                                   )
            return job_name

    def vqsr(self,
             prereq=None,
             ):
        """
        Variant Quality Score Recalibration
        """

        job_name = self.dataset_name + "_vqsr"
        job_script = VQSR_SCRIPT
        output_file = self.vqsr_vcf

        job_params = " -i " + self.get_genotyped_gvcfs_file_name()
        job_params += " -o " + output_file
        if self.project_code is None:
            self.exec_sh(job_script + job_params)
            return output_file
        else:
            self._submit_slurm_job(job_name,
                                   "2",
                                   job_script,
                                   job_params,
                                   email=self.gatk_params.dataset_usage_mail,
                                   prereq=prereq,
                                   )
            return job_name

    def process_dataset(self):
        """
        An aggregate flow to DNASeq workflow consisting of
          - process each sample
              - bwa-mem alignment
              - reads sorting
              - marking duplicates
              - haplotype caller
          - if needed, split gvcf
          - combine gvcf
          - genotype gvcf
          - if needed, concat vcf
        """

        prereq = []
        for sample_id in self.samples:
            if self.samples[sample_id].preprocess_sample:
                prereq.append(self.preprocess_sample(sample_id))
        if ((self.gatk_params.variants_calling) and
            (self.gatk_params.split_regions_txt_list is not None)
            ):
            job_name_split_gvcfs = self.split_gvcfs(prereq=prereq)
            job_name_genotype_gvcfs_list = []
            for region_txt in self.gatk_params.split_regions_txt_list:
                job_name_combine_gvcfs = self.combine_gvcfs(region_txt=region_txt,
                                                            prereq=[job_name_split_gvcfs])
                job_name_genotype_gvcfs_list.append(self.genotype_gvcfs(region_txt=region_txt,
                                                                        prereq=[job_name_combine_gvcfs]))
            job_name_concat_vcfs = self.concat_vcfs(prereq=job_name_genotype_gvcfs_list) 
            return job_name_concat_vcfs
        elif self.gatk_params.variants_calling:
            job_name_combine_gvcfs = self.combine_gvcfs(prereq=prereq) 
            job_name_genotype_gvcfs = self.genotype_gvcfs(prereq=[job_name_combine_gvcfs]) 
            return job_name_genotype_gvcfs

def create_jobs_setup_file(*args, **kwargs):
    job_setup_document, stream = init_jobs_setup_file(*args, **kwargs)
    gatk_params = {}
    known_indels_file = get_func_arg('known_indels_file', kwargs)
    if known_indels_file is not None:
        gatk_params[JOBS_SETUP_KNOWN_INDELS_KEY] = known_indels_file
    dbsnp_file = get_func_arg('dbsnp_file', kwargs)
    if dbsnp_file is not None:
        gatk_params[JOBS_SETUP_DBSNP_KEY] = dbsnp_file
    project_out_dir = kwargs['project_out_dir']
    gatk_params[JOBS_SETUP_REFFERENCE_KEY] = kwargs['reference_file']
    gatk_params[JOBS_SETUP_VARIANTS_CALLING_KEY] = get_func_arg('variants_calling',
                                                                kwargs,
                                                                default_val=False)
    targets_interval_list = get_func_arg('targets_interval_list', kwargs)
    if targets_interval_list is not None:
        gatk_params[JOBS_SETUP_TARGETS_INTERVAL_LIST_KEY] = targets_interval_list
    split_regions_file = get_func_arg('split_regions_file', kwargs)
    if split_regions_file is not None:
        gatk_params[JOBS_SETUP_SPLIT_REGIONS_FILE_KEY] = split_regions_file
    options = []
    dataset_usage_mail = get_func_arg('dataset_usage_mail',
                                      kwargs,
                                      default_val=False)
    if dataset_usage_mail:
        options.append(JOBS_SETUP_USAGE_MAIL)
    if len(options) > 0:
        gatk_params[JOBS_SETUP_OPTIONS_KEY] = options
    job_setup_document[JOBS_SETUP_GATK_PARAMS_SECTION] = gatk_params

    samples_root_dir = kwargs['samples_root_dir']
    sample_usage_mail = get_func_arg('sample_usage_mail',
                                      kwargs,
                                      default_val={})
    preprocess_sample = get_func_arg('preprocess_sample',
                                      kwargs,
                                      default_val=True)
    samples_dict = {}
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
                    fastq_files = []
                    for sample_item in listdir(sample_path):
                        if not sample_item.endswith('.fastq.gz'):
                            continue
                        file_path = join_path(sample_path, sample_item)
                        if is_r2(file_path):
                            continue
                        if sample_id not in samples_dict.keys():
                            # if sample information was not there, create one
                            sample = {}
                            sample[JOBS_SETUP_SAMPLE_ID_KEY] = '"' + sample_id + '"'
                            sample[JOBS_SETUP_PLATFORM_KEY] = DFLT_PLATFORM
                            sample[JOBS_SETUP_SAMPLE_GROUP_KEY] = sample_group
                            if sample_id in sample_usage_mail:
                                sample[JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY] = 'YES'
                            else:
                                sample[JOBS_SETUP_SAMPLE_USAGE_MAIL_KEY] = 'NO'
                            sample[JOBS_SETUP_PREPROCESS_SAMPLE_KEY] = preprocess_sample
                        # parse infor for R1,R2 pairs or R1 alone
                        if is_r1(file_path):
                            pair = {}
                            pair[JOBS_SETUP_R1_KEY] = file_path
                            R2_file = r1tor2(file_path)
                            if (R2_file is not None and
                                isfile(R2_file)):
                                pair[JOBS_SETUP_R2_KEY] = R2_file
                            fastq_files.append(pair)
                        else:
                            # Here it's supposed to be fastq.gz files but not R1 or R2
                            pair = {}
                            pair[JOBS_SETUP_R1_KEY] = file_path
                            fastq_files.append(pair)
                    if sample is not None:
                        sample[JOBS_SETUP_FASTQ_PAIRS_KEY] = fastq_files
                        # check fastq encoding version
                        fastq = Fastq(fastq_files[0][JOBS_SETUP_R1_KEY])
                        sample[JOBS_SETUP_FASTQ_ENCODING_KEY] = fastq.encoding
                        samples_dict[sample_id] = sample
    # dummy family (no family)
    families = []
    family_info = {}
    family_info[JOBS_SETUP_FAMILY_ID_KEY] = '"' + NO_FAMILY + '"'
    family_info[JOBS_SETUP_MEMBERS_LIST_KEY] = samples_dict.values()
    families.append(family_info)
    job_setup_document[JOBS_SETUP_SAMPLES_INFOS_KEY] = families

    pyaml.dump(job_setup_document, stream)
