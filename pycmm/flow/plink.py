# This Python file uses the following encoding: utf-8
import pyaml
from os.path import join as join_path
from pycmm.settings import PLINK_HAP_ASSOCS_SLURM_BIN
from pycmm.template import pyCMMBase
from pycmm.utils import exec_sh
from pycmm.utils.dnalib import ALL_CHROMS
from pycmm.utils.dnalib import DNARegion
from pycmm.flow import CMMPipeline
from pycmm.flow import init_jobs_setup_file
from pycmm.cmmlib.plinklib import merge_hap_assocs
from pycmm.cmmlib.plinklib import merge_lmiss_map

PLINK_DUMMY_SCRIPT = "$PYCMM/bash/plink_dummy.sh"

DFLT_CUTOFF_PVALUE = "0.05"
DFLT_HAP_WINDOW_SIZES = "1"

# *************** plink parameters section ***************
JOBS_SETUP_PLINK_PARAMS_SECTION = "PLINK_PARAMS"
JOBS_SETUP_INPUT_FILE_PREFIX_KEY = "INPUT_FILE_PREFIX"
JOBS_SETUP_INPUT_BINARY_KEY = "INPUT_BINARY"
JOBS_SETUP_INPUT_DNA_REGIONS_KEY = "INPUT_DNA_REGIONS"
JOBS_SETUP_PHENOTYPE_FILE_KEY = "PHENOTYPE_FILE"
JOBS_SETUP_CUTOFF_PVALUE_KEY = "CUTOFF_PVALUE"
JOBS_SETUP_HAP_WINDOW_SIZES_KEY = "HAP_WINDOW_SIZES"
JOBS_SETUP_HAP_ASSOC_FILTER_CRITERIA_KEY = "FILTERING_CRITERIA"
JOBS_SETUP_HAP_ASSOC_FILTER_PVALUE005 = "PVALUE005"
JOBS_SETUP_HAP_ASSOC_FILTER_DISEASE_SNP = "DISEASE_SNP"


class FilterCriteria(pyCMMBase):
    """  To handle and parse PLINK parameters  """

    def __init__(self, criterias, **kwargs):
        self.__criterias = criterias
        super(FilterCriteria, self).__init__(**kwargs)

    @property
    def filter_pvalue005(self):
        return JOBS_SETUP_HAP_ASSOC_FILTER_PVALUE005 in self.__criterias

    @property
    def filter_disease_snp(self):
        return JOBS_SETUP_HAP_ASSOC_FILTER_DISEASE_SNP in self.__criterias

class PlinkParams(pyCMMBase):
    """  To handle and parse PLINK parameters  """

    def __init__(self, params, **kwargs):
        self.__params = params
        super(PlinkParams, self).__init__(**kwargs)

    def get_raw_repr(self):
        return {"input file prefix": self.input_file_prefix,
                "input binary": self.input_binary,
                "phenotype file": self.phenotype_file,
                }

    @property
    def input_file_prefix(self):
        return self.__params[JOBS_SETUP_INPUT_FILE_PREFIX_KEY]

    @property
    def input_binary(self):
        return self.__params[JOBS_SETUP_INPUT_BINARY_KEY]

    @property
    def input_dna_regions(self):
        if JOBS_SETUP_INPUT_DNA_REGIONS_KEY in self.__params:
            input_dna_regions = map(lambda x: str(x),
                              self.__params[JOBS_SETUP_INPUT_DNA_REGIONS_KEY])
        else:
            input_dna_regions = ALL_CHROMS
        return map(lambda x: DNARegion(x), input_dna_regions)

    @property
    def phenotype_file(self):
        return self.__params[JOBS_SETUP_PHENOTYPE_FILE_KEY]

    @property
    def cutoff_pvalue(self):
        if JOBS_SETUP_CUTOFF_PVALUE_KEY in self.__params:
            return str(self.__params[JOBS_SETUP_CUTOFF_PVALUE_KEY])
        return DFLT_CUTOFF_PVALUE

    @property
    def hap_window_sizes(self):
        if JOBS_SETUP_HAP_WINDOW_SIZES_KEY in self.__params:
            return self.__params[JOBS_SETUP_HAP_WINDOW_SIZES_KEY]
        return DFLT_HAP_WINDOW_SIZES

    @property
    def filter_criteria(self):
        if JOBS_SETUP_HAP_ASSOC_FILTER_CRITERIA_KEY in self.__params:
            return FilterCriteria(self.__params[JOBS_SETUP_HAP_ASSOC_FILTER_CRITERIA_KEY])

class PlinkPipeline(CMMPipeline):
    """ To control PLINK execution pipeline """

    def __init__(self, **kwargs):
        self.__plink_params = None
        self.__scratch_file_prefix = None
        self.__hap_assocs_out = None
        super(PlinkPipeline, self).__init__(**kwargs)

    @property
    def plink_params(self):
        if self.__plink_params is None:
            self.__plink_params = PlinkParams(self._jobs_info[JOBS_SETUP_PLINK_PARAMS_SECTION])
        return self.__plink_params

    @property
    def hap_assoc_alloc_time(self):
        return self.flow_alloc_time

    @property
    def plink_out_dir(self):
        return self.data_out_dir

    @property
    def scratch_file_prefix(self):
        if self.__scratch_file_prefix is None:
            self.__scratch_file_prefix = join_path(self.local_scratch_dir,
                                                   self.get_tmp_file_name())
        return self.__scratch_file_prefix

    @property
    def hap_assocs_out(self):
        return self.__hap_assocs_out

    def get_plink_hap_assoc_paramss(self, input_file_prefix=None):
        params = " --noweb"
        if self.plink_params.input_binary:
            params += " --bfile "
        else:
            params += " --file "
        if input_file_prefix is not None:
            params += input_file_prefix
        else:
            params += self.plink_params.input_file_prefix
        if self.plink_params.phenotype_file is not None:
            params += " --pheno " + self.plink_params.phenotype_file
        params += " --hap-assoc"
        params += " --geno 0.1"
        params += " --maf 0.01"
        self.__hap_assocs_out = []
        for region_idx in xrange(len(self.plink_params.input_dna_regions)):
            input_dna_region = self.plink_params.input_dna_regions[region_idx]
            region_params = params + " --chr " + input_dna_region.chrom
            region_key = self.project_name
            region_key += "_chr" + input_dna_region.chrom
            if input_dna_region.start_pos is not None:
                region_key += "_" + input_dna_region.start_pos
                region_params += " --from-bp " + input_dna_region.start_pos
                region_params += " --to-bp " + input_dna_region.end_pos
            self.__hap_assocs_out.append([])
            for window_size in self.plink_params.hap_window_sizes:
                out_session_key = region_key
                out_session_key += "_WIN_" + str(window_size)
                out_file_prefix = join_path(self.plink_out_dir,
                                            out_session_key)
                self.__hap_assocs_out[region_idx].append(out_file_prefix + ".assoc.hap")
                out_params = region_params + " --out " + out_file_prefix
                out_params += " --hap-window " + str(window_size)
                yield out_params, out_session_key

    def __submit_hap_assoc_jobs(self):
        for params, session_key in self.get_plink_hap_assoc_paramss():
            job_name = session_key
            slurm_log_file = join_path(self.slurm_log_dir,
                                       job_name+".log")
            job_script = PLINK_DUMMY_SCRIPT
            self.submit_job(job_name,
                            self.project_code,
                            "core",
                            "1",
                            self.flow_alloc_time,
                            slurm_log_file,
                            job_script,
                            params,
                            )

    def monitor_init(self, **kwargs):
        super(PlinkPipeline, self).monitor_init(**kwargs)
        self.__submit_hap_assoc_jobs()

## *************************************************************** keep this part of code until I'm certain that there is no way to specify nodelist *************************************************************** 
#    def copy_plink_file_to_scratch(self, suffix):
#        src_file = self.plink_params.input_file_prefix + "." + suffix
#        dst_file = self.scratch_file_prefix + "." + suffix
#        self.copy_file(src_file, dst_file)
#
#    def copy_plink_files_to_scratch(self):
#        if self.plink_params.input_binary:
#            self.copy_plink_file_to_scratch("fam")
#            self.copy_plink_file_to_scratch("bed")
#            self.copy_plink_file_to_scratch("bim")
#        else:
#            self.copy_plink_file_to_scratch("ped")
#            self.copy_plink_file_to_scratch("map")
#
#    def monitor_init(self):
#        CMMPipeline.monitor_init(self)
#        # Assume that the analysis will heavily kill IO so the process
#        # need to copy the data to scratch first, then run the analys
#        # in split windows 
#        self.copy_plink_files_to_scratch()
#        for params, out_session_key in self.get_plink_hap_assoc_paramss(input_file_prefix=self.scratch_file_prefix):
#            job_name = out_session_key
#            slurm_log_file = join_path(self.slurm_log_dir,
#                                       job_name+".log")
#            job_script = PLINK_DUMMY_SCRIPT
#            self.submit_job(job_name,
#                            self.project_code,
#                            "core",
#                            "1",
#                            self.flow_alloc_time,
#                            slurm_log_file,
#                            job_script,
#                            params,
#                            nodelist="m85",
##                            nodelist=self.job_nodelist,
#                            )
## *************************************************************** keep this part of code until I'm certain that there is no way to specify nodelist *************************************************************** 

    def run_hap_assocs_slurm(self,
                             jobs_setup_file,
                             log_file=None,
                             ):
        job_name = self.project_name + "_hap_assocs"
        slurm_log_file = join_path(self.slurm_log_dir,
                                   job_name+".log")
        job_script = PLINK_HAP_ASSOCS_SLURM_BIN
        job_params = " -j " + jobs_setup_file
        if log_file is not None:
            job_params += " -l " + log_file
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        "1",
                        self.flow_alloc_time,
                        slurm_log_file,
                        job_script,
                        job_params,
                        )

    def extract_snps_stat(self, input_dna_region, input_file_prefix=None):
        params = " --noweb"
        if self.plink_params.input_binary:
            params += " --bfile "
        else:
            params += " --file "
        if input_file_prefix is not None:
            params += input_file_prefix
        else:
            params += self.plink_params.input_file_prefix
        params += " --missing"
        if self.plink_params.phenotype_file is not None:
            params += " --within " + self.plink_params.phenotype_file
        params += " --chr " + input_dna_region.chrom
        session_key = self.project_name
        session_key += "_chr" + input_dna_region.chrom
        if input_dna_region.start_pos is not None:
            session_key += "_" + input_dna_region.start_pos
            params += " --from-bp " + input_dna_region.start_pos
            params += " --to-bp " + input_dna_region.end_pos
        out_file_prefix = join_path(self.plink_out_dir,
                                    session_key)
        params += " --out " + out_file_prefix
        cmd = PLINK_DUMMY_SCRIPT + params
        exec_sh(cmd)
        return out_file_prefix + ".lmiss"

    def extract_snps_pos(self, input_dna_region, input_file_prefix=None):
        params = " --noweb"
        if self.plink_params.input_binary:
            params += " --bfile "
        else:
            params += " --file "
        if input_file_prefix is not None:
            params += input_file_prefix
        else:
            params += self.plink_params.input_file_prefix
        params += " --recode"
        params += " --tab"
        params += " --chr " + input_dna_region.chrom
        session_key = self.project_name
        session_key += "_chr" + input_dna_region.chrom
        if input_dna_region.start_pos is not None:
            session_key += "_" + input_dna_region.start_pos
            params += " --from-bp " + input_dna_region.start_pos
            params += " --to-bp " + input_dna_region.end_pos
        out_file_prefix = join_path(self.plink_out_dir,
                                    session_key)
        params += " --out " + out_file_prefix
        cmd = PLINK_DUMMY_SCRIPT + params
        exec_sh(cmd)
        return out_file_prefix + ".map"

    def merge_hap_assocs(self,
                         hap_assoc_files,
                         out_file,
                         locus_prefixs=None,
                         ):
        merge_hap_assocs(hap_assoc_files,
                         out_file,
                         locus_prefixs,
                         filter_criteria=self.plink_params.filter_criteria,
                         )
        self.info(" >>>> merged file: " + out_file)

    def run_hap_assocs_offline(self):
        for input_dna_region in self.plink_params.input_dna_regions:
            lmiss_file = self.extract_snps_stat(input_dna_region)
            map_file = self.extract_snps_pos(input_dna_region)
            out_file = join_path(self.plink_out_dir,
                                 self.project_name)
            out_file += "_chr" + input_dna_region.chrom
            if input_dna_region.start_pos is not None:
                out_file += "_" + input_dna_region.start_pos
            out_file += ".snp.info"
            merge_lmiss_map(lmiss_file, map_file, out_file)
            self.info(" >>>> snp info file: " + out_file)
        for params, session_key in self.get_plink_hap_assoc_paramss():
            cmd = PLINK_DUMMY_SCRIPT + params
            exec_sh(cmd)
        return self.hap_assocs_out

    def __garbage_collecting(self):
        pass

    def monitor_action(self, **kwargs):
        super(PlinkPipeline, self).monitor_action(**kwargs)
        self.__garbage_collecting()

def create_jobs_setup_file(project_name,
                           project_out_dir,
                           input_file_prefix,
                           phenotype_file,
                           input_binary=True,
                           input_dna_regions=None,
                           cutoff_pvalue=None,
                           hap_window_sizes=None,
                           project_code=None,
                           flow_alloc_time=None,
                           rpt_alloc_time=None,
                           filter_criteria=None,
                           jobs_report_file=None,
                           out_jobs_setup_file=None,
                           ):
    job_setup_document, stream = init_jobs_setup_file(project_name=project_name,
                                                      project_out_dir=project_out_dir,
                                                      project_code=project_code,
                                                      flow_alloc_time=flow_alloc_time,
                                                      rpt_alloc_time=rpt_alloc_time,
                                                      jobs_report_file=jobs_report_file,
                                                      out_jobs_setup_file=out_jobs_setup_file,
                                                      )
    plink_params = {}
    plink_params[JOBS_SETUP_INPUT_FILE_PREFIX_KEY] = input_file_prefix
    plink_params[JOBS_SETUP_INPUT_BINARY_KEY] = input_binary
    if input_dna_regions is not None:
        plink_params[JOBS_SETUP_INPUT_DNA_REGIONS_KEY] = input_dna_regions.split(",")
    plink_params[JOBS_SETUP_PHENOTYPE_FILE_KEY] = phenotype_file
    if cutoff_pvalue is not None:
        plink_params[JOBS_SETUP_CUTOFF_PVALUE_KEY] = cutoff_pvalue
    if hap_window_sizes is not None:
        plink_params[JOBS_SETUP_HAP_WINDOW_SIZES_KEY] = hap_window_sizes.split(",")
    if filter_criteria is not None:
        criteria_list = []
        for criteria in filter_criteria.split(","):
            if criteria == JOBS_SETUP_HAP_ASSOC_FILTER_PVALUE005:
                criteria_list.append(JOBS_SETUP_HAP_ASSOC_FILTER_PVALUE005)
            if criteria == JOBS_SETUP_HAP_ASSOC_FILTER_DISEASE_SNP:
                criteria_list.append(JOBS_SETUP_HAP_ASSOC_FILTER_DISEASE_SNP)
        plink_params[JOBS_SETUP_HAP_ASSOC_FILTER_CRITERIA_KEY] = criteria_list
    job_setup_document[JOBS_SETUP_PLINK_PARAMS_SECTION] = plink_params

    pyaml.dump(job_setup_document, stream)
