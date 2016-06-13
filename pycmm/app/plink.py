import sys
from collections import OrderedDict
from pycmm.settings import PLINK_CREATE_JOB_SETUP_FILE_DESCRIPTION
from pycmm.settings import PLINK_HAP_ASSOCS_DESCRIPTION
from pycmm.settings import PLINK_HAP_ASSOCS_SLURM_DESCRIPTION
from pycmm.settings import PLINK_MERGE_HAP_ASSOCS_DESCRIPTION
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm.utils import set_log_file
from pycmm.flow.plink import PlinkPipeline
from pycmm.flow.plink import create_jobs_setup_file

# *********************************************************************************************** Need refactoring ***********************************************************************************************
# Need to be replaced with the ones form pycmm.app
def display_plink_config(func_name,
                         app_description,
                         kwargs,
                         pl,
                         ):
    log_file = set_log_file(kwargs['log_file'])
    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['jobs setup file (-j)'] = kwargs['jobs_setup_file']
    if 'hap_assoc_files' in kwargs:
        required_params['raw haplotype association study files (--hap_assoc_files)'] = kwargs['hap_assoc_files']
    if 'hap_assoc_file' in kwargs:
        required_params['haplotype association study file (--hap_assoc)'] = kwargs['hap_assoc_file']
    if 'snp_info_file' in kwargs:
        required_params['SNP information file (--snp_info)'] = kwargs['snp_info_file']
    if 'out_file' in kwargs:
        required_params['output file name (--out)'] = kwargs['out_file']
    optional_params = None
    if log_file is not None:
        optional_params = OrderedDict()
        optional_params['log file (-l)'] = log_file
    disp.show_config(app_description=app_description,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.disp_params_set("Pipeline parameters", pl.get_params())
    if hasattr(pl, "plink_params"):
        plink_params = OrderedDict()
        plink_params['input file prefix'] = pl.plink_params.input_file_prefix
        plink_params['input binary'] = pl.plink_params.input_binary
        plink_params['phenotype file'] = pl.plink_params.phenotype_file
        input_dna_regions_param = []
        for input_dna_region in pl.plink_params.input_dna_regions:
            param = "chr" + input_dna_region.chrom
            if input_dna_region.start_pos is not None:
                param += ":" + input_dna_region.start_pos
                param += "-" + input_dna_region.end_pos
            input_dna_regions_param.append(param)
        plink_params['dna regions'] = input_dna_regions_param
        plink_params['haplotype windows'] = pl.plink_params.hap_window_sizes
        if pl.plink_params.filter_criteria is not None:
            filter_criteria = {}
            filter_criteria["p-value 0.05"] = pl.plink_params.filter_criteria.filter_pvalue005
            filter_criteria["disease snps"] = pl.plink_params.filter_criteria.filter_disease_snp
            plink_params['filter criteria'] = filter_criteria
        else:
            plink_params['filter criteria'] = None
        disp.disp_params_set("Plink parameters", plink_params)
    if hasattr(pl, "rpt_params"):
        rpt_params = OrderedDict()
        rpt_params['cutoff p-value'] = pl.rpt_params.cutoff_pvalue
        rpt_params['cutoff odds ratio'] = pl.rpt_params.cutoff_ors
        rpt_params['families haplotype file prefix'] = pl.rpt_params.fam_hap_prefix
        disp.disp_params_set("Report parameters", rpt_params)
    disp.new_section_txt(" . . . E X E C U T I N G . . . ")

def app_pycmm_plink_hap_assocs_slurm(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = PlinkPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_plink_config(func_name,
                         PLINK_HAP_ASSOCS_SLURM_DESCRIPTION,
                         kwargs,
                         pl,
                         )
    # check if itself is in job mode 
    if len(pl.job_nodelist) != 0:
        pl.monitor_jobs()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_plink_hap_assocs(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = PlinkPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_plink_config(func_name,
                         PLINK_HAP_ASSOCS_DESCRIPTION,
                         kwargs,
                         pl,
                         )
    if pl.project_code is not None:
        pl.run_hap_assocs_slurm(jobs_setup_file=kwargs['jobs_setup_file'],
                                log_file=kwargs['log_file'])
    else:
        pl.run_hap_assocs_offline()
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_pycmm_plink_merge_hap_assocs(*args, **kwargs):
    mylogger.getLogger(__name__)

    func_name = sys._getframe().f_code.co_name
    pl = PlinkPipeline(jobs_setup_file=kwargs['jobs_setup_file'])
    mylogger.getLogger(__name__)
    display_plink_config(func_name,
                         PLINK_MERGE_HAP_ASSOCS_DESCRIPTION,
                         kwargs,
                         pl,
                         )
    pl.merge_hap_assocs(hap_assoc_files=kwargs['hap_assoc_files'],
                        out_file=kwargs['out_file']
                        )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
# *********************************************************************************************** Need refactoring ***********************************************************************************************

def app_pycmm_plink_create_jobs_setup_file(*args, **kwargs):
    mylogger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params = OrderedDict()
    required_params['project name (-d)'] = kwargs['project_name']
    required_params['project output directory (-O)'] = kwargs['project_out_dir']
    input_binary = kwargs['input_binary']
    input_file_prefix = kwargs['input_file_prefix']
    if kwargs['input_binary']:
        required_params['input file prefix (--bfile)'] = input_file_prefix
    else:
        required_params['input file prefix (--file)'] = input_file_prefix
    optional_params = OrderedDict()
    if kwargs['input_dna_regions'] is not None:
        optional_params['input dna regions (-R)'] = kwargs['input_dna_regions'].split(",")
    if kwargs['filter_criteria'] is not None:
        optional_params['filter criteria (--filter_criteria)'] = kwargs['filter_criteria'].split(",")
    required_params['phenotype file (--pheno)'] = kwargs['phenotype_file']
    if kwargs['project_code'] is not None:
        optional_params['project code (-p)'] = kwargs['project_code']
    if kwargs['sample_info'] is not None:
        optional_params['sample information (-p)'] = kwargs['sample_info'].split(",")
    optional_params['flow allocation time (--flow_alloc_time)'] = kwargs['flow_alloc_time']
    optional_params['cut-off p-value (--cutoff_pvalue)'] = kwargs['cutoff_pvalue']
    optional_params['cut-off odds ratio (--cutoff_ors)'] = kwargs['cutoff_ors']
    optional_params['families haplotyep file prefix (--fam_hap)'] = kwargs['fam_hap_prefix']
    optional_params['haplotype windows (--hap_window)'] = kwargs['hap_window_sizes'].split(",")
    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=PLINK_CREATE_JOB_SETUP_FILE_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    create_jobs_setup_file(project_name=kwargs['project_name'],
                           project_out_dir=kwargs['project_out_dir'],
                           input_file_prefix=kwargs['input_file_prefix'],
                           phenotype_file=kwargs['phenotype_file'],
                           input_binary=kwargs['input_binary'],
                           input_dna_regions=kwargs['input_dna_regions'],
                           cutoff_pvalue=kwargs['cutoff_pvalue'],
                           cutoff_ors=kwargs['cutoff_ors'],
                           hap_window_sizes=kwargs['hap_window_sizes'],
                           project_code=kwargs['project_code'],
                           sample_info=kwargs['sample_info'],
                           flow_alloc_time=kwargs['flow_alloc_time'],
                           filter_criteria=kwargs['filter_criteria'],
                           fam_hap_prefix=kwargs['fam_hap_prefix'],
                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
                           )
    mylogger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")
