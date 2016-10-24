import yaml
import pyaml
from os.path import join as join_path
from collections import OrderedDict
from pycmm.settings import DFLT_JOB_ALLOC_TIME
from pycmm.utils import get_dict_val
from pycmm.utils.jobman import JobManager
from pycmm.cmmlib.familylib import params_to_yaml_doc
from pycmm.cmmlib.familylib import SamplesInfo
from pycmm.cmmlib.familylib import Family
from pycmm.cmmlib.familylib import JOBS_SETUP_SAMPLES_INFOS_KEY
from pycmm.cmmlib.familylib import NO_FAMILY

JOBS_SETUP_JOBS_REPORT_FILE_KEY = "JOBS_REPORT_FILE"
JOBS_SETUP_PROJECT_NAME_KEY = "PROJECT_NAME"
JOBS_SETUP_PROJECT_OUT_DIR_KEY = "PROJECT_OUT_DIR"
JOBS_SETUP_PROJECT_CODE_KEY = "PROJECT_CODE"
JOBS_SETUP_JOB_ALLOC_TIME_KEY = "JOB_ALLOC_TIME"


class CMMPipeline(JobManager):
    """ A basic pipeline control """

    def __init__(self,
                 jobs_setup_file,
                 family_template=Family,
                 **kwargs
                 ):
        self.__load_jobs_info(jobs_setup_file)
        self.__init_properties(family_template)
        kwargs['jobs_report_file'] = self.jobs_report_file
        super(CMMPipeline, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["project name"] = self.project_name
        raw_repr["project out dir"] = self.project_out_dir
        raw_repr["data out dir"] = self.data_out_dir
        if self.families_info is None:
            samples_info = None
        elif ((len(self.families_info) == 1) and
              (NO_FAMILY in self.families_info)
              ):
            samples_info = self.samples_id
        else:
            samples_info = OrderedDict()
            for fam_id in self.families_info:
                fam_info = self.families_info[fam_id]
                fam_info_dict = {}
                members_id = map(lambda x: x.sample_id,
                                 fam_info.members)
                samples_info["family "+fam_id+" - members"] = members_id
        raw_repr['sample information'] = samples_info
        raw_repr["reports out dir"] = self.rpts_out_dir
        raw_repr["project code"] = self.project_code
        raw_repr["job allocation time"] = self.job_alloc_time
        raw_repr["jobs report file"] = self.jobs_report_file
        return raw_repr

    def __load_jobs_info(self, jobs_setup_file):
        self.__jobs_setup_file = jobs_setup_file
        stream = file(jobs_setup_file, "r")
        self.__jobs_info = yaml.safe_load(stream)

    def __init_properties(self, family_template):
        self.__project_out_dir = None
        self.__samples_info = SamplesInfo(self._get_job_config(JOBS_SETUP_SAMPLES_INFOS_KEY),
                                          family_template=family_template)

    def _get_sub_project_dir(self, sub_dir):
        attr_name = "_"
        attr_name += self.__class__.__name__
        attr_name += "__"
        attr_name += sub_dir
        attr_name += "_dir"
        if not hasattr(self, attr_name):
            full_sub_dir = join_path(self.project_out_dir,
                                     sub_dir)
            setattr(self, attr_name, full_sub_dir)
            self.create_dir(full_sub_dir)
        return getattr(self, attr_name)

    def _get_job_config(self, key, required=False, default_val=None):
        return get_dict_val(self.__jobs_info,
                            key,
                            required=required,
                            default_val=default_val,
                            )

    def run_slurm_monitor_pipeline(self,
                                   class_slurm_bin,
                                   alloc_time,
                                   log_file=None,
                                   ):
        job_name = self.project_name + "_mgr"
        job_script = class_slurm_bin
        job_params = " -j " + self.jobs_setup_file
        if log_file is not None:
            job_params += " -l " + log_file
        self._submit_slurm_job(job_name,
                               "1",
                               job_script,
                               job_params,
                               alloc_time=alloc_time,
                               )

    def run_offline_pipeline(self):
        txt = " "
        txt += self.__class__.__name__
        txt += " doesn't support running offline"
        txt += " "
        self.info()
        self.info()
        self.info()
        self.info(txt.center(140, "!"))
        self.info()

    def _submit_slurm_job(self,
                          job_name,
                          ntasks,
                          job_script,
                          job_params,
                          alloc_time=None,
                          email=False,
                          prereq=None,
                          ):
        if alloc_time is None:
            alloc_time = self.job_alloc_time
        slurm_log_prefix = join_path(self.slurm_log_dir,
                                     job_name)
        self.submit_job(job_name,
                        self.project_code,
                        "core",
                        ntasks,
                        alloc_time,
                        slurm_log_prefix,
                        job_script,
                        job_params,
                        email=email,
                        prereq=prereq,
                        )

    def get_params(self):
        params = OrderedDict()
        params['project name'] = self.project_name
        params['project code'] = self.project_code
        if self.families_info is None:
            samples_info = None
        elif ((len(self.families_info) == 1) and
              (NO_FAMILY in self.families_info)
              ):
            samples_info = self.samples_id
        else:
            samples_info = OrderedDict()
            for fam_id in self.families_info:
                fam_info = self.families_info[fam_id]
                fam_info_txt = "fam_id: " + fam_id
                fam_info_txt += ", members: "
                fam_info_txt += ",".join(map(lambda x: x.sample_id,
                                             fam_info.members))
                samples_info[fam_id] = fam_info_txt
        params['sample information'] = samples_info
        params['job allocation time'] = self.job_alloc_time
        params['project output directory'] = self.project_out_dir
        return params

    @property
    def jobs_setup_file(self):
        return self.__jobs_setup_file

    @property
    def jobs_report_file(self):
        return join_path(self.project_out_dir,
                         self.project_name+"_rpt.txt")

    @property
    def dataset_name(self):
        return self.project_name

    @property
    def project_name(self):
        return self._get_job_config(JOBS_SETUP_PROJECT_NAME_KEY, required=True)

    @property
    def project_out_dir(self):
        if self.__project_out_dir is None:
            self.__project_out_dir = self._get_job_config(JOBS_SETUP_PROJECT_OUT_DIR_KEY,
                                                          required=True)
            self.create_dir(self.__project_out_dir)
        return self.__project_out_dir

    @property
    def slurm_log_dir(self):
        return self._get_sub_project_dir("slurm_log")

    @property
    def data_out_dir(self):
        return self._get_sub_project_dir("data_out")

    @property
    def rpts_out_dir(self):
        return self._get_sub_project_dir("rpts")

# *********************************************************************************************** Need refactoring ***********************************************************************************************
# need to replace it with scratch diretory somehow    
    @property
    def working_dir(self):
        return self._get_sub_project_dir("tmp")
# *********************************************************************************************** Need refactoring ***********************************************************************************************

    @property
    def project_code(self):
        return self._get_job_config(JOBS_SETUP_PROJECT_CODE_KEY)

    @property
    def job_alloc_time(self):
        return self._get_job_config(JOBS_SETUP_JOB_ALLOC_TIME_KEY)

    @property
    def families_info(self):
        return self.__samples_info.families

    @property
    def samples_dict(self):
        return self.__samples_info.samples_dict

    @property
    def samples_id(self):
        return self.__samples_info.samples_id

    @property
    def samples_id_w_fam_pref(self):
        return self.__samples_info.samples_id_w_fam_pref

    def monitor_init(self, **kwargs):
        super(CMMPipeline, self).monitor_init(**kwargs)

    def monitor_action(self, **kwargs):
        super(CMMPipeline, self).monitor_action(**kwargs)

    def monitor_finalize(self, **kwargs):
        super(CMMPipeline, self).monitor_finalize(**kwargs)

# a note on parameters of init_jobs_setup_file and create_jobs_setup_file
# The reason that I explicitly list them is for documentation purpose. It is to
# see the list of parameters that can be passed to these functions
def get_func_arg(var_name, jobs_setup_kwargs, default_val=None):
    if ((var_name in jobs_setup_kwargs) and
        (jobs_setup_kwargs[var_name] is not None)):
        return jobs_setup_kwargs[var_name]
    return default_val
    
def init_jobs_setup_file(*args, **kwargs):
    project_out_dir = kwargs['project_out_dir']
    project_name = kwargs['project_name']
    jobs_report_file = get_func_arg('jobs_report_file',
                                    kwargs,
                                    join_path(project_out_dir,
                                              project_name+"_rpt.txt"),
                                    )
    out_jobs_setup_file = get_func_arg('out_jobs_setup_file',
                                       kwargs,
                                       join_path(project_out_dir,
                                                 project_name+"_jobs_setup.txt"),
                                       )
    stream = file(out_jobs_setup_file, 'w')
    job_setup_document = {}
    job_setup_document[JOBS_SETUP_PROJECT_NAME_KEY] = kwargs['project_name']
    job_setup_document[JOBS_SETUP_PROJECT_OUT_DIR_KEY] = kwargs['project_out_dir']
    project_code = get_func_arg('project_code', kwargs)
    if project_code is not None:
        job_setup_document[JOBS_SETUP_PROJECT_CODE_KEY] = project_code
        job_alloc_time = get_func_arg('job_alloc_time',
                                      kwargs,
                                      DFLT_JOB_ALLOC_TIME,
                                      )
        job_setup_document[JOBS_SETUP_JOB_ALLOC_TIME_KEY] = '"' + job_alloc_time + '"'
    yaml = params_to_yaml_doc(sample_info=get_func_arg('sample_info', kwargs))
    for key in yaml:
        job_setup_document[key] = yaml[key]
    jobs_report_file = get_func_arg('jobs_report_file', kwargs)
    job_setup_document[JOBS_SETUP_JOBS_REPORT_FILE_KEY] = jobs_report_file
    return job_setup_document, stream

def create_jobs_setup_file(*args, **kwargs):
    job_setup_document, stream = init_jobs_setup_file(*args, **kwargs)
    pyaml.dump(job_setup_document, stream)
