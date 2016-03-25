import yaml
import pyaml
from os.path import join as join_path
from pycmm.settings import DFLT_FLOW_ALLOC_TIME
from pycmm.settings import DFLT_RPT_ALLOC_TIME
from pycmm.utils.jobman import JobManager

JOBS_SETUP_JOBS_REPORT_FILE_KEY = "JOBS_REPORT_FILE"
JOBS_SETUP_PROJECT_NAME_KEY = "PROJECT_NAME"
JOBS_SETUP_PROJECT_OUT_DIR_KEY = "PROJECT_OUT_DIR"
JOBS_SETUP_PROJECT_CODE_KEY = "PROJECT_CODE"
JOBS_SETUP_FLOW_ALLOC_TIME_KEY = "FLOW_ALLOC_TIME"
JOBS_SETUP_RPT_ALLOC_TIME_KEY = "RPT_ALLOC_TIME"


class CMMPipeline(JobManager):
    """ A basic pipeline control """

    def __init__(self,
                 jobs_setup_file,
                 ):
        self.__load_jobs_info(jobs_setup_file)
        self.__project_out_dir = None
        JobManager.__init__(self,
                            jobs_report_file=self.jobs_report_file)
        self.__slurm_log_dir = None
        self.__data_out_dir = None
        self.__rpts_out_dir = None

    def get_raw_repr(self):
        return {"project name": self.project_name,
                "project out dir": self.project_out_dir,
                "data out dir": self.data_out_dir,
                "reports out dir": self.rpts_out_dir,
                "project code": self.project_code,
                "flow allocation time": self.flow_alloc_time,
                "report allocation time": self.rpt_alloc_time,
                "jobs report file": self.jobs_report_file,
                }

    def __load_jobs_info(self, jobs_setup_file):
        stream = file(jobs_setup_file, "r")
        self._jobs_info = yaml.safe_load(stream)

    @property
    def jobs_report_file(self):
        return join_path(self.project_out_dir,
                         self.project_name+"_rpt.txt")

    @property
    def project_name(self):
        return self._jobs_info[JOBS_SETUP_PROJECT_NAME_KEY]

    @property
    def project_out_dir(self):
        if self.__project_out_dir is None:
            self.__project_out_dir = self._jobs_info[JOBS_SETUP_PROJECT_OUT_DIR_KEY]
            self.create_dir(self.__project_out_dir)
        return self.__project_out_dir

    @property
    def slurm_log_dir(self):
        if self.__slurm_log_dir is None:
            self.__slurm_log_dir = join_path(self.project_out_dir,
                                             "slurm_log")
            self.create_dir(self.__slurm_log_dir)
        return self.__slurm_log_dir

    @property
    def data_out_dir(self):
        if self.__data_out_dir is None:
            self.__data_out_dir = join_path(self.project_out_dir,
                                            "data_out")
            self.create_dir(self.__data_out_dir)
        return self.__data_out_dir

    @property
    def rpts_out_dir(self):
        if self.__rpts_out_dir is None:
            self.__rpts_out_dir = join_path(self.project_out_dir,
                                            "rpts")
            self.create_dir(self.__rpts_out_dir)
        return self.__rpts_out_dir

    @property
    def project_code(self):
        if JOBS_SETUP_PROJECT_CODE_KEY in self._jobs_info:
            return self._jobs_info[JOBS_SETUP_PROJECT_CODE_KEY]
        return None

    @property
    def flow_alloc_time(self):
        if JOBS_SETUP_FLOW_ALLOC_TIME_KEY in self._jobs_info:
            return self._jobs_info[JOBS_SETUP_FLOW_ALLOC_TIME_KEY]
        return None

    @property
    def rpt_alloc_time(self):
        if JOBS_SETUP_RPT_ALLOC_TIME_KEY in self._jobs_info:
            return self._jobs_info[JOBS_SETUP_RPT_ALLOC_TIME_KEY]
        return None

    def monitor_init(self):
        JobManager.monitor_init(self)

    def monitor_action(self):
        JobManager.monitor_action(self)

    def monitor_finalize(self):
        JobManager.monitor_finalize(self)

def init_jobs_setup_file(project_name,
                         project_out_dir,
                         project_code=None,
                         flow_alloc_time=None,
                         rpt_alloc_time=None,
                         jobs_report_file=None,
                         out_jobs_setup_file=None,
                         ):
    if jobs_report_file is None:
        jobs_report_file = join_path(project_out_dir,
                                     project_name+"_rpt.txt")
    if out_jobs_setup_file is None:
        out_jobs_setup_file = join_path(project_out_dir,
                                        project_name+"_jobs_setup.txt")
    stream = file(out_jobs_setup_file, 'w')
    job_setup_document = {}
    job_setup_document[JOBS_SETUP_PROJECT_NAME_KEY] = project_name
    job_setup_document[JOBS_SETUP_PROJECT_OUT_DIR_KEY] = project_out_dir
    if project_code is not None:
        job_setup_document[JOBS_SETUP_PROJECT_CODE_KEY] = project_code
        if flow_alloc_time is None:
            flow_alloc_time = DFLT_FLOW_ALLOC_TIME
        job_setup_document[JOBS_SETUP_FLOW_ALLOC_TIME_KEY] = '"' + flow_alloc_time + '"'
        if rpt_alloc_time is None:
            rpt_alloc_time = DFLT_RPT_ALLOC_TIME
        job_setup_document[JOBS_SETUP_RPT_ALLOC_TIME_KEY] = '"' + rpt_alloc_time + '"'
    job_setup_document[JOBS_SETUP_JOBS_REPORT_FILE_KEY] = jobs_report_file
    return job_setup_document, stream

def create_jobs_setup_file(project_name,
                           project_out_dir,
                           project_code=None,
                           flow_alloc_time=None,
                           rpt_alloc_time=None,
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
    pyaml.dump(job_setup_document, stream)
