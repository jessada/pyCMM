import time
import sys
from pycmm.template import pyCMMBase
from pycmm.utils import exec_sh
from pycmm.utils import mylogger

JOB_STATUS_PENDING="PENDING"
JOB_STATUS_RUNNING="RUNNING"
JOB_STATUS_COMPLETED="COMPLETED"
JOB_STATUS_FAILED="FAILED"

class JobRecord(pyCMMBase):
    """ to keep UPPMAX SLURM job information """

    def __init__(self):
        self.job_name = None
        self.project_code = None
        self.partition_type = None
        self.ntasks = None
        self.alloc_time = None
        self.slurm_log_file = None
        self.job_script = None
        self.job_params = None
        self.job_id = None
        self.prerequisite = None
    
    def get_raw_repr(self):
        return {"job name": self.job_name,
                "project code": self.project_code,
                "partition type": self.partition_type,
                "number of cpus": self.ntasks,
                "allocation time": self.alloc_time,
                "slurm log file": self.slurm_log_file,
                "job script": self.job_script,
                "job paramters": self.job_params,
                "job id": self.job_id,
                "pre-requisite": self.prerequisite,
                }

class JobManager(pyCMMBase):
    """ A class to manage UPPMAX SLURM job """

    def __init__(self):
        self.__job_dict = {}

    def get_raw_repr(self):
        return {"NA1": "NA1",
                "NA2": "NA2",
                }

    def __exec_sh(self, cmd):
        p, err_code = exec_sh(cmd)
        out, err_msg = p.communicate()
        if err_code:
            raise Exception("[exit " + str(err_code) + "]: " + err_msg)
        return out

    def __get_sbatch_cmd(self, job_rec):
        cmd = "sbatch"
        cmd += " -A " + job_rec.project_code
        cmd += " -p " + job_rec.partition_type
        cmd += " -n " + job_rec.ntasks
        cmd += " -t " + job_rec.alloc_time
        cmd += " -J " + job_rec.job_name
        cmd += " -o " + job_rec.slurm_log_file
        if (job_rec.prerequisite is not None) and (type(job_rec.prerequisite) is list):
            cmd += " --dependency=afterok"
            mylogger.debug(job_rec.prerequisite)
            for job_name in job_rec.prerequisite:
                job_id = self.get_job_id(job_name)
                cmd += ":" + job_id
        cmd += " " + job_rec.job_script
        cmd += " " + job_rec.job_params
        return cmd

    def submit_job(self,
                   job_name,
                   project_code,
                   partition_type,
                   ntasks,
                   alloc_time,
                   slurm_log_file,
                   job_script,
                   job_params,
                   prerequisite=None,
                   ):
        mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
        job_rec = JobRecord()
        job_rec.job_name = job_name
        job_rec.project_code = project_code
        job_rec.partition_type = partition_type
        job_rec.ntasks = ntasks
        job_rec.alloc_time = alloc_time
        job_rec.slurm_log_file = slurm_log_file
        job_rec.job_script = job_script
        job_rec.job_params = job_params
        job_rec.prerequisite =  prerequisite
        cmd = self.__get_sbatch_cmd(job_rec)
        out = self.__exec_sh(cmd)
        job_rec.job_id = out.strip().split()[-1]
        self.__job_dict[job_name] = job_rec
        return None

    def get_job_status(self,
                       job_name,
                       ):
        mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
        if job_name.isdigit():
            job_id = job_name
        else:
            job_id = self.get_job_id(job_name)
        cmd = "sacct -j"
        cmd += " " + job_id
        cmd += " | grep " + job_id
        cmd += " | head -1"
        out = ""
        icount = 0
        while (len(out) == 0) and (icount < 100):
            time.sleep(1)
            icount += 1
            out = self.__exec_sh(cmd)
        if len(out) == 0 :
            raise Exception(job_name + " status cannot be found")
        return out.strip().split()[5]

    def get_job_id(self,
                   job_name,
                   ):
        return self.__job_dict[job_name].job_id

    def cancel_job(self,
                   job_name,
                   ):
        if job_name.isdigit():
            job_id = job_name
        else:
            job_id = self.get_job_id(job_name)
        job_status = self.get_job_status(job_id)
        mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
        cmd = "scancel -b " + job_id
        out = self.__exec_sh(cmd)
