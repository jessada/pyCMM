import time
import sys
import datetime
from pycmm.template import pyCMMBase
from pycmm.utils import exec_sh

JOB_STATUS_PENDING = "PENDING"
JOB_STATUS_RUNNING = "RUNNING"
JOB_STATUS_COMPLETED = "COMPLETED"
JOB_STATUS_CANCELLED = "CANCELLED"
JOB_STATUS_FAILED = "FAILED"

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
        self.job_status = "NA"
        self.email = None
        self.prereq = None
    
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
                "job status": self.job_status,
                "report usage email": self.email,
                "pre-requisite": self.prereq,
                }

class JobManager(pyCMMBase):
    """ A class to manage UPPMAX SLURM job """

    def __init__(self,
                 jobs_report_file=None,
                 ):
        self.job_dict = {}
        self.__job_rpt_fmt = "{job_id}"
        self.__job_rpt_fmt += "\t{partition}"
        self.__job_rpt_fmt += "\t{job_name}"
        self.__job_rpt_fmt += "\t{project_code}"
        self.__job_rpt_fmt += "\t{alloc_time}"
        self.__job_rpt_fmt += "\t{cpus}"
        self.__job_rpt_fmt += "\t{usage_mail}"
        self.__job_rpt_fmt += "\t{dependency}"
        self.__job_rpt_fmt += "\t{job_status}"
        self.__job_rpt_file = jobs_report_file

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
        if job_rec.email:
            cmd += " -C usage_mail"
        cmd += " -o " + job_rec.slurm_log_file
        if (job_rec.prereq is not None) and (len(job_rec.prereq) > 0) and (type(job_rec.prereq) is list):
            cmd += " --dependency=afterok"
            self.debug(job_rec.prereq)
            for job_name in job_rec.prereq:
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
                   email=False,
                   prereq=None,
                   ):
#        self.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
        job_rec = JobRecord()
        job_rec.job_name = job_name
        job_rec.project_code = project_code
        job_rec.partition_type = partition_type
        job_rec.ntasks = ntasks
        job_rec.alloc_time = alloc_time
        job_rec.slurm_log_file = slurm_log_file
        job_rec.job_script = job_script
        job_rec.job_params = job_params
        job_rec.email = email
        job_rec.prereq =  prereq
        cmd = self.__get_sbatch_cmd(job_rec)
        out = self.__exec_sh(cmd)
        job_rec.job_id = out.strip().split()[-1]
        self.job_dict[job_name] = job_rec
        return None

    def __write_job_report(self):
        f_rpt = open(self.__job_rpt_file, "w")
        f_rpt.write("# Last update: " + str(datetime.datetime.now()) + "\n")
        f_rpt.write(self.__job_rpt_fmt.format(job_id="#JOBID",
                                              partition="PARTITION",
                                              job_name="NAME",
                                              project_code="ACCOUNT",
                                              alloc_time="ALLOC_TIME",
                                              cpus="CPUS",
                                              usage_mail="USAGE_EMAIL",
                                              dependency="DEPENDENCY",
                                              job_status="STATUS",
                                              )+"\n")
        for job_name in self.job_dict:
            job_rec = self.job_dict[job_name]
            f_rpt.write(self.__job_rpt_fmt.format(job_id=job_rec.job_id,
                                                  partition=job_rec.partition_type,
                                                  job_name=job_rec.job_name,
                                                  project_code=job_rec.project_code,
                                                  alloc_time=job_rec.alloc_time,
                                                  cpus=job_rec.ntasks,
                                                  usage_mail=str(job_rec.email),
                                                  dependency=str(job_rec.prereq),
                                                  job_status=job_rec.job_status,
                                                  )+"\n")
        f_rpt.close()

    def monitor_init(self):
        """
        this function will be executed at the begining of
        monitor_jobs process
        """
        # virtual function
        pass

    def monitor_action(self):
        """
        this function will be executed every interval during
        monitor_jobs process
        """
        # virtual function
        if self.__job_rpt_file is not None:
            self.__write_job_report()

    def monitor_finalize(self):
        """
        this function will be executed after
        monitor_jobs process
        """
        # virtual function
        pass
    
    def monitor_jobs(self, interval=10):
        self.monitor_init()
        while True:
            self.update_job_status()
            self.monitor_action()
            if self.all_job_done:
                break
            time.sleep(interval)
        self.monitor_finalize()

    def update_job_status(self):
        for job_name in self.job_dict:
            job_rec = self.job_dict[job_name]
            if job_rec.job_status == JOB_STATUS_COMPLETED:
                continue
            if job_rec.job_status == JOB_STATUS_FAILED:
                continue
            if job_rec.job_status.startswith(JOB_STATUS_CANCELLED):
                continue
            # refresh job status
            job_rec.job_status = self.get_job_status(job_name)

    @property
    def all_job_done(self):
        """
        If there is at least one job running or pending,
        return False (not yet done)
        otherwise return True 
        """
        for job_name in self.job_dict:
            job_rec = self.job_dict[job_name]
            if job_rec.job_status == JOB_STATUS_COMPLETED:
                continue
            if job_rec.job_status == JOB_STATUS_FAILED:
                continue
            if job_rec.job_status.startswith(JOB_STATUS_CANCELLED):
                continue
            return False
        return True

    def get_job_status(self,
                       job_name,
                       ):
#        self.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
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
        return self.job_dict[job_name].job_id

    def cancel_job(self,
                   job_name,
                   ):
        if job_name.isdigit():
            job_id = job_name
        else:
            job_id = self.get_job_id(job_name)
        job_status = self.get_job_status(job_id)
#        self.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
        cmd = "scancel -b " + job_id
        out = self.__exec_sh(cmd)
