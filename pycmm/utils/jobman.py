import time
import datetime
import sys
from pycmm.template import pyCMMBase

JOB_STATUS_PENDING = "PENDING"
JOB_STATUS_RUNNING = "RUNNING"
JOB_STATUS_COMPLETED = "COMPLETED"
JOB_STATUS_CANCELLED = "CANCELLED"
JOB_STATUS_FAILED = "FAILED"

class JobRecord(pyCMMBase):
    """ to keep UPPMAX SLURM job information """

    def __init__(self, **kwargs):
        self.job_name = None
        self.project_code = None
        self.partition_type = None
        self.ntasks = None
        self.alloc_time = None
        self.__slurm_log_prefix = None
        self.__slurm_log_file = None
        self.job_script = None
        self.job_params = None
        self.job_id = None
        self.job_status = "NA"
        self.email = None
        self.prereq = None
        self.nodelist = None
        super(JobRecord, self).__init__(**kwargs)
    
    def get_raw_repr(self):
        return {"job name": self.job_name,
                "project code": self.project_code,
                "partition type": self.partition_type,
                "number of cpus": self.ntasks,
                "allocation time": self.alloc_time,
                "slurm log file": self.slurm_log_prefix,
                "job script": self.job_script,
                "job paramters": self.job_params,
                "job id": self.job_id,
                "job status": self.job_status,
                "report usage email": self.email,
                "pre-requisite": self.prereq,
                "node list": self.nodelist,
                }

    @property
    def slurm_log_prefix(self):
        return self.__slurm_log_prefix

    @slurm_log_prefix.setter
    def slurm_log_prefix(self, value):
        self.__slurm_log_prefix = value
        self.__slurm_log_file = value
        self.__slurm_log_file += "_"
        self.__slurm_log_file += self.time_stamp.strftime("%Y%m%d%H%M%S")
        self.__slurm_log_file += ".log"

    @property
    def slurm_log_file(self):
        return self.__slurm_log_file

class JobManager(pyCMMBase):
    """ A class to manage UPPMAX SLURM job """

    def __init__(self,
                 jobs_report_file=None,
                 **kwargs
                 ):
        self.__job_dict = {}
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
        self.__job_id = None
        self.__job_name = None
        self.__job_nodelist = None
        super(JobManager, self).__init__(**kwargs)

    @property
    def job_dict(self):
        return self.__job_dict

    @property
    def job_id(self):
        if self.__job_id is None:
            p, stdout_data = self.exec_sh("echo $SLURM_JOB_ID")
            self.__job_id = stdout_data.strip()
        return self.__job_id

    @property
    def job_name(self):
        if self.__job_name is None:
            p, stdout_data = self.exec_sh("echo $SLURM_JOB_NAME")
            self.__job_name = stdout_data.strip()
        return self.__job_name

    @property
    def job_nodelist(self):
        if self.__job_nodelist is None:
            p, stdout_data = self.exec_sh("echo $SLURM_JOB_NODELIST")
            self.__job_nodelist = stdout_data.strip()
        return self.__job_nodelist

    def get_raw_repr(self):
        return None

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
            for job_name in job_rec.prereq:
                job_id = self.get_job_id(job_name)
                cmd += ":" + job_id
        if (job_rec.nodelist is not None) and (len(job_rec.nodelist) > 0):
            cmd += " --nodelist=" + job_rec.nodelist
        cmd += " " + job_rec.job_script
        cmd += " " + job_rec.job_params
        return cmd

    def submit_job(self,
                   job_name,
                   project_code,
                   partition_type,
                   ntasks,
                   alloc_time,
                   slurm_log_prefix,
                   job_script,
                   job_params,
                   email=False,
                   prereq=None,
                   nodelist=None,
                   ):
        job_rec = JobRecord()
        job_rec.job_name = job_name
        job_rec.project_code = project_code
        job_rec.partition_type = partition_type
        job_rec.ntasks = ntasks
        job_rec.alloc_time = alloc_time
        job_rec.slurm_log_prefix = slurm_log_prefix
        job_rec.job_script = job_script
        job_rec.job_params = job_params
        job_rec.email = email
        job_rec.prereq =  prereq
        job_rec.nodelist =  nodelist
        cmd = self.__get_sbatch_cmd(job_rec)
        p, stdout_data = self.exec_sh(cmd)
        job_rec.job_id = stdout_data.strip().split()[-1]
        self.__job_dict[job_name] = job_rec
        return job_rec.job_id

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
        for job_name in self.__job_dict:
            job_rec = self.__job_dict[job_name]
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
    
    def monitor_jobs(self, interval=30):
        self.monitor_init()
        while True:
            self.update_job_status()
            self.monitor_action()
            if self.all_job_done:
                break
            time.sleep(interval)
        self.monitor_finalize()

    def update_job_status(self):
        for job_name in self.__job_dict:
            job_rec = self.__job_dict[job_name]
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
        for job_name in self.__job_dict:
            job_rec = self.__job_dict[job_name]
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
        if job_name.isdigit():
            job_id = job_name
        else:
            job_id = self.get_job_id(job_name)
        cmd = "sacct -j"
        cmd += " " + job_id
        cmd += " | grep " + job_id
        cmd += " | head -1"
        stdout_data = ""
        icount = 0
        while (len(stdout_data) == 0) and (icount < 100):
            time.sleep(1)
            icount += 1
            p, stdout_data = self.exec_sh(cmd)
        if len(stdout_data) == 0 :
            raise Exception(job_name + " status cannot be found")
        return stdout_data.strip().split()[5]

    def get_job_node(self,
                     job_name,
                     ):
        if job_name.isdigit():
            job_id = job_name
        else:
            job_id = self.get_job_id(job_name)
        cmd = "squeue -j"
        cmd += " " + job_id
        cmd += " | grep " + job_id
        cmd += " | head -1"
        stdout_data = ""
        icount = 0
        while (len(stdout_data) == 0) and (icount < 100):
            time.sleep(1)
            icount += 1
            p, stdout_data = self.exec_sh(cmd)
        if len(stdout_data) == 0 :
            raise Exception(job_name + " status cannot be found")
        return stdout_data.strip().split()[7]

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
        cmd = "scancel -b " + job_id
        self.exec_sh(cmd)
