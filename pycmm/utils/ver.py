import pkg_resources
from pycmm.template import pyCMMBase
from pycmm.utils import exec_sh

class VersionManager(pyCMMBase):
    """ to semi-manaully handle version control of all modules used in pyCMM """

    def __init__(self, **kwargs):
        super(VersionManager, self).__init__(**kwargs)
    
    def get_raw_repr(self):
        return "None"
#        return {"job name": self.job_name,
#                "project code": self.project_code,
#                }

    @property
    def pyaml_version(self):
        return pkg_resources.get_distribution("pyaml").version
