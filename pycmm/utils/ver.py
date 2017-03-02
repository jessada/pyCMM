import os
import re
import pkg_resources
from pycmm.settings import DUMMY_TABLE_ANNOVAR_BASH
from pycmm.template import pyCMMBase
from pycmm.utils import exec_sh

class VersionManager(pyCMMBase):
    """ to semi-manaully handle version control of all modules used in pyCMM """

    def __init__(self, **kwargs):
        super(VersionManager, self).__init__(**kwargs)
    
    def get_raw_obj_str(self):
        return {"pycmm": self.pycmm_version,
                "pysam": self.pysam_version,
                "pyvcf": self.pyvcf_version,
                "pyaml": self.pyaml_version,
                "openpyxl": self.openpyxl_version,
                "xlsxwriter": self.xlsxwriter_version,
                "plink": self.plink_version,
                "table_annovar": self.table_annovar_version,
                "GATK": self.gatk_version,
                "vcftools": self.vcftools_version,
                }

    @property
    def pycmm_version(self):
        return pkg_resources.get_distribution("pycmm").version

    @property
    def pysam_version(self):
        return pkg_resources.get_distribution("pysam").version

    @property
    def pyvcf_version(self):
        return pkg_resources.get_distribution("pyvcf").version

    @property
    def pyaml_version(self):
        return pkg_resources.get_distribution("pyaml").version

    @property
    def openpyxl_version(self):
        return pkg_resources.get_distribution("openpyxl").version

    @property
    def xlsxwriter_version(self):
        return pkg_resources.get_distribution("xlsxwriter").version

    @property
    def plink_version(self):
        return "1.07"

    @property
    def table_annovar_version(self):
        p, data = exec_sh(DUMMY_TABLE_ANNOVAR_BASH + " -h", silent=True)
        search_result = re.search("\$(.*)\$", data)
        return search_result.group(1)

    @property
    def gatk_version(self):
        gatk_path = os.environ.get("GATK_dir")
        search_result = re.search("(\d+\.){0,3}(\d+)$", gatk_path)
        if search_result is None:
            return None
        return search_result.group(0)

    @property
    def vcftools_version(self):
        p, data = exec_sh("vcftools", silent=True)
        search_result = re.search("(\d+\.){0,3}(\d+)", data)
        return search_result.group(0)

