import sys
import glob
import pkgutil
import os
import fnmatch
from setuptools import setup
from pycmm.settings import DNASEQ_SLURM_MONITOR_PIPELINE_BIN
from pycmm.settings import DUMMY_TABLE_ANNOVAR_BIN
from pycmm.settings import MUTREP_SLURM_MONITOR_PIPELINE_BIN
from pycmm.settings import MUTREP_FAMILY_REPORT_BIN
from pycmm.settings import MUTREP_SUMMARY_REPORT_BIN
from pycmm.settings import PLINK_SLURM_MONITOR_PIPELINE_BIN
from pycmm.settings import PLINK_HAP_ASSOCS_REPORT_BIN
from pycmm.settings import PLINK_MERGE_HAP_ASSOCS_BIN

def opj(*args):
    path = os.path.join(*args)
    return os.path.normpath(path)

def find_data_files(srcdir, *wildcards, **kw):
    # get a list of all files under the srcdir matching wildcards,
    # returned in a format to be used for install_data
    def walk_helper(arg, dirname, files):
        if '.svn' in dirname:
            return
        names = []
        lst, wildcards = arg
        for wc in wildcards:
            wc_name = opj(dirname, wc)
            for f in files:
                filename = opj(dirname, f)

                if fnmatch.fnmatch(filename, wc_name) and not os.path.isdir(filename):
                    names.append(filename)
        if names:
            lst.append( (dirname, names ) )

    file_list = []
    recursive = kw.get('recursive', True)
    if recursive:
        os.path.walk(srcdir, walk_helper, (file_list, wildcards))
    else:
        walk_helper((file_list, wildcards),
                    srcdir,
                    [os.path.basename(f) for f in glob.glob(opj(srcdir, '*'))])
    return file_list

#csv_files = find_data_files('data/', '*.csv')
all_data_files = find_data_files('data/', '*.*')
#all_data_files = find_data_files('script/', '*.*')

setup(
    name='pyCMM',
    version='0.0.1',
    author='Jessada Thutkawkorapin',
    author_email='jessada.thutkawkorapin@gmail.com',
    packages=['pycmm',
              'pycmm.app',
              'pycmm.utils',
              'pycmm.cmmlib',
              'pycmm.flow',
              ],
    scripts=['bin/'+DNASEQ_SLURM_MONITOR_PIPELINE_BIN,
             'bin/pyCMM-dnaseq-pipeline',
             'bin/pyCMM-dnaseq-create-job-setup-file',
             'bin/pyCMM-cmmdb-cal-mut-stat',
             'bin/pyCMM-cmmdb-table-annovar',
             'bin/pyCMM-cmmdb-create-job-setup-file',
             'bin/'+DUMMY_TABLE_ANNOVAR_BIN,
             'bin/'+MUTREP_SLURM_MONITOR_PIPELINE_BIN,
             'bin/pyCMM-mutrep-pipeline',
             'bin/pyCMM-mutrep-mutation-reports',
             'bin/'+MUTREP_FAMILY_REPORT_BIN,
             'bin/'+MUTREP_SUMMARY_REPORT_BIN,
             'bin/pyCMM-mutrep-create-job-setup-file',
             'bin/pyCMM-plink-create-job-setup-file',
             'bin/pyCMM-plink-hap-assocs',
             'bin/'+PLINK_SLURM_MONITOR_PIPELINE_BIN,
             'bin/'+PLINK_HAP_ASSOCS_REPORT_BIN,
             'bin/'+PLINK_MERGE_HAP_ASSOCS_BIN,
             ],
    package=['pyCMM'],
#    package_data={'': ['data/CBV/*.cbv']
#                  },
    data_files=all_data_files,
    url='http://pypi.python.org/pypi/pyCMM/',
    license='LICENSE.txt',
    description='Python packages for my sequencing data analysis at Center of Molecular Medicine, Karolinska Institute, Stockholm, Sweden',
    long_description=open('README.md').read(),
    install_requires=[
        "pysam >= 0.7",
        "pyvcf >= 0.6.0",
        "pyaml >= 15.5.7",
        "openpyxl >= 2.3.3",
        "xlsxwriter >= 0.5.3",
        ],
)
