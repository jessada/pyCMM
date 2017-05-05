import yaml
import re
from os.path import isfile
from collections import OrderedDict
from collections import defaultdict
from pycmm.settings import PHENOTYPE_MISSING
from pycmm.settings import PHENOTYPE_UNAFFECTED
from pycmm.settings import PHENOTYPE_AFFECTED
from pycmm.template import pyCMMBase
from pycmm.cmmlib import CMMParams
from pycmm.utils import mylogger

JOBS_SETUP_SAMPLES_INFOS_KEY = "SAMPLES_INFOS"
JOBS_SETUP_DATASET_ID_KEY = "DATASET_ID"
JOBS_SETUP_FAMILIES_LIST_KEY = "FAMILIES"
JOBS_SETUP_FAMILY_ID_KEY = "FAMILY_ID"
JOBS_SETUP_MEMBERS_LIST_KEY = "MEMBERS"
JOBS_SETUP_SAMPLE_ID_KEY = "SAMPLE_ID"
JOBS_SETUP_SAMPLE_GROUP_KEY = "GROUP"
JOBS_SETUP_SAMPLE_PHENOTYPE_KEY = "PHENOTYPE"

SAMPLE_GROUP_MISSING = 0
UNKNOWN_DATASET = "UNKNOWN_DATASET"

NO_FAMILY = "NO_FAMILY"


class Sample(CMMParams):
    """  To parse and structure family member information """

    def __init__(self, member_info, fam_id=None, *args, **kwargs):
        self.__gts = None
        self.__fam_id = fam_id
        kwargs["entries"] = member_info
        super(Sample, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(Sample, self).get_raw_obj_str(*args, **kwargs)
        raw_repr["sample id"] = self.sample_id
        raw_repr["family id"] = self.fam_id
        return raw_repr

    @property
    def sample_id(self):
        return self._get_job_config(JOBS_SETUP_SAMPLE_ID_KEY, required=True)

    @property
    def sample_id_w_fam_pref(self):
        return re.sub(NO_FAMILY+r"[a-zA-Z_0-9]*-",
                      "",
                      str(self.fam_id)+"-"+self.sample_id)

    @property
    def sample_group(self):
        return self._get_job_config(JOBS_SETUP_SAMPLE_GROUP_KEY,
                                    default_val=SAMPLE_GROUP_MISSING)

    @property
    def phenotype(self):
        return self._get_job_config(JOBS_SETUP_SAMPLE_PHENOTYPE_KEY,
                                    default_val=PHENOTYPE_MISSING)

    @property
    def fam_id(self):
        return self.__fam_id

    @property
    def gts(self):
        return self.__gts

    # explicit setting gts to explain possible usage
    @gts.setter
    def gts(self, value):
        self.__gts = value

#class SamplesGroup(list, pyCMMBase):
#    """  To parse information of a group of sample """
#
#    def __init__(self,
#                 *args,
#                 **kwargs
#                 ):
#        list.__init__(self, *args)
#        pyCMMBase.__init__(self, *args, **kwargs)
#        self.__ids = None
#        self.__ids_w_fam_pref = None
##        self.__affected = None
##        self.__unaffected = None
#
#    def __parse_ids(self):
#        if self is None:
#            return []
#        return map(lambda x: x.sample_id, self)
#
#    def __parse_ids_w_fam_pref(self):
#        if self is None:
#            return []
#        return map(lambda x: re.sub(NO_FAMILY+r"[a-zA-Z_0-9]*-", "", x),
#                   map(lambda x: str(x.fam_id)+"-"+x.sample_id,
#                       self))
#
##    def __parse_affected_samples(self):
##        if self is None:
##            return []
##        return SamplesGroup(filter(lambda x: x.phenotype==PHENOTYPE_AFFECTED,
##                                   self))
##
##    def __parse_unaffected_samples(self):
##        if self is None:
##            return []
##        return SamplesGroup(filter(lambda x: x.phenotype==PHENOTYPE_UNAFFECTED,
##                                   self))
##
#    @property
#    def ids(self):
#        if self.__ids is None:
#            self.__ids = self.__parse_ids()
#        return self.__ids
#
#    @property
#    def ids_w_fam_pref(self):
#        if self.__ids_w_fam_pref is None:
#            self.__ids_w_fam_pref = self.__parse_ids_w_fam_pref()
#        return self.__ids_w_fam_pref
#
##    @property
##    def affected(self):
##        if self.__affected is None:
##            self.__affected = self.__parse_affected_samples()
##        return self.__affected
##
##    @property
##    def unaffected(self):
##        if self.__unaffected is None:
##            self.__unaffected = self.__parse_unaffected_samples()
##        return self.__unaffected
##
class Family(CMMParams):
    """  To parse and structure family information """

    def __init__(self, fam_info, *args, **kwargs):
        self.__members = None
        kwargs["entries"] = fam_info
        super(Family, self).__init__(*args, **kwargs)

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(Family, self).get_raw_obj_str(*args, **kwargs)
        raw_repr["family id"] = self.fam_id
        raw_repr["members"] = self.members
        return raw_repr

    @property
    def fam_id(self):
        return self._get_job_config(JOBS_SETUP_FAMILY_ID_KEY, required=True)

    @property
    def members(self):
        if self.__members is None:
            self.__members = map(lambda x: Sample(member_info=x,
                                                  fam_id=self.fam_id),
                                 self._get_job_config(JOBS_SETUP_MEMBERS_LIST_KEY,
                                                      required=True)
                                 )
        return self.__members

    @property
    def samples_list(self):
        return self.members

class Dataset(CMMParams):
    """  To parse and structure dataset information """

    def __init__(self, dataset_info, *args, **kwargs):
        self.__members = None
        kwargs["entries"] = dataset_info
        super(Dataset, self).__init__(*args, **kwargs)
        self.__parse_families(dataset_info)
        self.__parse_samples_list()

    def get_raw_obj_str(self, *args, **kwargs):
        raw_repr = super(Family, self).get_raw_obj_str(*args, **kwargs)
        raw_repr["daaset id"] = self.dataset_id
        raw_repr["samples list"] = self.samples_list
        return raw_repr

    def __parse_families(self, dataset_info):
        if dataset_info is None:
            self.__families = None
            return
        families_info = self._get_job_config(JOBS_SETUP_FAMILIES_LIST_KEY, required=True)
        self.__families = OrderedDict()
        for entry in families_info:
            family = Family(entry)
            self.__families[family.fam_id] = family

    def __parse_samples_list(self):
        if self.families is None:
            self.__samples_list = None
            return
        self.__samples_list = reduce(lambda x, y: x+y, 
                                     map(lambda x: self.__families[x].members,
                                         self.__families
                                         )
                                     )

    @property
    def _id(self):
        return self.dataset_id

    @property
    def dataset_id(self):
        return self._get_job_config(JOBS_SETUP_DATASET_ID_KEY, required=True)

    @property
    def families(self):
        return self.__families

    @property
    def samples_list(self):
        return self.__samples_list

class SamplesInfo(pyCMMBase):
    """  To parse and structure all samples information """

    def __init__(self,
                 samples_info,
                 item_template=Family,
                 *args,
                 **kwargs
                 ):
        super(SamplesInfo, self).__init__(*args, **kwargs)
        self.__parse_items(samples_info, item_template)
        self.__samples_list = self.__parse_samples_list()
#        self.__samples_groups = self.__parse_samples_groups()

    def __parse_items(self, samples_info, item_template):
        if samples_info is None:
            self.__items = None
            return
        self.__items = OrderedDict()
        for entry in samples_info:
            item = item_template(entry)
            self.__items[item._id] = item

    def __parse_samples_list(self):
        if self.__items is None:
            return None
        samples_list = reduce(lambda x, y: x+y, 
                              map(lambda x: x.samples_list,
                                  self.__items.values()
                                  )
                              )
        return samples_list

#    def __parse_samples_groups(self):
#        if self.samples_list is None:
#            return None
#        raw_samples_groups = defaultdict(list)
#        for sample in self.samples_list:
#            raw_samples_groups[sample.sample_group].append(sample)
#        samples_groups = defaultdict(list)
#        for group_no in raw_samples_groups:
#            samples_groups[group_no] = SamplesGroup(raw_samples_groups[group_no])
#        return samples_groups
#
    @property
    def has_info(self):
        return self.__items is not None

    @property
    def items(self):
        return self.__items

    @property
    def samples_list(self):
        return self.__samples_list

    @property
    def samples_dict(self):
        if self.__samples_list is None:
            return None
        return dict((x.sample_id, x) for x in self.samples_list)

#    @property
#    def samples_groups(self):
#        return self.__samples_groups
#
    @property
    def samples_id(self):
        if self.samples_list is None:
            return None
        return map(lambda x: x.sample_id, self.samples_list)

    @property
    def samples_id_w_fam_pref(self):
        if self.samples_list is None:
            return None
        return map(lambda x: x.sample_id_w_fam_pref, self.samples_list)

    @property
    def families(self):
        return self.__items

#    @property
#    def affected_samples(self):
#        return self.samples_list.affected
#
#    @property
#    def unaffected_samples(self):
#        return self.samples_list.unaffected
#
class MutRepSamplesInfo(SamplesInfo):
    """  To parse and structure all samples information for mutation report """

    def __init__(self,
                 samples_info,
                 *args,
                 **kwargs
                 ):
        super(MutRepSamplesInfo, self).__init__(samples_info,
                                                Dataset,
                                                *args,
                                                **kwargs)
        self.__parse_families()

    def __parse_families(self):
        if self.datasets is None:
            self.__families = None
            return
        self.__families = OrderedDict()
        for dataset in self.datasets.values():
            for fam_id in dataset.families:
                family = dataset.families[fam_id]
                self.__families[fam_id] = family

    @property
    def datasets(self):
        return self.items

    @property
    def families(self):
        return self.__families

def params_to_yaml_doc(*args, **kwargs):
    yaml_doc = {}
    if 'sample_info' not in kwargs:
        return yaml_doc
    sample_info = kwargs['sample_info']
    if sample_info is None:
        return yaml_doc
    if isfile(sample_info):
        s_stream = file(sample_info, "r")
        document = yaml.safe_load(s_stream)
        samples_info = document[JOBS_SETUP_SAMPLES_INFOS_KEY]
        # Currently, the parser below can only support sample.info that are 
        # compatible with MutRepSamplesInfo. So if there the files with other
        # format, the parser below need to be refactored
        for dataset_info in samples_info:
            dataset_id = dataset_info[JOBS_SETUP_DATASET_ID_KEY]
            dataset_info[JOBS_SETUP_DATASET_ID_KEY] = '"' + dataset_id + '"'
            families_info = dataset_info[JOBS_SETUP_FAMILIES_LIST_KEY]
            for family_info in families_info:
                fam_id = family_info[JOBS_SETUP_FAMILY_ID_KEY]
                family_info[JOBS_SETUP_FAMILY_ID_KEY] = '"' + str(fam_id) + '"'
                for member_info in family_info[JOBS_SETUP_MEMBERS_LIST_KEY]:
                    sample_id = member_info[JOBS_SETUP_SAMPLE_ID_KEY]
                    member_info[JOBS_SETUP_SAMPLE_ID_KEY] = '"' + sample_id + '"'
        yaml_doc[JOBS_SETUP_SAMPLES_INFOS_KEY] = samples_info
        return yaml_doc
    # Dataset ID cannnot be identified as Samples info are parsed as text
    # so Dataset ID will be identified as UNKNOWN_DATASET
    families = []
    samples_w_no_fam = []
    for raw_family_info in sample_info.split(","):
        if raw_family_info.find(":") == -1:
            member_info = {}
            member_info[JOBS_SETUP_SAMPLE_ID_KEY] = '"' + raw_family_info + '"'
            samples_w_no_fam.append(member_info)
        else:
            family_info = {}
            family_items = raw_family_info.split(":")
            family_info[JOBS_SETUP_FAMILY_ID_KEY] = '"' + family_items[0] + '"'
            members_info = []
            for sample_id in family_items[1:]:
                member_info = {}
                member_info[JOBS_SETUP_SAMPLE_ID_KEY] = '"' + sample_id + '"'
                members_info.append(member_info)
            family_info[JOBS_SETUP_MEMBERS_LIST_KEY] = members_info
            families.append(family_info)
    if len(samples_w_no_fam) > 0:
        family_info = {}
        family_info[JOBS_SETUP_FAMILY_ID_KEY] = '"' + NO_FAMILY + '"'
        family_info[JOBS_SETUP_MEMBERS_LIST_KEY] = samples_w_no_fam
        families.append(family_info)
    datasets = []
    dataset = {}
    dataset[JOBS_SETUP_DATASET_ID_KEY] = UNKNOWN_DATASET
    dataset[JOBS_SETUP_FAMILIES_LIST_KEY] = families
    datasets.append(dataset)
    yaml_doc[JOBS_SETUP_SAMPLES_INFOS_KEY] = datasets
    return yaml_doc
