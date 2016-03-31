from collections import OrderedDict
from pycmm.template import pyCMMBase

JOBS_SETUP_SAMPLE_INFOS_KEY = "SAMPLE_INFOS"
JOBS_SETUP_FAMILY_ID_KEY = "FAMILY_ID"
JOBS_SETUP_MEMBERS_LIST_KEY = "MEMBERS"
JOBS_SETUP_SAMPLE_ID_KEY = "SAMPLE_ID"


class MemberInfo(pyCMMBase):
    """  To structure family member information """

    def __init__(self, info, **kwargs):
        self.__info = info
        super(MemberInfo, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["sample id"] = self.sample_id
        return raw_repr

    @property
    def sample_id(self):
        return self.__info[JOBS_SETUP_SAMPLE_ID_KEY]

class FamilyInfo(pyCMMBase):
    """  To structure family information """

    def __init__(self, info, **kwargs):
        self.__info = info
        super(FamilyInfo, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["family"] = self.fam_id
        raw_repr["members"] = self.members
        return raw_repr

    @property
    def fam_id(self):
        return str(self.__info[JOBS_SETUP_FAMILY_ID_KEY])

    @property
    def members(self):
        return map(lambda x: MemberInfo(info=x),
                   self.__info[JOBS_SETUP_MEMBERS_LIST_KEY])

def params_to_yaml(**kwargs):
    yaml = {}
    if 'sample_info' not in kwargs:
        return yaml
    sample_info = kwargs['sample_info']
    if sample_info is None:
        return yaml
    if sample_info.find(":") == -1:
        yaml[JOBS_SETUP_SAMPLE_INFOS_KEY] = sample_info
    elif sample_info.find(":") > -1:
        families_info = []
        for raw_family_info in sample_info.split(","):
            family_info = {}
            family_items = raw_family_info.split(":")
            family_info[JOBS_SETUP_FAMILY_ID_KEY] = '"' + family_items[0] + '"'
            members_info = []
            for family_item in family_items[1:]:
                member_info = {}
                member_info[JOBS_SETUP_SAMPLE_ID_KEY] = '"' + family_item + '"'
                members_info.append(member_info)
            family_info[JOBS_SETUP_MEMBERS_LIST_KEY] = members_info
            families_info.append(family_info)
        yaml[JOBS_SETUP_SAMPLE_INFOS_KEY] = families_info
    return yaml

def extract_families_info(yaml):
    if JOBS_SETUP_SAMPLE_INFOS_KEY not in yaml:
        return None
    sample_infos = yaml[JOBS_SETUP_SAMPLE_INFOS_KEY]
    if type(sample_infos) is str:
        return None
    else:
        families_info = {}
        for entry in sample_infos:
            family_info = FamilyInfo(entry)
            families_info[family_info.fam_id] = family_info
        return families_info

def extract_samples_list(yaml):
    if JOBS_SETUP_SAMPLE_INFOS_KEY not in yaml:
        return None
    sample_infos = yaml[JOBS_SETUP_SAMPLE_INFOS_KEY]
    if type(sample_infos) is str:
        return sample_infos.split(",")
    else:
        samples_list = []
        families_info = extract_families_info(yaml)
        for fam_id in families_info:
            for member in families_info[fam_id].members:
                samples_list.append(member.sample_id)
        return samples_list
