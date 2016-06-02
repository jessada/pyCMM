import yaml
from os.path import isfile
from collections import OrderedDict
from pycmm.template import pyCMMBase
from pycmm.cmmlib import CMMParams

JOBS_SETUP_SAMPLES_INFOS_KEY = "SAMPLES_INFOS"
JOBS_SETUP_FAMILY_ID_KEY = "FAMILY_ID"
JOBS_SETUP_MEMBERS_LIST_KEY = "MEMBERS"
JOBS_SETUP_SAMPLE_ID_KEY = "SAMPLE_ID"

NO_FAMILY = "NO_FAMILY"


class Sample(CMMParams):
    """  To parse and structure family member information """

    def __init__(self, member_info, fam_id=None, **kwargs):
        self.__gts = None
        self.__fam_id = fam_id
        kwargs["entries"] = member_info
        super(Sample, self).__init__(**kwargs)

    def get_raw_repr(self, **kwargs):
        raw_repr = super(Sample, self).get_raw_repr(**kwargs)
        raw_repr["sample id"] = self.sample_id
        return raw_repr

    @property
    def sample_id(self):
        return self._get_job_config(JOBS_SETUP_SAMPLE_ID_KEY, required=True)

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

class Family(CMMParams):
    """  To parse and structure family information """

    def __init__(self, fam_info, **kwargs):
        self.__members = None
        kwargs["entries"] = fam_info
        super(Family, self).__init__(**kwargs)

    def get_raw_repr(self, **kwargs):
        raw_repr = super(Family, self).get_raw_repr(**kwargs)
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

class SamplesInfo(pyCMMBase):
    """  To parse and structure all samples information """

    def __init__(self,
                 samples_info,
                 family_template=Family,
                 **kwargs):
        super(SamplesInfo, self).__init__(**kwargs)
        self.__parse_families(samples_info, family_template)

    def __parse_families(self, samples_info, family_template):
        if samples_info is None:
            self.__families = None
            return
        self.__families = OrderedDict()
        for entry in samples_info:
            family = family_template(fam_info=entry)
            self.__families[family.fam_id] = family

    @property
    def families(self):
        return self.__families

    @property
    def samples_dict(self):
        if self.families is None:
            return None
        return dict((x.sample_id, x) for x in self.samples_list)

    @property
    def samples_list(self):
        if self.families is None:
            return None
        return reduce(lambda x, y: x+y, 
                      map(lambda x: self.__families[x].members,
                          self.__families
                          )
                      )

    @property
    def samples_id(self):
        if self.samples_list is None:
            return None
        return map(lambda x: x.sample_id, self.samples_list)

    @property
    def samples_id_w_fam_pref(self):
        if self.samples_list is None:
            return None
        return map(lambda x: x.fam_id+"-"+x.sample_id,
                   self.samples_list)

def params_to_yaml_doc(**kwargs):
    yaml_doc = {}
    if 'sample_info' not in kwargs:
        return yaml_doc
    sample_info = kwargs['sample_info']
    if sample_info is None:
        return yaml_doc
    if isfile(sample_info):
        s_stream = file(sample_info, "r")
        document = yaml.safe_load(s_stream)
        yaml_doc[JOBS_SETUP_SAMPLES_INFOS_KEY] = document[JOBS_SETUP_SAMPLES_INFOS_KEY]
        return yaml_doc
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
    yaml_doc[JOBS_SETUP_SAMPLES_INFOS_KEY] = families
    return yaml_doc
