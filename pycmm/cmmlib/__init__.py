from collections import OrderedDict
from pycmm.template import pyCMMBase
from pycmm.utils import get_dict_val


class CMMParams(pyCMMBase):
    """  A template to parse jobs parameters in dictionary format  """

    def __init__(self, entries, *args, **kwargs):
        self.__entries = entries
        super(CMMParams, self).__init__(*args, **kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        return raw_repr

    def _get_job_config(self, key, required=False, default_val=None):
        return get_dict_val(self.__entries,
                            key,
                            required=required,
                            default_val=default_val,
                            )
