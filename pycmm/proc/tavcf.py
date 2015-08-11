from vcf import Reader as VCFReader
from vcf import RESERVED_INFO
from collections import OrderedDict
from pycmm.utils import mylogger


class TableAnnovarVCFReader(VCFReader):
    """
    An encapsulated version of vcf.Reader from pyVCF package to fix
      - _parse_info  <see _parse_info function>
    """

    def __init__(self,
                 fsock=None,
                 filename=None,
                 compressed=None,
                 prepend_chr=False,
                 strict_whitespace=False,
                 ):
        super(TableAnnovarVCFReader, self).__init__(fsock=fsock,
                                                    filename=filename,
                                                    compressed=compressed,
                                                    prepend_chr=prepend_chr,
                                                    strict_whitespace=strict_whitespace,
                                                    )
        self.__parse_annovar_infos()

    @property
    def annovar_infos(self):
        return self.__annovar_infos

    def __parse_annovar_infos(self):
        """
        To get the list of infos field annotated by ANNOVAR
        by detecting any INFO in metaheader that end with "ANNOVAR"
        """
        self.__annovar_infos = OrderedDict()
        for info_key in self.infos:
            val = self.infos[info_key]
            if (val.desc is not None) and (val.desc.endswith("ANNOVAR")):
                self.__annovar_infos[info_key] = val
        return self.__annovar_infos

    def _parse_info(self, info_str):
        """
        Remake version of pyVCF.parser._parse_info

        Parse the INFO field of a VCF entry into a dictionary of Python
        types.

        ***** Remake part *****
        - convert list that unintentionally indicated by comma back to string
        - decode hex value coded for reserved character (";", "=", etc.)
        - detect if any 'info' entries are annotated by ANNOVAR more than once
          , if yes convert them into 'list'
        """
        if info_str == '.':
            return {}

        entries = info_str.split(';')
        retdict = OrderedDict()

        for entry in entries:
            entry = entry.split('=', 1)
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if (ID in self.annovar_infos) and (entry_type in ('String', 'Character')):
                val = entry[1]
                # decode hex character
                hex_idx = val.find("\\x")
                while hex_idx > 0:
                    hex_str = val[hex_idx:hex_idx+4]
                    decode_chr = val[hex_idx+2:hex_idx+4].decode('hex')
                    val = val.replace(hex_str, decode_chr)
                    hex_idx = val.find("\\x")
            elif entry_type == 'Integer':
                vals = entry[1].split(',')
                try:
                    val = self._map(int, vals)
                # Allow specified integers to be flexibly parsed as floats.
                # Handles cases with incorrectly specified header types.
                except ValueError:
                    val = self._map(float, vals)
            elif entry_type == 'Float':
                vals = entry[1].split(',')
                val = self._map(float, vals)
            elif entry_type == 'Flag':
                val = True
            elif entry_type in ('String', 'Character'):
                try:
                    vals = entry[1].split(',') # commas are reserved characters indicating multiple values
                    val = self._map(str, vals)
                except IndexError:
                    entry_type = 'Flag'
                    val = True

            try:
                if self.infos[ID].num == 1 and entry_type not in ( 'Flag', ):
                    val = val[0]
            except KeyError:
                pass

            # Check if there are infos of more than one alternate alleles annotated by ANNOVAR
            if (ID in self.annovar_infos) and (ID in retdict) and (type(retdict[ID]) is list):
                retdict[ID].append(val)
                val = retdict[ID]
            elif (ID in self.annovar_infos) and (ID in retdict):
                val = [retdict[ID], val]
            retdict[ID] = val

        return retdict

