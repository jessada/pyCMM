from collections import OrderedDict
from pycmm.template import pyCMMBase
from pycmm.utils import is_number
from pycmm.utils import DefaultOrderedDict

RAW_HAP_ASSOC_LOCUS_IDX = 0
RAW_HAP_ASSOC_HAPLOTYPE_IDX = 1
RAW_HAP_ASSOC_F_A_IDX = 2
RAW_HAP_ASSOC_F_U_IDX = 3
RAW_HAP_ASSOC_CHISQ_IDX = 4
RAW_HAP_ASSOC_DF_IDX = 5
RAW_HAP_ASSOC_P_IDX = 6
RAW_HAP_ASSOC_SNPS_IDX = 7

LMISS_CHR_IDX = 0
LMISS_SNP_IDX = 1
LMISS_CLST_IDX = 2
LMISS_NMISS_IDX = 3
LMISS_N_GENO_IDX = 4
LMISS_N_CLST_IDX = 5
LMISS_F_MISS_IDX = 6

MAP_CHR_IDX = 0
MAP_SNP_IDX = 1
MAP_DIST_IDX = 2
MAP_POS_IDX = 3

SNP_INFO_SNP_IDX = 0
SNP_INFO_F_MISS_A_IDX = 1
SNP_INFO_F_MISS_U_IDX = 2
SNP_INFO_CHROM_IDX = 3
SNP_INFO_POS_IDX = 4

FAM_FAM_ID_IDX = 0
FAM_INDV_ID_IDX = 1
FAM_PAT_ID_IDX = 2
FAM_MAT_ID_IDX = 3
FAM_SEX_IDX = 4
FAM_PHENOTYPE_IDX = 5

UNAFFECTED_PHENOTYPE = "1"
AFFECTED_PHENOTYPE = "2"


class Reader(pyCMMBase):
    """ a base class for read and parse file content """

    def __init__(self, 
                 file_name,
                 parser,
                 **kwargs
                 ):
        self.__file_name = file_name
        self.__reader = open(file_name, 'rt')
        self.__parser = parser
        super(Reader, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["file name"] = self.__file_name
        raw_repr["parser"] = self.__parser
        return raw_repr

    @property
    def reader(self):
        return self.__reader

    def __iter__(self):
        return self

    def next(self):
        return self.__parser(self.__reader.next())

class HapAssocUtils(pyCMMBase):
    """ an util to gather basic info from haplotype association study result """

    def __init__(self, hap_assoc_file, **kwargs):
        self.__hap_assoc_file = hap_assoc_file
        super(HapAssocUtils, self).__init__(**kwargs)

    @property
    def hap_assoc_file(self):
        return self.__hap_assoc_file

    @property
    def nlines(self):
        with open(self.hap_assoc_file, 'r') as f:
            nlines = 0
            for line in f:
                nlines += 1
        return nlines-1

    @property
    def ncols(self):
        with open(self.hap_assoc_file, 'r') as f:
            return len(f.next().strip().split())

class HapAssoc(pyCMMBase):
    """ To parse haplotype association study result """

    def __init__(self, rec, locus_prefix=None, **kwargs):
        self.__items = rec.strip().split()
        self.__ors = None
        if locus_prefix is None:
            self.__locus_prefix = ""
        else:
            self.__locus_prefix = locus_prefix + "_"
        self.__set_idxs()
        super(HapAssoc, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["locus"] = self.locus
        raw_repr["locus prefix"] = self.locus_prefix
        raw_repr["haplotype"] = self.haplotype
        raw_repr["frequency affected"] = self.f_a
        raw_repr["frequency unaffected"] = self.f_u
        raw_repr["chi square"] = self.chisq
        raw_repr["degree of freedom"] = self.df
        raw_repr["p-value"] = self.pvalue
        raw_repr["odd ratios"] = self.ors
        raw_repr["snps"] = self.snps
        return raw_repr

    def __set_idxs(self):
        self.__locus_idx = RAW_HAP_ASSOC_LOCUS_IDX
        self.__haplotype_idx = RAW_HAP_ASSOC_HAPLOTYPE_IDX
        self.__f_a_idx = RAW_HAP_ASSOC_F_A_IDX
        self.__f_u_idx = RAW_HAP_ASSOC_F_U_IDX
        self.__chisq_idx = RAW_HAP_ASSOC_CHISQ_IDX
        self.__df_idx = RAW_HAP_ASSOC_DF_IDX
        self.__p_idx = RAW_HAP_ASSOC_P_IDX
        if len(self.__items) == 8:
            self.__ors_idx = None
            self.__snps_idx = RAW_HAP_ASSOC_SNPS_IDX
        else:
            self.__ors_idx = RAW_HAP_ASSOC_SNPS_IDX
            self.__snps_idx = RAW_HAP_ASSOC_SNPS_IDX + 1

    @property
    def locus_idx(self):
        return self.__locus_idx

    @property
    def haplotype_idx(self):
        return self.__haplotype_idx

    @property
    def f_a_idx(self):
        return self.__f_a_idx

    @property
    def f_u_idx(self):
        return self.__f_u_idx

    @property
    def chisq_idx(self):
        return self.__chisq_idx

    @property
    def df_idx(self):
        return self.__df_idx

    @property
    def p_idx(self):
        return self.__p_idx

    @property
    def ors_idx(self):
        return self.__ors_idx

    @property
    def snps_idx(self):
        return self.__snps_idx

    @property
    def locus_prefix(self):
        return self.__locus_prefix

    @property
    def locus(self):
        locus = self.__items[self.locus_idx]
        if locus == "LOCUS":
            return locus
        return self.locus_prefix + locus

    @property
    def haplotype(self):
        return self.__items[self.haplotype_idx]

    @property
    def f_a(self):
        f_a = self.__items[self.f_a_idx]
        if is_number(f_a):
            return float(f_a)
        else:
            return f_a

    @property
    def f_u(self):
        f_u = self.__items[self.f_u_idx]
        if is_number(f_u):
            return float(f_u)
        else:
            return f_u

    @property
    def chisq(self):
        chisq = self.__items[self.chisq_idx]
        if is_number(chisq):
            return float(chisq)
        else:
            return chisq

    @property
    def df(self):
        return self.__items[self.df_idx]

    @property
    def pvalue(self):
        pvalue = self.__items[self.p_idx]
        if is_number(pvalue):
            return float(pvalue)
        else:
            return pvalue

    @property
    def ors(self):
        if self.__ors is not None:
            return self.__ors
        if self.ors_idx is not None:
            self.__ors = self.__items[self.ors_idx]
            if is_number(self.__ors):
                self.__ors = float(self.__ors)
            return self.__ors
        if not is_number(self.f_a):
            self.__ors = "NA"
            return self.__ors
        f_a = float(self.f_a)
        f_u = float(self.f_u)
        if ((f_a == 1) or
            (f_u == 1) or
            (f_a == 0) or
            (f_u == 0)
            ):
            self.__ors = "NA"
            return self.__ors
        self.__ors = (f_a/(1-f_a))/(f_u/(1-f_u))
        return self.__ors

    @property
    def snps(self):
        return self.__items[self.snps_idx]

class HapAssocReader(Reader):
    """ To read and parse haplotype association study result """

    def __init__(self, 
                 locus_prefix=None,
                 **kwargs
                 ):
        self.__locus_prefix = locus_prefix
        kwargs['parser'] = HapAssoc
        super(HapAssocReader, self).__init__(**kwargs)

    def next(self):
        while True:
            hap_assoc_rec = HapAssoc(self.reader.next(),
                                     locus_prefix=self.__locus_prefix)
            if hap_assoc_rec.haplotype == "OMNIBUS":
                continue
            break
        return hap_assoc_rec

def merge_hap_assocs(hap_assoc_files,
                     out_file,
                     locus_prefixs=None,
                     filter_criteria=None,
                     ):
    writer = open(out_file, 'w')
    header = "LOCUS"
    header += "\tHAPLOTYPE"
    header += "\tF_A"
    header += "\tF_U"
    header += "\tCHISQ"
    header += "\tDF"
    header += "\tP"
    header += "\tORS"
    header += "\tSNPS"
    header += "\n"
    writer.write(header)
    for file_idx in xrange(len(hap_assoc_files)):
        file_name = hap_assoc_files[file_idx]
        if locus_prefixs is None:
            locus_prefix = "File" + str(file_idx+1)
        else:
            locus_prefix = locus_prefixs[file_idx]
        hap_assoc_reader = HapAssocReader(file_name=file_name, locus_prefix=locus_prefix)
        hap_assoc_reader.next()
        for hap_assoc_rec in hap_assoc_reader:
            if filter_criteria is not None:
                if (filter_criteria.filter_pvalue005 and
                    (hap_assoc_rec.pvalue > 0.05)
                    ):
                    continue
                if (filter_criteria.filter_disease_snp and
                    (hap_assoc_rec.f_a < hap_assoc_rec.f_u)
                    ):
                    continue
            content = hap_assoc_rec.locus
            content += "\t" + hap_assoc_rec.haplotype
            content += "\t" + str(hap_assoc_rec.f_a)
            content += "\t" + str(hap_assoc_rec.f_u)
            content += "\t" + str(hap_assoc_rec.chisq)
            content += "\t" + hap_assoc_rec.df
            content += "\t" + str(hap_assoc_rec.pvalue)
            if type(hap_assoc_rec.ors) is str:
                content += "\t" + hap_assoc_rec.ors
            else:
                content += "\t" + "{:4.4f}".format(hap_assoc_rec.ors)
            content += "\t" + hap_assoc_rec.snps
            content += "\n"
            writer.write(content)
    writer.close()

class LMiss(pyCMMBase):
    """ To parse snps missing genotyping information """

    def __init__(self, rec, **kwargs):
        self.__items = rec.strip().split()
        if len(self.__items) != 7:
            raise IOError("Invalid lmiss file")
        super(LMiss, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["chrom"] = self.chrom
        raw_repr["snp"] = self.snp
        raw_repr["clst"] = self.clst
        raw_repr["nmiss"] = self.nmiss
        raw_repr["n_geno"] = self.n_geno
        raw_repr["n_clst"] = self.n_clst
        raw_repr["f_miss"] = self.f_miss
        return raw_repr

    @property
    def chrom(self):
        return self.__items[LMISS_CHR_IDX]

    @property
    def snp(self):
        return self.__items[LMISS_SNP_IDX]

    @property
    def clst(self):
        return self.__items[LMISS_CLST_IDX]

    @property
    def nmiss(self):
        return self.__items[LMISS_NMISS_IDX]

    @property
    def n_geno(self):
        return self.__items[LMISS_N_GENO_IDX]

    @property
    def n_clst(self):
        return self.__items[LMISS_N_CLST_IDX]

    @property
    def f_miss(self):
        f_miss = self.__items[LMISS_F_MISS_IDX]
        if is_number(f_miss):
            return float(f_miss)
        else:
            return f_miss

class LMissReader(Reader):
    """ To read and parse snps missing genotyping information """

    def __init__(self, 
                 **kwargs
                 ):
        kwargs['parser'] = LMiss
        super(LMissReader, self).__init__(**kwargs)

class Map(pyCMMBase):
    """ To parse a record in PLINK map file """

    def __init__(self, rec, **kwargs):
        self.__items = rec.strip().split()
        if len(self.__items) < 4:
            raise IOError("Invalid map file")
        super(Map, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["chrom"] = self.chrom
        raw_repr["snp"] = self.snp
        raw_repr["dist"] = self.dist
        raw_repr["pos"] = self.pos
        return raw_repr

    @property
    def chrom(self):
        return self.__items[MAP_CHR_IDX]

    @property
    def snp(self):
        return self.__items[MAP_SNP_IDX]

    @property
    def dist(self):
        dist = self.__items[MAP_DIST_IDX]
        if is_number(dist):
            return float(dist)
        else:
            return dist

    @property
    def pos(self):
        return self.__items[MAP_POS_IDX]

class MapReader(Reader):
    """ To read and parse file plink map file """

    def __init__(self, 
                 **kwargs
                 ):
        kwargs['parser'] = Map
        super(MapReader, self).__init__(**kwargs)

class SnpInfo(pyCMMBase):
    """ To parse snp information """

    def __init__(self, rec=None, **kwargs):
        if rec is not None:
            items = rec.strip().split()
            if len(items) != 5:
                raise IOError("Invalid snp.info file")
            self.__snp = items[SNP_INFO_SNP_IDX]
            f_miss_a = items[SNP_INFO_F_MISS_A_IDX]
            if is_number(f_miss_a):
                self.__f_miss_a = float(f_miss_a)
            else:
                self.__f_miss_a = f_miss_a
            f_miss_u = items[SNP_INFO_F_MISS_U_IDX]
            if is_number(f_miss_u):
                self.__f_miss_u = float(f_miss_u)
            else:
                self.__f_miss_u = f_miss_u
            self.__chrom = items[SNP_INFO_CHROM_IDX]
            self.__pos = items[SNP_INFO_POS_IDX]
        else:
            self.__snp = None
            self.__f_miss_a = None
            self.__f_miss_u = None
            self.__chrom = None
            self.__pos = None
        self.__snp_idx = None
        super(SnpInfo, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["snp"] = self.snp
        raw_repr["f_miss_a"] = self.f_miss_a
        raw_repr["f_miss_u"] = self.f_miss_u
        raw_repr["chrom"] = self.chrom
        raw_repr["pos"] = self.pos
        raw_repr["snp_idx"] = self.snp_idx
        return raw_repr

    @property
    def snp(self):
        return self.__snp

    @snp.setter
    def snp(self, value):
        self.__snp = value

    @property
    def f_miss_a(self):
        return self.__f_miss_a

    @f_miss_a.setter
    def f_miss_a(self, value):
        self.__f_miss_a = value

    @property
    def f_miss_u(self):
        return self.__f_miss_u

    @f_miss_u.setter
    def f_miss_u(self, value):
        self.__f_miss_u = value

    @property
    def chrom(self):
        return self.__chrom

    @chrom.setter
    def chrom(self, value):
        self.__chrom = value

    @property
    def pos(self):
        return self.__pos

    @pos.setter
    def pos(self, value):
        self.__pos = value

    @property
    def snp_idx(self):
        return self.__snp_idx

    @snp_idx.setter
    def snp_idx(self, value):
        self.__snp_idx = value

class SnpInfoReader(Reader):
    """ To read and parse snp information """

    def __init__(self, 
                 **kwargs
                 ):
        kwargs['parser'] = SnpInfo
        super(SnpInfoReader, self).__init__(**kwargs)

def merge_lmiss_map(lmiss_file, map_file, out_file):
    # read snp info
    snp_info = DefaultOrderedDict(SnpInfo)
    lmiss_reader = LMissReader(file_name=lmiss_file)
    for lmiss_rec in lmiss_reader:
        snp = lmiss_rec.snp
        if lmiss_rec.clst == UNAFFECTED_PHENOTYPE:
            snp_info[snp].snp = snp
            snp_info[snp].f_miss_u = lmiss_rec.f_miss
        elif lmiss_rec.clst == AFFECTED_PHENOTYPE:
            snp_info[snp].snp = snp
            snp_info[snp].f_miss_a = lmiss_rec.f_miss
    map_reader = MapReader(file_name=map_file)
    for map_rec in map_reader:
        snp = map_rec.snp
        snp_info[snp].snp = snp
        snp_info[snp].chrom = map_rec.chrom
        snp_info[snp].pos = map_rec.pos
    # write snp to out file
    writer = open(out_file, 'w')
    header = "SNP"
    header += "\tF_MISS_A"
    header += "\tF_MISS_U"
    header += "\tCHROM"
    header += "\tPOS"
    header += "\n"
    writer.write(header)
    for snp in snp_info:
        content = snp_info[snp].snp
        content += "\t" + str(snp_info[snp].f_miss_a)
        content += "\t" + str(snp_info[snp].f_miss_u)
        content += "\t" + snp_info[snp].chrom
        content += "\t" + snp_info[snp].pos
        content += "\n"
        writer.write(content)
    writer.close()

class Fam(pyCMMBase):
    """ To parse a record in PLINK fam and tfam file """

    def __init__(self, rec, **kwargs):
        self.__items = rec.strip().split()
        if len(self.__items) != 6:
            raise IOError("Invalid fam file")
        super(Fam, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = OrderedDict()
        raw_repr["Family ID"] = self.fam_id
        raw_repr["Individual ID"] = self.indv_id
        raw_repr["Paternal ID"] = self.pat_id
        raw_repr["Maternal ID"] = self.mat_id
        raw_repr["Sex"] = self.sex
        raw_repr["Phenotype"] = self.phenotype
        return raw_repr

    @property
    def fam_id(self):
        return self.__items[FAM_FAM_ID_IDX]

    @property
    def indv_id(self):
        return self.__items[FAM_INDV_ID_IDX]

    @property
    def pat_id(self):
        return self.__items[FAM_PAT_ID_IDX]

    @property
    def mat_id(self):
        return self.__items[FAM_MAT_ID_IDX]

    @property
    def sex(self):
        return self.__items[FAM_SEX_IDX]

    @property
    def phenotype(self):
        return self.__items[FAM_PHENOTYPE_IDX]

class FamReader(Reader):
    """ To read and parse Fam information """

    def __init__(self, 
                 **kwargs
                 ):
        kwargs['parser'] = Fam
        super(FamReader, self).__init__(**kwargs)

class TPed(Map):
    """ To parse a record in PLINK tped file """

    def __init__(self, rec, **kwargs):
        self.__items = rec.strip().split("\t")
        if (len(self.__items) < 5) and (len(self.__items) != 0):
            self.dbg(len(self.__items))
            raise IOError("Invalid TPed file at line: " + rec)
        kwargs['rec'] = rec
        super(TPed, self).__init__(**kwargs)

    def get_raw_repr(self):
        raw_repr = super(TPed, self).get_raw_repr()
        raw_repr["genotypes"] = self.gts
        return raw_repr

    @property
    def gts(self):
        return self.__items[4:]

class TPedReader(Reader):
    """ To read and parse Fam information """

    def __init__(self, 
                 **kwargs
                 ):
        kwargs['parser'] = TPed
        super(TPedReader, self).__init__(**kwargs)
