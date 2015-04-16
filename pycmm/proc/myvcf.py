import sys
import vcf
from random import random
from vcf.parser import _Info as VcfInfo
from vcf.parser import field_counts as vcf_field_counts
from pycmm.utils import mylogger
from collections import OrderedDict

def cal_zygo(vcf_in_file,
             zg_out_file,
             ):
    pass

def write_vcf(vcf_in_file,
              vcf_out_file,
              ):
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    vcf_reader = vcf.Reader(filename=vcf_in_file)
    vcf_reader.infos['EXAMPLE_INFO'] = VcfInfo(
                'EXAMPLE_INFO', vcf_field_counts['A'], 'Float',
                    'Example info value for each allele')
    vcf_writer = vcf.Writer(open(vcf_out_file, 'w'), vcf_reader)
    f_keys = vcf_reader.formats.keys()
#    for info in vcf_reader.infos:
#        mylogger.debug(info)
#        mylogger.debug(vcf_reader.infos[info])
    #for record in vcf_reader:
    record = vcf_reader.next()
    record.add_info('EXAMPLE_INFO', [random() for _ in record.ALT])
    #record.genotype('1052/05') = CallData()
    vcf_writer.write_record(record)

def read_vcf(vcf_in_file,
             ):
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    if vcf_in_file.endswith('.vcf.gz'):
        mylogger.debug('endswith .vcf.gz')
        vcf_reader = vcf.Reader(filename=vcf_in_file)
    else:
        mylogger.debug('does not endswith .vcf.gz')
        vcf_reader = vcf.Reader(open(vcf_in_file, 'r'))
    #for record in vcf_reader.fetch('18', 12512300, 12512375):
    for record in vcf_reader.fetch('18', 12512360, 14513530):
        mylogger.debug(record)
        rec_prop = OrderedDict()
        rec_prop['QUAL'] = record.QUAL
        rec_prop['INFO'] = record.INFO
        rec_prop['aaf'] = record.aaf
#        rec_prop['affected_end'] = record.affected_end
#        rec_prop['affected_start'] = record.affected_start
        rec_prop['alleles'] = record.alleles
        rec_prop['call_rate'] = record.call_rate
        rec_prop['start'] = record.start
        rec_prop['end'] = record.end
        rec_prop['num_called'] = record.num_called
        rec_prop['num_het'] = record.num_het
        rec_prop['list(get_hets())'] = list(record.get_hets())[0:2]
        rec_prop['num_hom_alt'] = record.num_hom_alt
        rec_prop['list(get_hom_alts())[0:2]'] = list(record.get_hom_alts())[0:2]
        rec_prop['num_hom_ref'] = record.num_hom_ref
        rec_prop['list(get_hom_refs())[0:2]'] = list(record.get_hom_refs())[0:2]
        rec_prop['num_unknown'] = record.num_unknown
        rec_prop['list(get_unknowns())[0:2]'] = list(record.get_unknowns())[0:2]
        rec_prop['heterozygosity'] = record.heterozygosity
        rec_prop['is_deletion'] = record.is_deletion
        rec_prop['is_indel'] = record.is_indel
        rec_prop['is_monomorphic'] = record.is_monomorphic
        rec_prop['is_snp'] = record.is_snp
        rec_prop['is_sv'] = record.is_sv
        rec_prop['is_sv_precise'] = record.is_sv_precise
        rec_prop['is_transition'] = record.is_transition
        rec_prop['nucl_diversity'] = record.nucl_diversity
        rec_prop['sv_end'] = record.sv_end
        rec_prop['var_subtype'] = record.var_subtype
        rec_prop['var_subtype'] = record.var_subtype
        rec_prop['vcf_reader.samples[64:67]'] = vcf_reader.samples[64:67]
        rec_prop['record.genotype(vcf_reader.samples[0])'] = record.genotype(vcf_reader.samples[0])
        rec_prop['..reader.samples[0]).gt_bases'] = record.genotype(vcf_reader.samples[0]).gt_bases
        rec_prop['..reader.samples[0]).gt_type'] = record.genotype(vcf_reader.samples[0]).gt_type
        rec_prop['..reader.samples[0]).is_het'] = record.genotype(vcf_reader.samples[0]).is_het
        rec_prop['..reader.samples[0]).is_variant'] = record.genotype(vcf_reader.samples[0]).is_variant
        rec_prop['..reader.samples[0]).phased'] = record.genotype(vcf_reader.samples[0]).phased
        rec_prop['..reader.samples[0]).sample'] = record.genotype(vcf_reader.samples[0]).sample
        rec_prop['..reader.samples[0]).site'] = record.genotype(vcf_reader.samples[0]).site
        rec_prop['record.genotype(vcf_reader.samples[64])'] = record.genotype(vcf_reader.samples[64])
        rec_prop['..reader.samples[64]).gt_bases'] = record.genotype(vcf_reader.samples[64]).gt_bases
        rec_prop['..reader.samples[64]).gt_type'] = record.genotype(vcf_reader.samples[64]).gt_type
        rec_prop['..reader.samples[64]).is_het'] = record.genotype(vcf_reader.samples[64]).is_het
        rec_prop['..reader.samples[64]).is_variant'] = record.genotype(vcf_reader.samples[64]).is_variant
        rec_prop['..reader.samples[64]).phased'] = record.genotype(vcf_reader.samples[64]).phased
        rec_prop['..reader.samples[64]).sample'] = record.genotype(vcf_reader.samples[64]).sample
        rec_prop['..reader.samples[64]).site'] = record.genotype(vcf_reader.samples[64]).site
        rec_prop['record.genotype(vcf_reader.samples[66])'] = record.genotype(vcf_reader.samples[66])
        rec_prop['..reader.samples[66]).gt_bases'] = record.genotype(vcf_reader.samples[66]).gt_bases
        rec_prop['..reader.samples[66]).gt_type'] = record.genotype(vcf_reader.samples[66]).gt_type
        rec_prop['..reader.samples[66]).is_het'] = record.genotype(vcf_reader.samples[66]).is_het
        rec_prop['..reader.samples[66]).is_variant'] = record.genotype(vcf_reader.samples[66]).is_variant
        rec_prop['..reader.samples[66]).phased'] = record.genotype(vcf_reader.samples[66]).phased
        rec_prop['..reader.samples[66]).sample'] = record.genotype(vcf_reader.samples[66]).sample
        rec_prop['..reader.samples[66]).site'] = record.genotype(vcf_reader.samples[66]).site
#        mylogger.debug(rec_prop)
        fmt_rec_prop = " >> {name:<45}: {value}"
        prop_data = []
        for prop_key in rec_prop.keys():
            prop_data.append(fmt_rec_prop.format(name=prop_key,
                                                 value=str(rec_prop[prop_key])))
        mylogger.debug("\n".join(prop_data))
        mylogger.debug(record.genotype('1052/05'))

def load_vcf(vcf_in_file,
             ):
    mylogger.getLogger(__name__ + "." + sys._getframe().f_code.co_name)
    if vcf_in_file.endswith('.vcf.gz'):
        mylogger.debug('endswith .vcf.gz')
        vcf_reader = vcf.Reader(filename=vcf_in_file)
    else:
        mylogger.debug('does not endswith .vcf.gz')
        vcf_reader = vcf.Reader(open(vcf_in_file, 'r'))
#    mylogger.debug(vcf_reader.metadata)
#    for metadata_key in vcf_reader.metadata:
#        mylogger.debug(metadata_key)
#        mylogger.debug(vcf_reader.metadata[metadata_key])
    for record in vcf_reader:
        mylogger.debug(record)
        mylogger.debug(record.QUAL)
        mylogger.debug(record.INFO)
        mylogger.debug(record.genotype('1052/05'))
#    mylogger.debug('')
#    mylogger.debug('')
#    mylogger.debug('                 I   N   F   O                 ')
#    mylogger.debug('')
#    mylogger.debug('')
#    for info in vcf_reader.infos:
#        mylogger.debug(info)
#        mylogger.debug(vcf_reader.infos[info])
#    mylogger.debug('')
#    mylogger.debug('')
#    mylogger.debug('                 F   O   R   M   A   T                 ')
#    mylogger.debug('')
#    mylogger.debug('')
#    for fmt in vcf_reader.formats:
#        mylogger.debug(fmt)
#        mylogger.debug(vcf_reader.formats[fmt])
