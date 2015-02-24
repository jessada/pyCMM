import sys
import datetime
from pycmm.utils import mylogger
from pycmm.utils import disp
from pycmm import settings
#import numpy as np
#import combivep.settings as cbv_const
#from combivep.preproc.dataset import DataSetManager
#from combivep.engine.wrapper import Trainer
#from combivep.engine.wrapper import Predictor
#from combivep.refdb.control import UcscController
#from combivep.refdb.control import LjbController


def vcf2avdb(*args, **kwargs):
    func_name = sys._getframe().f_code.co_name
    time_stamp = datetime.datetime.now()
    log_file = kwargs["log_file"]
    mylogger.init(log_file=log_file,
                  dbg_mode=kwargs["debug_mode"],
                  )

    disp.new_section_txt("S T A R T <" + func_name + ">")
    required_params={}
    required_params['vcf input file'] = kwargs['vcf_input_file']
    required_params['avdb output file'] = kwargs['avdb_output_file']
    optional_params={}
    optional_params['log file'] = log_file
    disp.show_config(app_description=settings.VCF2AVDB_DESCRIPTION,
                     time_stamp=time_stamp,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.new_section_txt("F I N I S H <" + func_name + ">")

#def app_combivep_trainer():
#    descriptiion = "An application to train the consensus predictor using "
#    decscription += "Multilayer Perceptron Learning technique."
#    argp = argparse.ArgumentParser(description=description)
#    tmp_help = []
#    tmp_help.append("A file with list of SNPs in CBV (CombiVEP) format.")
#    tmp_help.append("CBV format consists of five fields:")
#    tmp_help.append("CHROM, POS, REF, ALT, ACTUAL_DELETERIOUS_EFFECT.")
#    tmp_help.append("Each field is separated by a tab.")
#    tmp_help.append("SNP Position(POS) is 1-based index.")
#    argp.add_argument('training_data_file', help=' '.join(tmp_help))
#    args = argp.parse_args()
#    train_combivep_using_cbv_data(args.training_data_file)
#
#
#def app_combivep_predictor():
#    description = "An application to predict how likely "
#    description += "the given SNPs will have deleterious effect"
#    argp = argparse.ArgumentParser(description=description)
#    tmp_help = []
#    tmp_help.append("A file with list of SNPs.")
#    tmp_help.append("It can be either in standard VCF format")
#    tmp_help.append("or CBV (CombiVEP) format.")
#    tmp_help.append("CBV format consists of five fields:")
#    tmp_help.append("CHROM, POS, REF, ALT, ACTUAL_DELETERIOUS_EFFECT.")
#    tmp_help.append("Each field is separated by a tab.")
#    tmp_help.append("SNP Position(POS) is 1-based index.")
#    argp.add_argument('input_file', help=' '.join(tmp_help))
#    argp.add_argument('-F',
#                      help='Input format (%(choices)s). Default is in %(default)s format.',
#                      choices=[cbv_const.FILE_TYPE_VCF,
#                               cbv_const.FILE_TYPE_CBV],
#                      metavar='FORMAT',
#                      default=cbv_const.FILE_TYPE_VCF,
#                      dest='input_format'
#                      )
#    args = argp.parse_args()
#    predict_deleterious_probability(args.input_file,
#                                    file_type=args.input_format)
#
#
#def train_combivep_using_cbv_data(training_data_file,
#                                  params_out_file=cbv_const.USER_PARAMS_FILE,
#                                  random_seed=cbv_const.DFLT_SEED,
#                                  n_hidden_nodes=cbv_const.DFLT_HIDDEN_NODES,
#                                  figure_dir=cbv_const.DFLT_FIGURE_DIR,
#                                  iterations=cbv_const.DFLT_ITERATIONS,
#                                  cfg_file=cbv_const.CBV_CFG_FILE,
#                                  ):
#    """
#
#    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
#    CBV has 5 fields:
#    - CHROM
#    - POS
#    - REF
#    - ALT
#    - EFFECT(1=deleterious, 0=neutral)
#    All are tab separated.
#    Required arguments
#    - neutral_data_file : list of SNPs with no harmful effect, CBV format
#    - pathognice_data_file : list of SNPs with deleterious effect, CBV format
#
#    """
#    #pre-processing dataset
#    print >> sys.stderr, 'pre-processing dataset, this may take a while (around 750 SNPs/mins). . .'
#    dm = DataSetManager(cfg_file=cfg_file)
#    dm.load_data(training_data_file, file_type=cbv_const.FILE_TYPE_CBV)
#    dm.validate_data()
#    dm.calculate_scores()
#    dm.set_shuffle_seed(random_seed)
#    dm.shuffle_data()
#    dm.partition_data()
#
#    #partition data
#    training_data   = dm.get_training_data()
#    validation_data = dm.get_validation_data()
#
#    #train !!!
#    print >> sys.stderr, 'Training CombiVEP, please wait (around 500 SNPs/mins) . . .'
#    trainer = Trainer(training_data,
#                      validation_data,
#                      random_seed,
#                      n_hidden_nodes,
#                      figure_dir)
#    trainer.train(iterations)
#    if not os.path.exists(cbv_const.USER_PARAMS_DIR):
#        os.makedirs(cbv_const.USER_PARAMS_DIR)
#    trainer.export_best_parameters(params_out_file)
#
#
#def predict_deleterious_probability(SNPs_file,
#                                    params_file=cbv_const.USER_PARAMS_FILE,
#                                    file_type=cbv_const.FILE_TYPE_VCF,
#                                    output_file=None,
#                                    cfg_file=cbv_const.CBV_CFG_FILE,
#                                    ):
#    """
#
#    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
#    CBV has 5 fields:
#    - CHROM
#    - POS
#    - REF
#    - ALT
#    - EFFECT(1=deleterious, 0=neutral)
#    All are tab separated.
#    Required arguments
#    - SNPs_file : list of SNPs to be predicted, can be either VCF or
#                  CBV (default is VCF)
#
#    """
#    #pre-processing test dataset
#    print >> sys.stderr, 'pre-processing dataset, this may take a while (around 750 SNPs/mins). . .'
#    dm = DataSetManager(cfg_file=cfg_file)
#    dm.load_data(SNPs_file, file_type=file_type)
#    dm.validate_data()
#    dm.calculate_scores()
#
#    #predict
#    predictor = Predictor()
#    predictor.import_parameters(params_file=params_file)
#    out = (np.array(predictor.predict(dm.dataset)).reshape(-1,))
#
#    #print output
#    if output_file is not None:
#        sys.stdout = open(output_file, 'w')
#    tmp_rec = []
#    tmp_rec.append("CHROM")
#    tmp_rec.append("POS")
#    tmp_rec.append("REF")
#    tmp_rec.append("ALT")
#    tmp_rec.append("ACTUAL_DELETERIOUS_EFFECT")
#    tmp_rec.append("PREDICTED_DELETERIOUS_PROBABILITY")
#    tmp_rec.append("PHYLOP_SCORE")
#    tmp_rec.append("SIFT_SCORE")
#    tmp_rec.append("PP2_SCORE")
#    tmp_rec.append("LRT_SCORT")
#    tmp_rec.append("MT_SCORE")
#    tmp_rec.append("GERP_SCORE")
#    print "#" + "\t".join(tmp_rec)
#    for i in xrange(len(dm.dataset)):
#        del tmp_rec[:]
#        snp_data = dm.dataset[i][cbv_const.KW_SNP_DATA]
#        scores = dm.dataset[i][cbv_const.KW_SCORES]
#        tmp_rec.append(snp_data.chrom)
#        tmp_rec.append(snp_data.pos)
#        tmp_rec.append(snp_data.ref)
#        tmp_rec.append(snp_data.alt)
#        tmp_rec.append("%s" % snp_data.target)
#        tmp_rec.append("%6.4f" % out[i])
#        tmp_rec.append(scores.phylop_score)
#        tmp_rec.append(scores.sift_score)
#        tmp_rec.append(scores.pp2_score)
#        tmp_rec.append(scores.lrt_score)
#        tmp_rec.append(scores.mt_score)
#        tmp_rec.append(scores.gerp_score)
#        print "\t".join(tmp_rec)
#    sys.stdout = sys.__stdout__
