from magiconfig import MagiConfig
from builtins import dict
from builtins import float
from builtins import int
from builtins import list
from builtins import str
from builtins import tuple

config = MagiConfig()
config.dataset = MagiConfig()
config.dataset.background = {'QCD': ['tree_QCD_Pt_300to470_NN_M0', 'tree_QCD_Pt_470to600_NN', 'tree_QCD_Pt_600to800_NN', 'tree_QCD_Pt_800to1000_NN', 'tree_QCD_Pt_1000to1400_NN', 'tree_QCD_Pt_1400to1800_NN', 'tree_QCD_Pt_1800to2400_NN', 'tree_QCD_Pt_2400to3200_NN'], 'TTJets': ['tree_TTJets_DiLept_genMET-150_NN', 'tree_TTJets_DiLept_NN_M0', 'tree_TTJets_HT-1200to2500_NN', 'tree_TTJets_HT-2500toInf_NN', 'tree_TTJets_HT-600to800_NN', 'tree_TTJets_HT-800to1200_NN', 'tree_TTJets_Incl_NN_M0', 'tree_TTJets_SingleLeptFromTbar_genMET-150_NN_M0', 'tree_TTJets_SingleLeptFromTbar_NN_M0', 'tree_TTJets_SingleLeptFromT_genMET-150_NN_M0', 'tree_TTJets_SingleLeptFromT_NN_M0']}
config.dataset.path = '/work1/cms_svj/keane/trainingFilesEventLevel/trainingFiles/'
config.dataset.sample_fractions = [0.7, 0.15, 0.15]
config.dataset.signal = {'signal': ['tree_SVJ_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1_NN', 'tree_SVJ_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1_NN']}
config.features = MagiConfig()
config.features.eventVariables = ['met', 'mT', 'METrHT_pt30']
config.features.jetVariables = ['jPtAK8', 'jEtaAK8', 'jPhiAK8']
config.hyper = MagiConfig()
config.hyper.batchSize = 512
config.hyper.dropout = 0.3
config.hyper.epochs = 10
config.hyper.lambdaDC = 0.0
config.hyper.lambdaGR = 1.0
config.hyper.lambdaReg = 0.0
config.hyper.lambdaTag = 1.0
config.hyper.learning_rate = 0.001
config.hyper.num_classes = 3
config.hyper.num_of_layers = 1
config.hyper.num_of_nodes = 200
config.hyper.rseed = 100
config.model = 'net.pth'
config.outf = 'logTestEL'