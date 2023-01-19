import particlenet_pf
from magiconfig import MagiConfig
import torch
import copy
from utils.python.svjgnntagger import SVJGNNTagger
from utils.data.GNNTagger.SVJTagger import SVJTagger

# load in local model
gnn_local = SVJGNNTagger(score_tag='score',
                local_path='./particleNetModel.pth',
                model_structure='utils.data.GNNTagger.SVJTagger',
                model_inputs='./utils/data/GNNTagger/svj.yaml',
                dec_thresh=0.999)
local_model = gnn_local.load_model()

# convert model with jit and save to file
jit_model = torch.jit.script(local_model)
JIT_PATH = 'model.pt'
torch.jit.save(jit_model, JIT_PATH)