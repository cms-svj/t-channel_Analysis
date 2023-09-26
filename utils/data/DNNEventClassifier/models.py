import torch
import torch.nn as nn
import torch.nn.functional as f
from torch.autograd import Variable
from utils.data.DNNEventClassifier.GradientReversal import GradientReversalFunction
from utils.data.DNNEventClassifier.GradientReversal import GradientReversal as GR
#from torchsummary import summary

class DNN(nn.Module):
    def __init__(self, n_var=10, n_layers=1, n_nodes=100, n_outputs=2, drop_out_p=0.3):
        super(DNN, self).__init__()
        layers = []
        layers.append(nn.Linear(n_var, n_nodes))
        layers.append(nn.ReLU())

        for n in list(n_nodes for x in range(n_layers)):
            layers.append(nn.Linear(n, n))
            layers.append(nn.ReLU())

        layers.append(nn.Dropout(p=drop_out_p))
        layers.append(nn.Linear(n_nodes, n_outputs))

        self.dnn = nn.Sequential(*layers)

    def forward(self, x):
        return self.dnn(x)

class DNN_GRF(nn.Module):
    def __init__(self, n_var=10, n_layers_features=1, n_layers_tag=1, n_layers_pT=1, n_nodes=100, n_outputs=2, n_pTBins=10, drop_out_p=0.3):
        super(DNN_GRF, self).__init__()

        # Input and feature layers
        self.feature = nn.Sequential()
        self.feature.add_module('i_linear1', nn.Linear(n_var, n_nodes))
        self.feature.add_module('i_relu1',   nn.ReLU())
        self.feature.add_module('i_batchNorm1', nn.BatchNorm1d(n_nodes))
        for i, n in enumerate(list(n_nodes for x in range(n_layers_features))):
            self.feature.add_module('f_linear{}'.format(i+1), nn.Linear(n_nodes, n_nodes))
            self.feature.add_module('f_relu{}'.format(i+1),   nn.ReLU())
            self.feature.add_module('f_batchNorm{}'.format(i+1), nn.BatchNorm1d(n_nodes))

        # Jet tagger classifier
        self.tagger = nn.Sequential()
        for i, n in enumerate(list(n_nodes for x in range(n_layers_tag))):
            self.tagger.add_module('t_linear{}'.format(i+1), nn.Linear(n_nodes, n_nodes))
            self.tagger.add_module('t_relu{}'.format(i+1),   nn.ReLU())
            self.feature.add_module('t_batchNorm{}'.format(i+1), nn.BatchNorm1d(n_nodes))
        self.tagger.add_module('t_dropout',   nn.Dropout(p=drop_out_p))
        self.tagger.add_module('t_linearOut', nn.Linear(n_nodes, n_outputs))

        # Jet pT classifier
        self.pTClass = nn.Sequential()
        # for i, n in enumerate(list(n_nodes for x in range(n_layers_pT))):
        #     self.pTClass.add_module('p_linear{}'.format(i+1), nn.Linear(n_nodes, n_nodes))
        #     self.pTClass.add_module('p_relu{}'.format(i+1),   nn.ReLU())
        # self.pTClass.add_module('p_dropout',   nn.Dropout(p=drop_out_p))
        self.pTClass.add_module('p_linearOut', nn.Linear(n_nodes, n_pTBins))

    def forward(self, input_data, lambdaGR=1.0):
        feature = self.feature(input_data)
        tagger_output = self.tagger(feature)
        # reverse_feature = feature
        reverse_feature = GradientReversalFunction.apply(feature, lambdaGR)
        pTClass_output = self.pTClass(reverse_feature)
        return tagger_output, pTClass_output

if __name__=="__main__":
    # Test CrossEntropyLoss function
    # input is 5 events with 2 outputs
    # target is 5 numbers that correspond to the correct index number for each event
    # Note can not have softmax in model since CrossEntropyLoss does this tep for you. Will need to apply softmax when using the model after training.
    loss = nn.CrossEntropyLoss()
    input = torch.randn(5, 2, requires_grad=True)
    target = torch.empty(5, dtype=torch.long).random_(2)
    output = loss(input, target)
    output.backward()

    # Test DNN model
    # input is 5 events with 10 input variables
    # out is the non-softmax output of the model
    # outFinal is the softmax final version of the model output
    # target is 5 numbers that correspond to the correct index number for each event
    net = DNN()
    loss = nn.CrossEntropyLoss()
    input = torch.randn(5, 10)
    out = net(input)
    outFinal = f.softmax(out,dim=1)
    target = torch.empty(5, dtype=torch.long).random_(2)
    output = loss(out, target)
    output.backward()

    # Print model info
    dnn = DNN(n_var=10, n_layers=1)
    #summary(dnn, (1,10))
    out = dnn(input)
    print(out)

    dnn_grf = DNN_GRF(n_var=10, n_layers=1)
    summary(dnn_grf, [(1,10), (1,1)])
    out, reg = dnn_grf(input_data=input, lambdaGR=1.0)
    print(out, reg)
