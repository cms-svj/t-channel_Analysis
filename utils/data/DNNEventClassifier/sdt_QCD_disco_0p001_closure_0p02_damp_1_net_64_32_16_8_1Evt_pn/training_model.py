
import torch
import torch.nn as nn



hl1 = 64
hl2 = 32
hl3 = 16
hl4 = 8
features = 42

# Define the neural network
class SimpleNet(nn.Module):
    def __init__(self):
        super(SimpleNet, self).__init__()
        
        self.fc1 = nn.Linear(features, hl1)
        self.batchnorm_x_l1 = nn.BatchNorm1d(hl1)
        
        self.fc2 = nn.Linear(hl1, hl2)
        self.batchnorm_x_l2 = nn.BatchNorm1d(hl2)
        
        self.fc3 = nn.Linear(hl2, hl3)
        self.batchnorm_x_l3 = nn.BatchNorm1d(hl3)
        
        self.fc4 = nn.Linear(hl3, hl4)
        self.batchnorm_x_l4 = nn.BatchNorm1d(hl4)
        
        self.fc5 = nn.Linear(hl4, 1)

    def forward(self, x):
        
        x = self.fc1(x)
        x = self.batchnorm_x_l1(x)
        x = torch.relu(x)

        x = self.fc2(x)
        x = self.batchnorm_x_l2(x)
        x = torch.relu(x)
        
        x = self.fc3(x)
        x = self.batchnorm_x_l3(x)
        x = torch.relu(x)
        
        x = self.fc4(x)
        x = self.batchnorm_x_l4(x)
        x = torch.relu(x)
        
        x = self.fc5(x)
        x = torch.sigmoid(x)
        
        return x

# Create a model instance from the network
model = SimpleNet()


