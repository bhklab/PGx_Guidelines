import torch
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

def Network(args, X):
    IE_dim = X.shape[1]

    class Net1(nn.Module):
        def __init__(self, args):
            super(Net1, self).__init__()

            self.features = torch.nn.Sequential(
                nn.Dropout(args.idr),
                nn.Linear(IE_dim, 512),
                nn.Tanh(),
                nn.Dropout(args.ldr),
                nn.Linear(512, 1))

        def forward(self, x):
            out = self.features(x)
            return out  

    class Net2(nn.Module):
        def __init__(self, args):
            super(Net2, self).__init__()

            self.features = torch.nn.Sequential(
                nn.Dropout(args.idr),
                nn.Linear(IE_dim, 256),
                nn.Tanh(),
                nn.Dropout(args.ldr),
                nn.Linear(256, 256),
                nn.Tanh(),
                nn.Dropout(args.ldr),
                nn.Linear(256, 1)) 

        def forward(self, x):
            out = self.features(x)
            return out            

    class Net3(nn.Module):
        def __init__(self, args):
            super(Net3, self).__init__()    

            self.features = torch.nn.Sequential(
                nn.Dropout(args.idr),
                nn.Linear(IE_dim, 128),
                nn.Tanh(),
                nn.Dropout(args.ldr),
                nn.Linear(128, 128),
                nn.Tanh(),
                nn.Dropout(args.ldr),
                nn.Linear(128, 128),
                nn.Tanh(),
                nn.Dropout(args.ldr),
                nn.Linear(128, 1))                        

        def forward(self, x):
            out = self.features(x)
            return out
    if args.hd == 1:
        Model = Net1(args)
    elif args.hd == 2:
        Model = Net2(args)
    elif args.hd == 3:
        Model = Net3(args)

        
    return Model
