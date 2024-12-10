#!/root/miniconda3/bin/python3
import sys
print(sys.version)
import pandas as pd
import numpy as np
import pickle
import torch.utils.data as data_utils
from sklearn.preprocessing import LabelEncoder
import torch
from torch import nn
from torch import optim
import os
from torch.autograd import Variable
from torch.utils.data import DataLoader
import torch.nn.functional as F
import warnings

class NN_Model(nn.Module):
    # Model
    def __init__(self, n_input, n_output):
        super(NN_Model, self).__init__()
        self.layer_out = nn.Linear(n_input, n_output, bias=False)

    def forward(self, x):
        x = self.layer_out(x)
        return x

class NN_classifier():
    def __init__(self, model_files_path):
        # Load model and auxiliary files from model_files_path
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        with open(model_files_path,'rb') as f:
            model_files = pickle.load(f)
            self.model = model_files[0]
            self.enc = model_files[1]
            self.input_features = model_files[2]

        self.DM = NN_Model(self.model[list(self.model)[-1]].size()[1], self.model[list(self.model)[-1]].size()[0])
        self.DM.load_state_dict(self.model)
        self.DM.to(self.device)

    def predict(self, sample):
        ''' Make prediction based on binarized sample vector
        Predicts methylation class from a binarized vector of CpG sites

        Argument:
        sample: a data frame with columns probe_id and methylation_call:
                probe_id        methylation_call
        0       cg25742403      1.0
        1       cg04856205      1.0

        Returns:
         - scores for each class in descending order, i.e. majority vote first
         - class labels in descending order
         - number of features overlapping model
        '''

        input_dnn = self.input_features.merge(sample,how='left')
        input_dnn['methylation_call']=input_dnn['methylation_call'].fillna(0)
        torch_tensor = torch.tensor(input_dnn['methylation_call'].values)
        y_val_pred_masked = self.DM(torch_tensor.float().to(self.device))

        predictions = torch.softmax( (y_val_pred_masked - y_val_pred_masked.mean().item())/y_val_pred_masked.std(unbiased=False).item(), dim = 0)
        class_labels = self.enc.inverse_transform(torch.topk(predictions, len(predictions)).indices.tolist()).tolist()

        return torch.topk(predictions, len(predictions)).values.tolist(), class_labels, len(input_dnn[input_dnn['methylation_call']!=0])

    def predict_from_bedMethyl(self, path):
        '''Preprocess annotated bedMethyl file specified by path (name column should hold Illumina probe ID) and predict'''
        bed = pd.read_csv(path, delim_whitespace=True, header=None, usecols=range(11)) # read BED9+2
        bed.columns = ['chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','coverage','MAF']

        sample = bed.rename(columns = {"name": "probe_id", "MAF": "methylation_call"})[['probe_id', 'methylation_call']]

        sample.loc[sample["methylation_call"] < 60, "methylation_call"] = -1
        sample.loc[sample["methylation_call"] >= 60, "methylation_call"] = 1

        return self.predict(sample)
