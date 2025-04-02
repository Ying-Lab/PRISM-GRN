import sys
import argparse
from lib2to3.pytree import Base
import torch
import numpy as np
import pandas as pd
import random
import yaml
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
import pyro
from pyro.infer import SVI,  JitTraceEnum_ELBO, Trace_ELBO, config_enumerate
from pyro.optim import Adam, ExponentialLR
import torchmetrics
import torch.optim as optim
import torch.nn.functional as F
import torch.utils.data as Data
import scanpy as sc

# import utils
import prism.model as model
from prism.utils import set_rng_seed
from prism.utils import load_yaml_config
from prism.utils import load_sc_data, load_sc_causal_data
from prism.utils import load_sc_data_clean, load_sc_causal_data_clean

## set random seed


## import parameters    
def ImportArgs(arg_path):
    config = load_yaml_config(arg_path) 
    args = { }
    for key in config.keys():
        name = key
        if name not in args:
            args[name] = config[key]
    ## using gpu only when gpu is available and required by users  
    if args['cuda']:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device("cpu")
    ## if gpu is unavailable, using cpu    
    if device == torch.device("cpu"):
        args['cuda'] = False
    args['device'] = device
    return args
# print(args)


def Get_metrics(predicted_y, y_prob, y_true):
    correct_prediction = torch.eq(torch.topk(predicted_y, 1)[1].squeeze(), y_true)
    accuracy = torch.mean(correct_prediction.type(torch.FloatTensor))
    AUC = roc_auc_score(y_true.cpu().numpy(), y_prob[:,1].cpu().detach().numpy())
    AUPRC = average_precision_score(y_true.cpu().numpy(), y_prob[:,1].cpu().detach().numpy())
    return accuracy.item(),AUC, AUPRC


## For evaluation: splitting known GRN into traninig and testing
class Trainer:
    def __init__(self, args, Expression_data_path, Genescore_data_path, label_path):
        self.args = args
        self.device = args['device']
        self.onehot_num = 3 if args['flag'] else 2
        set_rng_seed(args['seed'])


        if args['flag']:
            self.Eval_acc = torchmetrics.Accuracy(task='multiclass', num_classes=3).to(self.device)
            self.Eval_auc = torchmetrics.AUROC(task='multiclass', num_classes=3).to(self.device)
            self.Eval_ap = torchmetrics.AveragePrecision(task='multiclass', num_classes=3).to(self.device)


        self.load_data(Expression_data_path, Genescore_data_path, label_path)
        self.genes = pd.read_csv(Expression_data_path, index_col=0).index.values

        # 归一化
        self.adj_train = F.normalize(self.adj_train, p=1, dim=1).to(self.device)
        self.feature = self.feature.to(self.device)
        self.feature_atac = self.feature_atac.to(self.device)
        if args['pairing']:
            self.scc = model.PRISM(nfeat=self.feature.shape[1], nhid=args['hidden'], dropout=args['dropout'],
                               ns=args['ns'], alpha=args['alpha'], flag=args['flag'], use_cuda=args['cuda']).to(self.device)
        else:
            self.scc = model.PRISM_UP(nfeat=self.feature.shape[1], nfeat_atac=self.feature_atac.shape[1],nhid=args['hidden'], dropout=args['dropout'],
                               ns=args['ns'], alpha=args['alpha'], flag=args['flag'], use_cuda=args['cuda']).to(self.device)
        self.setup_vae_loss()

    def setup_vae_loss(self, adam_params = None, decayRate = None):
        optimizer = torch.optim.Adam
        if adam_params is None:
            adam_params = {'lr': self.args['lr'], 'betas':(0.99, 0.999), 'weight_decay': self.args['weight_decay']} # default
        if decayRate is None:
            decayRate = self.args['decayrate']
        scheduler = ExponentialLR({'optimizer': optimizer, 'optim_args': adam_params, 'gamma': decayRate})
        pyro.clear_param_store()
        guide = config_enumerate(self.scc.guide, expand=True)
        elbo = JitTraceEnum_ELBO(max_plate_nesting=1, strict_enumeration_warning=False)
        loss_basic = SVI(self.scc.model, guide, scheduler, loss=elbo)
        loss_aux = SVI(self.scc.model_GRNrecon, self.scc.guide_GRNrecon, scheduler, loss=Trace_ELBO())
        self.losses = [loss_basic, loss_aux]
        self.scheduler = scheduler

    def load_data(self, Expression_data_path, Genescore_data_path, label_path):
        """加载数据并转换格式"""
        if self.args['flag']:
            self.adj_train, self.feature, self.feature_atac, self.train_ids, self.val_ids, self.test_ids, \
            self.train_labels, self.val_labels, self.test_labels = load_sc_causal_data(Expression_data_path, Genescore_data_path, label_path)
        else:
            self.adj_train, self.feature, self.feature_atac, self.train_ids, self.val_ids, self.test_ids, \
            self.train_labels, self.val_labels, self.test_labels = load_sc_data(Expression_data_path, Genescore_data_path, label_path)

        self.train_labels = self.train_labels.to(self.device).long()
        self.val_labels = self.val_labels.to(self.device).long()
        self.test_labels = self.test_labels.to(self.device).long()

        self.train_labels_onehot = F.one_hot(self.train_labels, self.onehot_num)
        self.val_labels_onehot = F.one_hot(self.val_labels, self.onehot_num)
        self.test_labels_onehot = F.one_hot(self.test_labels, self.onehot_num)

    def train(self):
        """执行训练过程"""
        best_acc_val = 0
        loss_reconRNA, loss_reconGRN = [], []

        for i in range(self.args['epoch']):
            loss1 = self.losses[0].step(self.feature, self.feature_atac, self.adj_train, self.train_ids, self.train_labels_onehot)
            loss2 = self.losses[1].step(self.feature, self.feature_atac, self.adj_train, self.train_ids, self.train_labels_onehot)

            loss_reconRNA.append(loss1)
            loss_reconGRN.append(loss2)

            if (i + 1) % 100 == 0:
                self.validate(i, best_acc_val)

        self.evaluate()
    
    def validate(self, epoch, best_acc_val):
        """执行验证过程"""
        val_loss1 = self.losses[0].step(self.feature, self.feature_atac, self.adj_train, self.val_ids, self.val_labels_onehot)
        val_loss2 = self.losses[1].step(self.feature, self.feature_atac, self.adj_train, self.val_ids, self.val_labels_onehot)
        val_y, val_y_prob = self.scc.classifier(self.feature, self.adj_train, self.val_ids)

        if self.args['flag']:
            val_acc = self.Eval_acc(val_y_prob, self.val_labels)
        else:
            val_acc, val_AUC, val_AUPRC = Get_metrics(val_y, val_y_prob, self.val_labels)

        print(f'On validation epoch {epoch + 1}: RNA recon loss {val_loss1}, GRN recon loss {val_loss2}')

        if best_acc_val < val_acc:
            best_acc_val = val_acc
        print(f'Validation Accuracy: {best_acc_val}')

    def evaluate(self):
        """执行测试评估"""
        with torch.no_grad():
            test_y, test_y_prob = self.scc.classifier(self.feature, self.adj_train, self.test_ids)

            if self.args['flag']:
                test_acc = self.Eval_acc(test_y_prob, self.test_labels).item()
                test_AUC = self.Eval_auc(test_y_prob, self.test_labels).item()
                test_AUPRC = self.Eval_ap(test_y_prob, self.test_labels).item()
            else:
                test_acc, test_AUC, test_AUPRC = Get_metrics(test_y, test_y_prob, self.test_labels)

            print(f'On test set, Accuracy is {test_acc}, AUROC is {test_AUC}, AUPRC is {test_AUPRC}.')

    def Get_GRN(self):
        """执行测试评估"""
        with torch.no_grad():
            test_y, test_y_prob = self.scc.classifier(self.feature, self.adj_train, self.test_ids)
            test_set = pd.DataFrame(self.test_ids.cpu().numpy(),columns=['Gene1','Gene2'])
            test_set['Gene1'] = self.genes[test_set['Gene1'].values]
            test_set['Gene2'] = self.genes[test_set['Gene2'].values]
            if self.args['flag']:  
                test_set['Prob_0'] = test_y_prob.detach().cpu().numpy()[:,0]
                test_set['Prob_1'] = test_y_prob.detach().cpu().numpy()[:,1]
                test_set['Prob_2'] = test_y_prob.detach().cpu().numpy()[:,2]
                test_set['Pred_Label'] = test_y.detach().cpu().numpy()[:,1]                                                                                                                                                       
            else:
                test_set['Prob'] = test_y_prob.detach().cpu().numpy()[:,1]
                test_set['Pred_Label'] = test_y.detach().cpu().numpy()[:,1]
        return test_set


    def save_model(self, path="prism_model.pth"):
        """保存模型"""
        torch.save(self.scc.state_dict(), path)
        print(f"Model saved at {path}")

    def load_model(self, path="prism_model.pth"):
        """加载模型"""
        self.scc.load_state_dict(torch.load(path, map_location=self.device))
        self.scc.to(self.device)
        print(f"Model loaded from {path}")



## All prior GRN are used for training, the test sets are all rest edges
class Trainer_allprior:
    def __init__(self, args, Expression_data_path, Genescore_data_path, label_path):
        self.args = args
        self.device = args['device']
        self.onehot_num = 3 if args['flag'] else 2
        set_rng_seed(args['seed'])


        if args['flag']:
            self.Eval_acc = torchmetrics.Accuracy(task='multiclass', num_classes=3).to(self.device)
            self.Eval_auc = torchmetrics.AUROC(task='multiclass', num_classes=3).to(self.device)
            self.Eval_ap = torchmetrics.AveragePrecision(task='multiclass', num_classes=3).to(self.device)


        self.load_data(Expression_data_path, Genescore_data_path, label_path)
        self.genes = pd.read_csv(Expression_data_path, index_col=0).index.values

        # 归一化
        self.adj_train = F.normalize(self.adj_train, p=1, dim=1).to(self.device)
        self.feature = self.feature.to(self.device)
        self.feature_atac = self.feature_atac.to(self.device)
        if args['pairing']:
            self.scc = model.PRISM(nfeat=self.feature.shape[1], nhid=args['hidden'], dropout=args['dropout'],
                               ns=args['ns'], alpha=args['alpha'], flag=args['flag'], use_cuda=args['cuda']).to(self.device)
        else:
            self.scc = model.PRISM_UP(nfeat=self.feature.shape[1], nfeat_atac=self.feature_atac.shape[1],nhid=args['hidden'], dropout=args['dropout'],
                               ns=args['ns'], alpha=args['alpha'], flag=args['flag'], use_cuda=args['cuda']).to(self.device)
        self.setup_vae_loss()

    def setup_vae_loss(self, adam_params = None, decayRate = None):
        optimizer = torch.optim.Adam
        if adam_params is None:
            adam_params = {'lr': self.args['lr'], 'betas':(0.99, 0.999), 'weight_decay': self.args['weight_decay']} # default
        if decayRate is None:
            decayRate = self.args['decayrate']
        scheduler = ExponentialLR({'optimizer': optimizer, 'optim_args': adam_params, 'gamma': decayRate})
        pyro.clear_param_store()
        guide = config_enumerate(self.scc.guide, expand=True)
        elbo = JitTraceEnum_ELBO(max_plate_nesting=1, strict_enumeration_warning=False)
        loss_basic = SVI(self.scc.model, guide, scheduler, loss=elbo)
        loss_aux = SVI(self.scc.model_GRNrecon, self.scc.guide_GRNrecon, scheduler, loss=Trace_ELBO())
        self.losses = [loss_basic, loss_aux]
        self.scheduler = scheduler

    def load_data(self, Expression_data_path, Genescore_data_path, label_path):
        """加载数据并转换格式"""
        if self.args['flag']:
            self.adj_train, self.feature, self.feature_atac, self.train_ids,  self.test_ids, \
            self.train_labels = load_sc_causal_data_clean(Expression_data_path, Genescore_data_path, label_path)
        else:
            self.adj_train, self.feature, self.feature_atac, self.train_ids, self.test_ids, \
            self.train_labels = load_sc_data_clean(Expression_data_path, Genescore_data_path, label_path)

        self.train_labels = self.train_labels.to(self.device).long()
        self.train_labels_onehot = F.one_hot(self.train_labels, self.onehot_num)

    def train(self):
        """执行训练过程"""
        best_acc_val = 0
        loss_reconRNA, loss_reconGRN = [], []

        for i in range(self.args['epoch']):
            loss1 = self.losses[0].step(self.feature, self.feature_atac, self.adj_train, self.train_ids, self.train_labels_onehot)
            loss2 = self.losses[1].step(self.feature, self.feature_atac, self.adj_train, self.train_ids, self.train_labels_onehot)

            loss_reconRNA.append(loss1)
            loss_reconGRN.append(loss2)

            if (i + 1) % 100 == 0:
                print(f'On tranining epoch {i + 1}: RNA recon loss {loss1}, GRN recon loss {loss2}')
    

    def Get_GRN(self):
        """执行测试评估"""
        with torch.no_grad():
            test_y, test_y_prob = self.scc.classifier(self.feature, self.adj_train, self.test_ids)
            test_set = pd.DataFrame(self.test_ids.cpu().numpy(),columns=['Gene1','Gene2'])
            test_set['Gene1'] = self.genes[test_set['Gene1'].values]
            test_set['Gene2'] = self.genes[test_set['Gene2'].values]
            if self.args['flag']:  
                test_set['Prob_0'] = test_y_prob.detach().cpu().numpy()[:,0]
                test_set['Prob_1'] = test_y_prob.detach().cpu().numpy()[:,1]
                test_set['Prob_2'] = test_y_prob.detach().cpu().numpy()[:,2]
                test_set['Pred_Label'] = test_y.detach().cpu().numpy()[:,1]                                                                                                                                                       
            else:
                test_set['Prob'] = test_y_prob.detach().cpu().numpy()[:,1]
                test_set['Pred_Label'] = test_y.detach().cpu().numpy()[:,1]
        return test_set


    def save_model(self, path="prism_model.pth"):
        """保存模型"""
        torch.save(self.scc.state_dict(), path)
        print(f"Model saved at {path}")

    def load_model(self, path="prism_model.pth"):
        """加载模型"""
        self.scc.load_state_dict(torch.load(path, map_location=self.device))
        self.scc.to(self.device)
        print(f"Model loaded from {path}")