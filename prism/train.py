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

from prism import utils
from prism import model
from prism.utils import set_rng_seed
from prism.utils import load_yaml_config
from prism.utils import load_sc_data, load_sc_causal_data

## set random seed
# set_rng_seed(2222)

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
    return args,device
# print(args)

def SetupVAELoss(Model, args, adam_params = None, decayRate = None):
    optimizer = torch.optim.Adam
    if adam_params is None:
        adam_params = {'lr': args['lr'], 'betas':(0.99, 0.999), 'weight_decay': args['weight_decay']} # default
    if decayRate is None:
        decayRate = args['decayrate']
    scheduler = ExponentialLR({'optimizer': optimizer, 'optim_args': adam_params, 'gamma': decayRate})
    pyro.clear_param_store()
    guide = config_enumerate(Model.guide, expand = True)
    elbo = JitTraceEnum_ELBO(max_plate_nesting = 1, strict_enumeration_warning = False)
    loss_basic = SVI(Model.model, guide, scheduler, loss = elbo)
    losses = [loss_basic]
    elbo =  Trace_ELBO()
    loss_aux = SVI(Model.model_GRNrecon, Model.guide_GRNrecon, scheduler, loss = elbo)
    losses.append(loss_aux)
    return losses, scheduler



def Get_metrics(predicted_y, y_prob, y_true):
    correct_prediction = torch.eq(torch.topk(predicted_y, 1)[1].squeeze(), y_true)
    accuracy = torch.mean(correct_prediction.type(torch.FloatTensor))
    AUC = roc_auc_score(y_true.cpu().numpy(), y_prob[:,1].cpu().detach().numpy())
    AUPRC = average_precision_score(y_true.cpu().numpy(), y_prob[:,1].cpu().detach().numpy())
    return accuracy.item(),AUC, AUPRC


def train(arg_path, Expression_data_path, Genescore_data_path, label_path):
    set_rng_seed(2222)  
    args,device = ImportArgs(arg_path)
    
    if args['flag']:
        Eval_acc = torchmetrics.Accuracy(task='multiclass', num_classes = 3).to(device)
        Eval_auc = torchmetrics.AUROC(task='multiclass', num_classes = 3).to(device)
        Eval_ap = torchmetrics.AveragePrecision(task="multiclass", num_classes = 3).to(device)
        onehot_num = 3
    else:
        onehot_num = 2
        
    if args['flag']:
        adj_train, feature, feature_atac, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels = load_sc_causal_data(Expression_data_path, Genescore_data_path, label_path)
    else:
        adj_train, feature, feature_atac, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels = load_sc_data(Expression_data_path, Genescore_data_path, label_path)
    adj_train = F.normalize(adj_train, p=1, dim=1)
    scc = model.PRISM(nfeat=feature.shape[1],     ## the size of feature -> cell num
                    nhid=args['hidden'],         ## hidden layer size
                    dropout=args['dropout'],     ## hyperparameter
                    ns=args['ns'],               ## the size of VAE node embedding 
                    alpha=args['alpha'],         ## hyperparameter
                    flag=args['flag'],
                    use_cuda= args['cuda']).to(device)
    # BCE_loss = torch.nn.BCELoss(reduction='mean')
    # BCEWL_loss = torch.nn.BCEWithLogitsLoss(reduction = 'mean')

    adj_train = adj_train.to(device)
    feature = feature.to(device)
    feature_atac = feature_atac.to(device) # 后续再修改
    train_ids = train_ids.to(device)
    val_ids = val_ids.to(device)
    test_ids = test_ids.to(device)
    train_labels = train_labels.to(device).long()
    val_labels = val_labels.to(device).long()
    test_labels = test_labels.to(device).long()

    losses,scheduler = SetupVAELoss(scc,args)
    train_labels_onehot = F.one_hot(train_labels, onehot_num)
    val_labels_onehot = F.one_hot(val_labels, onehot_num)
    test_labels_onehot = F.one_hot(test_labels, onehot_num)
    loss_reconRNA = []
    loss_reconGRN = []
    best_acc_val = 0

    for i in range(args['epoch']):
        loss1 = losses[0].step(feature, feature_atac, adj_train, train_ids, train_labels_onehot)
        loss_reconRNA.append(loss1)
        loss2 = losses[1].step(feature, feature_atac, adj_train, train_ids, train_labels_onehot)
        loss_reconGRN.append(loss2)
        if (i+1)%100 ==0:
            val_loss1 = losses[0].step(feature, feature_atac, adj_train, val_ids, val_labels_onehot)
            val_loss2 = losses[1].step(feature, feature_atac, adj_train, val_ids, val_labels_onehot)
            val_y, val_y_prob = scc.classifier(feature, adj_train, val_ids)
            if args['flag']:
                val_acc = Eval_acc(val_y_prob , val_labels)
            else:
                val_acc, val_AUC, val_AUPRC = Get_metrics(val_y, val_y_prob,val_labels)
            print(f'On validation epoch {i+1}: RNA recon loss {val_loss1}, GRN recon loss {val_loss2}')
            if best_acc_val <  val_acc:
                best_acc_val = val_acc
            print(f'Validation Accuracy: {best_acc_val}')
    test_y, test_y_prob = scc.classifier(feature, adj_train, test_ids)
    if args['flag']:
        test_acc = Eval_acc(test_y_prob , test_labels).item() 
        test_AUC = Eval_auc(test_y_prob, test_labels).item()
        test_AUPRC = Eval_ap(test_y_prob, test_labels).item()
    else:
        test_acc, test_AUC, test_AUPRC = Get_metrics(test_y, test_y_prob,test_labels)
    print(f'On test set, Accuracy is {test_acc}, AUROC is {test_AUC}, AUPRC is {test_AUPRC}.')


