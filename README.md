# PRISM
A Probabilistic Bayesian Model to Recover Gene Regulatory Networks by a Biologically Interpretable Structure on single-cell Multi-Omics Data
![The framework of PRISM](https://github.com/Ying-Lab/PRISM/blob/main/Figure1.jpg)

Prerequisites
-----

- Python 3.8.10

Installation
-----

```bash
git clone https://github.com/Ying-Lab/PRISM
cd PRISM
pip install -r requirements.txt 
python setup.py install
```
or
```bash
pip install prism-grn
```

Parameters
-----
flag: different tasks, True is causality prediction; False is undirected GRN reconstruction. Default is True.

cuda: Whether using GPU. Default is True.

Paring: Whether the scRNA-seq and scATAC-seq data are paired. Default is True.

epoch: Training epoches. Default is 2000.

lr: Initial learning rate. Default is 0.0003


Example
-----
for GRN Reconstruction
```bash
from prism import model
from prism import utils
from utils import load_sc_data

args['flag'] = False
adj_train, feature, feature_ATAC, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels = load_sc_data(Expression_data_path, Genescore_data_path, label_path)
adj_train = F.normalize(adj_train, p=1, dim=1)

scc = model.PRISM( nfeat=feature.shape[1],     ## the size of feature -> cell num
                    nhid=args['hidden'],         ## hidden layer size
                    dropout=args['dropout'],     
                    ns=args['ns'],               ## the size of VAE node embedding 
                    alpha=args['alpha'],         
                    flag=args['flag'],           ## causal or not
                    use_cuda= args['cuda']).to(device)

```

for Causality prediction
```bash
from prism import model
from prism import utils
from utils import load_sc_causal_data

args['flag'] = True
adj_train, feature, feature_ATAC, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels = load_sc_causal_data(Expression_data_path, Genescore_data_path, label_path)
adj_train = F.normalize(adj_train, p=1, dim=1)

scc = model.PRISM( nfeat=feature.shape[1],     ## the size of feature -> cell num
                    nhid=args['hidden'],         ## hidden layer size
                    dropout=args['dropout'],     
                    ns=args['ns'],               ## the size of VAE node embedding 
                    alpha=args['alpha'],         
                    flag=args['flag'],           ## causal or not
                    use_cuda= args['cuda']).to(device)

```
The more detailed usage is exemplified in demo.


Citation
-----
