# PRISM
A Probabilistic Bayesian Model to Recover Gene Regulatory Networks by Incorporating a Biologically Interpretable Structure and Effectively Utilizing Multi-Omics Data
![The framework of PRISM](https://github.com/Ying-Lab/PRISM/blob/main/Figure1.jpg)

Installation
-----

```bash
git clone https://github.com/Ying-Lab/PRISM
cd PRISM
pip install -r requirements.txt 
python setup.py install

```

Example
-----
```bash
from prism import model
from prism import utils

args['flag'] = False
adj_train, feature, feature_ATAC, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels = load_sc_data(Expression_data_path, Genescore_data_path, label_path)
adj_train = F.normalize(adj_train, p=1, dim=1)

scc = model4.BioGRN(nfeat=feature.shape[1],     ## the size of feature -> cell num
                    nhid=args['hidden'],         ## hidden layer size
                    dropout=args['dropout'],     ## hyperparameter
                    ns=args['ns'],               ## the size of VAE node embedding 
                    alpha=args['alpha'],         ## hyperparameter
                    flag=args['flag']).to(device)

```


Citation
-----
