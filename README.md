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
import scibd as si
import scanpy as sc
dataset = sc.read_h5ad(./data/Forebrain.h5ad)
KNNITER = si.KNNIter(dataset,strategy = 'PCoA')
result = KNNITER.IterCall()

```


Citation
-----
