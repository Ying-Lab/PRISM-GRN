from cgi import test
from pyexpat import features
import dgl.data
import torch
import networkx as nx 
import scipy
import scipy.sparse as sp
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import roc_auc_score,average_precision_score
import random
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import LabelBinarizer
import yaml

def set_rng_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True


def load_yaml_config(file_path):
    with open(file_path, encoding='utf-8') as f:
        return yaml.load(f, Loader=yaml.FullLoader)['VAE']

def Get_metrics(predicted_y, y_prob, y_true):
    correct_prediction = torch.eq(torch.topk(predicted_y, 1)[1].squeeze(), y_true)
    accuracy = torch.mean(correct_prediction.type(torch.FloatTensor))
    AUC = roc_auc_score(y_true.cpu().numpy(), y_prob[:,1].cpu().detach().numpy())
    AUPRC = average_precision_score(y_true.cpu().numpy(), y_prob[:,1].cpu().detach().numpy())
    return accuracy.item(),AUC, AUPRC

def Get_metrics_Multi(predicted_y, y_prob, y_true):
    correct_prediction = torch.eq(torch.topk(predicted_y, 1)[1].squeeze(), y_true)
    accuracy = torch.mean(correct_prediction.type(torch.FloatTensor))
    AUC = roc_auc_score(y_true.cpu().numpy(), y_prob[:,1].cpu().detach().numpy())
    AUPRC = average_precision_score(y_true.cpu().numpy(), y_prob[:,1].cpu().detach().numpy())
    return accuracy,AUC, AUPRC

def setup_seed(seed):
     torch.manual_seed(seed)
     torch.cuda.manual_seed_all(seed)
     np.random.seed(seed)
     #random.seed(seed)
     torch.backends.cudnn.deterministic = True


def GetClusters(mat, cm = 'leiden'):
    adata = sc.AnnData(mat, obs= pd.DataFrame(range(mat.shape[0])), var= pd.DataFrame(range(mat.shape[1])))
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    if cm == 'leiden':
        sc.tl.leiden(adata,resolution=1)
    if cm == 'louvain':
        sc.tl.louvain(adata,resolution=1)
    # sc.pl.umap(adata,color='cell type')
    # sc.pl.umap(adata,color='leiden')
    return adata

def GetClusterMetrics(a,b):
    RI = metrics.rand_score(a,b)
    ARI = adjusted_rand_score(a,b)
    AMI = metrics.adjusted_mutual_info_score(a,b)
    NMI = normalized_mutual_info_score(a,b)
    print(f"RI:{RI:.2f}",f"ARI:{ARI:.2f}",f"AMI:{AMI:.2f}",f"NMI:{NMI:.2f}")


def Move2CPU(model):
    model.to('cpu')
    for param in model.parameters():
        param.data = param.data.cpu()
        if param.grad is not None:
            param.grad.data = param.grad.data.cpu()
    return model

def TF_IDF(mat):
    ## cell * peak
    if not scipy.sparse.issparse(mat):
        mat = scipy.sparse.coo_matrix(mat.T)
        
    nfreqs = mat.multiply(1.0 / mat.sum(axis=0))
    tfidf_mat = nfreqs.multiply(np.log(1 + 1.0 * mat.shape[1] / mat.sum(axis=1)).reshape(-1,1)).tocoo()
    
    return tfidf_mat.T


def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(
        np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse.FloatTensor(indices, values, shape).to_dense()

def load_data(name):
    dataset = dgl.data.CoraGraphDataset()
    g = dataset[0]
    
    # feature normalize
    feature = g.ndata['feat']
    label = g.ndata['label']
    # idx_train = g.ndata['train_mask']
    # idx_val = g.ndata['val_mask']
    # idx_test = g.ndata['test_mask']

    # get edges feature
    edges = g.edges()
    g = dgl.graph(edges)
    g = dgl.to_simple(g)
    g = dgl.remove_self_loop(g)
    g = dgl.to_bidirected(g)
    g = g.to_networkx()

    adj = nx.adjacency_matrix(g)
    adj = adj + adj.T.multiply(adj.T > adj) - adj.multiply(adj.T > adj)
    adj = adj + sp.eye(adj.shape[0])
    adj = sparse_mx_to_torch_sparse_tensor(adj)

    idx_train = range(140)
    idx_val = range(200, 500)
    idx_test = range(500, 1500)
    train_mask = np.zeros((adj.size()[0]), dtype=bool)
    val_mask = np.zeros((adj.size()[0]), dtype=bool)
    test_mask = np.zeros((adj.size()[0]), dtype=bool)
    train_mask[idx_train] = True
    val_mask[idx_val] = True
    test_mask[idx_test] = True
    

    idx_train = torch.from_numpy(train_mask)
    idx_val = torch.from_numpy(val_mask)
    idx_test = torch.from_numpy(test_mask)


    return adj, feature, label, idx_train, idx_val, idx_test

def accuracy(output, labels):
    preds = output.max(1)[1].type_as(labels)
    correct = preds.eq(labels).double()
    correct = correct.sum()
    return correct / len(labels)

def accuracy_LP(output, labels):
    output = torch.ge(torch.sigmoid(output), 0.5).type(torch.IntTensor)
    labels = labels.type(torch.IntTensor)
    correct_prediction = torch.eq(output, labels)
    accuracy = torch.mean(correct_prediction.type(torch.FloatTensor))
    return accuracy


def load_sc_data(data_path1, data_path2, label_path):
    label_file = pd.read_csv(label_path, header = 0, sep = ",")
    data = pd.read_csv(data_path1, header = 0, index_col = 0).T ## cell * gene 
    data = data.transform(lambda x: np.log(x + 1))
    data_atac = pd.read_csv(data_path2, header = 0, index_col = 0).T ## cell * gene

    u = []
    v = []
    var_names = list(data.columns) # genes
    ## locate gene index for TF-gene 
    for row_index, row in label_file.iterrows(): 
        u.append(var_names.index(row[0]))
        v.append(var_names.index(row[1]))
    # gene1 index list
    u = np.array(u)
    u = torch.LongTensor(u)
    # gene2 index list
    v = np.array(v)
    v = torch.LongTensor(v)
    ## permutate
    eids = np.arange(label_file.shape[0])
    eids = np.random.permutation(eids)
    # split data
    test_size = int(len(eids) * 0.1)
    val_size = test_size
    train_size = label_file.shape[0] - test_size - val_size
    
    test_pos_u, test_pos_v = u[eids[:test_size]], v[eids[:test_size]]
    val_pos_u, val_pos_v = u[eids[test_size:(test_size + val_size)]], v[eids[test_size:(test_size + val_size)]]
    train_pos_u, train_pos_v = u[eids[(test_size + val_size):]], v[eids[(test_size + val_size):]]

    #find all negative edges and split them for training and testing 
    #the edges those not supported by ChIP-seq are negative
    graph = dgl.graph((u, v), num_nodes = len(var_names))
    
    adj = sp.coo_matrix((np.ones(u.shape), (u, v)),
                        shape=(len(var_names), len(var_names)),
                        dtype=np.float32)
    #adj = graph.adj(scipy_fmt = 'coo')
    adj_neg = 1 - adj.todense() - np.eye(len(var_names))
    neg_u, neg_v = np.where(adj_neg != 0)

    ##For 1:1 Pos-Neg
    # neg_eids = np.random.choice(len(neg_u), label_file.shape[0])
    # 
    ##For 1:All Pos-Neg
    ##  split negative ones
    neg_eids = np.arange(len(neg_u))
    np.random.shuffle(neg_eids)

    test_neg_size = int(len(neg_eids) * 0.1)
    val_neg_size = test_neg_size
    train_neg_size = neg_eids.shape[0] - test_neg_size - val_neg_size

    test_neg_u, test_neg_v = (
        neg_u[neg_eids[:test_neg_size]],
        neg_v[neg_eids[:test_neg_size]],
    )
    val_neg_u, val_neg_v = (
        neg_u[neg_eids[test_neg_size:(test_neg_size + val_neg_size)]],
        neg_v[neg_eids[test_neg_size:(test_neg_size + val_neg_size)]],
    )
    train_neg_u, train_neg_v = (
        neg_u[neg_eids[(test_neg_size + val_neg_size):]],
        neg_v[neg_eids[(test_neg_size + val_neg_size):]],
    )

    train_u = np.concatenate((train_pos_u, train_neg_u), axis = 0)
    train_v = np.concatenate((train_pos_v, train_neg_v), axis = 0)

    val_u = np.concatenate((val_pos_u, val_neg_u), axis = 0)
    val_v = np.concatenate((val_pos_v, val_neg_v), axis = 0)

    test_u = np.concatenate((test_pos_u, test_neg_u), axis = 0)
    test_v = np.concatenate((test_pos_v, test_neg_v), axis = 0)

    #Create train and test mask
    train_ids = np.stack((train_u, train_v), axis = 1)
    train_labels = np.concatenate([np.ones(train_pos_u.shape[0]), np.zeros(train_neg_u.shape[0])], axis = 0)

    val_ids = np.stack((val_u, val_v), axis = 1)
    val_labels = np.concatenate([np.ones(val_pos_u.shape[0]), np.zeros(val_neg_u.shape[0])], axis = 0)

    test_ids = np.stack((test_u, test_v), axis = 1)
    test_labels = np.concatenate([np.ones(test_pos_u.shape[0]), np.zeros(test_neg_u.shape[0])], axis = 0)

    ## the prior graph is spliited corresponding to the splited sets
    train_g = dgl.remove_edges(graph, eids[:(test_size + val_size)])
    # val_g = dgl.graph((val_pos_u, val_pos_v), num_nodes = len(var_names))
    # test_g = dgl.remove_edges(graph, eids[test_size:])
    adj_train = train_g.adj(scipy_fmt = 'coo')


    features = data.T.values ##  gene * cell
    features = normalize_features(features) ## normalize per gene
    features = torch.FloatTensor(features) 

    features_atac = data_atac.T.values ##  gene * cell
    features_atac = normalize_features(features_atac) ## normalize per gene
    features_atac = torch.FloatTensor(features_atac) 

    #Convert the matrix to torch sparse tensor
    adj_train = sparse_mx_to_torch_sparse_tensor(adj_train)
    train_ids = torch.LongTensor(train_ids)
    val_ids = torch.LongTensor(val_ids)
    test_ids = torch.LongTensor(test_ids)
    train_labels = torch.FloatTensor(train_labels)
    val_labels = torch.FloatTensor(val_labels)
    test_labels = torch.FloatTensor(test_labels)

    return adj_train, features, features_atac, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels


def load_sc_data_downsampling(data_path1, data_path2, label_path):
    label_file = pd.read_csv(label_path, header=0, sep=",")
    data = pd.read_csv(data_path1, header=0, index_col=0).T  # cell * gene
    data = data.transform(lambda x: np.log(x + 1))
    data_atac = pd.read_csv(data_path2, header=0, index_col=0).T  # cell * gene

    u = []
    v = []
    var_names = list(data.columns)  # genes
    # locate gene index for TF-gene 
    for row_index, row in label_file.iterrows():
        u.append(var_names.index(row[0]))
        v.append(var_names.index(row[1]))
    # gene1 index list
    u = np.array(u)
    u = torch.LongTensor(u)
    # gene2 index list
    v = np.array(v)
    v = torch.LongTensor(v)
    # permutate
    eids = np.arange(label_file.shape[0])
    eids = np.random.permutation(eids)
    # split data
    test_size = int(len(eids) * 0.1)
    val_size = test_size
    train_size = label_file.shape[0] - test_size - val_size

    test_pos_u, test_pos_v = u[eids[:test_size]], v[eids[:test_size]]
    val_pos_u, val_pos_v = u[eids[test_size:(test_size + val_size)]], v[eids[test_size:(test_size + val_size)]]
    train_pos_u, train_pos_v = u[eids[(test_size + val_size):]], v[eids[(test_size + val_size):]]

    # Create graph and adjacency matrix
    graph = dgl.graph((u, v), num_nodes=len(var_names))
    adj = sp.coo_matrix((np.ones(u.shape), (u, v)), shape=(len(var_names), len(var_names)), dtype=np.float32)

    # Create the negative adjacency matrix (all pairs except the positive edges)
    adj_neg = 1 - adj.todense() - np.eye(len(var_names))
    neg_u, neg_v = np.where(adj_neg != 0)

    # Split negative samples
    neg_eids = np.arange(len(neg_u))
    np.random.shuffle(neg_eids)

    test_neg_size = int(len(neg_eids) * 0.1)
    val_neg_size = test_neg_size
    train_neg_size = neg_eids.shape[0] - test_neg_size - val_neg_size

    # Generate test and validation negative samples (unchanged)
    test_neg_u, test_neg_v = (
        neg_u[neg_eids[:test_neg_size]],
        neg_v[neg_eids[:test_neg_size]],
    )
    val_neg_u, val_neg_v = (
        neg_u[neg_eids[test_neg_size:(test_neg_size + val_neg_size)]],
        neg_v[neg_eids[test_neg_size:(test_neg_size + val_neg_size)]],
    )

    # Train negative samples: downsample to match the number of positive samples
    train_neg_u, train_neg_v = (
        neg_u[neg_eids[(test_neg_size + val_neg_size):]],
        neg_v[neg_eids[(test_neg_size + val_neg_size):]],
    )
    # Downsample negative samples to match the positive sample size
    train_neg_indices = np.random.choice(len(train_neg_u), size=train_pos_u.shape[0], replace=False)
    train_neg_u = train_neg_u[train_neg_indices]
    train_neg_v = train_neg_v[train_neg_indices]

    # Combine positive and negative samples for training, validation, and testing
    train_u = np.concatenate((train_pos_u, train_neg_u), axis=0)
    train_v = np.concatenate((train_pos_v, train_neg_v), axis=0)

    val_u = np.concatenate((val_pos_u, val_neg_u), axis=0)
    val_v = np.concatenate((val_pos_v, val_neg_v), axis=0)

    test_u = np.concatenate((test_pos_u, test_neg_u), axis=0)
    test_v = np.concatenate((test_pos_v, test_neg_v), axis=0)

    # Create train and test mask
    train_ids = np.stack((train_u, train_v), axis=1)
    train_labels = np.concatenate([np.ones(train_pos_u.shape[0]), np.zeros(train_neg_u.shape[0])], axis=0)

    val_ids = np.stack((val_u, val_v), axis=1)
    val_labels = np.concatenate([np.ones(val_pos_u.shape[0]), np.zeros(val_neg_u.shape[0])], axis=0)

    test_ids = np.stack((test_u, test_v), axis=1)
    test_labels = np.concatenate([np.ones(test_pos_u.shape[0]), np.zeros(test_neg_u.shape[0])], axis=0)

    # The prior graph is split corresponding to the split sets
    train_g = dgl.remove_edges(graph, eids[:(test_size + val_size)])
    adj_train = train_g.adj(scipy_fmt='coo')

    features = data.T.values  # gene * cell
    features = normalize_features(features)  # normalize per gene
    features = torch.FloatTensor(features)

    features_atac = data_atac.T.values  # gene * cell
    features_atac = normalize_features(features_atac)  # normalize per gene
    features_atac = torch.FloatTensor(features_atac)

    # Convert the matrix to torch sparse tensor
    adj_train = sparse_mx_to_torch_sparse_tensor(adj_train)
    train_ids = torch.LongTensor(train_ids)
    val_ids = torch.LongTensor(val_ids)
    test_ids = torch.LongTensor(test_ids)
    train_labels = torch.FloatTensor(train_labels)
    val_labels = torch.FloatTensor(val_labels)
    test_labels = torch.FloatTensor(test_labels)

    return adj_train, features, features_atac, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels



def load_sc_causal_data(data_path1, data_path2, label_path):
    if data_path1.split(".")[-1] == "h5":
        store = pd.HDFStore(data_path1)
        data = store['RPKMs']
        store.close()

        ##preprocess the raw expression data
        cellinfo = pd.DataFrame(data.index,index=data.index,columns=['sample_index'])
        geneinfo = pd.DataFrame(data.columns,index=data.columns,columns=['gene_index'])
        adata = sc.AnnData(data.values, obs = cellinfo, var = geneinfo)
        sc.pp.filter_cells(adata, min_genes = 200)
        sc.pp.filter_genes(adata, min_cells = 100)
        data = pd.DataFrame(adata.X, index = adata.obs.index, columns = adata.var.index)
        # 
        ##filter genes using sc_gene_list
        # gene_list = pd.read_csv("~/src_codes/CNNC-master/data/sc_gene_list.txt", sep = "\s+", header = None)
        # data = data[gene_list[1]][:]

        label_file = pd.read_csv(label_path, header = None, sep = "\t")
        print("read data complete!")
    else:
        label_file = pd.read_csv(label_path, header = 0, sep = ",")
        data = pd.read_csv(data_path1, header = 0, index_col = 0).T
        data_atac = pd.read_csv(data_path2, header = 0, index_col = 0).T ## cell * gene
        print("read data complete!")
    
    data = data.transform(lambda x: np.log(x + 1))
    print("log the expression data")

    u = []
    v = []
    d = []
    var_names = list(data.columns)
    for row_index, row in label_file.iterrows():
        u.append(var_names.index(row[0]))
        v.append(var_names.index(row[1]))
        d.append(1)

        u.append(var_names.index(row[1]))
        v.append(var_names.index(row[0]))
        d.append(2)

    print("process the ground truth!")

    u = np.array(u)
    u = torch.LongTensor(u)
    v = np.array(v)
    v = torch.LongTensor(v)
    d = np.array(d)
    d = torch.LongTensor(d)

    eids = np.arange(len(u))
    eids = np.random.permutation(eids)
    test_size = int(len(eids) * 0.1)
    val_size = test_size
    train_size = len(u) - test_size - val_size

    test_pos_u, test_pos_v, test_pos_d = u[eids[:test_size]], v[eids[:test_size]], d[eids[:test_size]]
    val_pos_u, val_pos_v, val_pos_d = u[eids[test_size:(test_size + val_size)]], v[eids[test_size:(test_size + val_size)]], d[eids[test_size:(test_size + val_size)]]
    train_pos_u, train_pos_v, train_pos_d = u[eids[(test_size + val_size):]], v[eids[(test_size + val_size):]], d[eids[(test_size + val_size):]]

    #find all negative edges and split them for training and testing
    graph = dgl.graph((u, v), num_nodes = len(var_names))
    
    adj = sp.coo_matrix((np.ones(u.shape), (u, v)),
                        shape=(len(var_names), len(var_names)),
                        dtype=np.float32)
    #adj = graph.adj(scipy_fmt = 'coo')
    adj_neg = 1 - adj.todense() - np.eye(len(var_names))
    neg_u, neg_v = np.where(adj_neg != 0)
    print("find all negative edges and split them for training and testing!")

    ##For 1:1 Pos-Neg
    # neg_eids = np.random.choice(len(neg_u), label_file.shape[0])
    # 
    ##For 1:All Pos-Neg
    neg_eids = np.arange(len(neg_u))
    np.random.shuffle(neg_eids)

    test_neg_size = int(len(neg_eids) * 0.1)
    val_neg_size = test_neg_size
    train_neg_size = neg_eids.shape[0] - test_neg_size - val_neg_size

    test_neg_u, test_neg_v = (
        neg_u[neg_eids[:test_neg_size]],
        neg_v[neg_eids[:test_neg_size]],
    )
    val_neg_u, val_neg_v = (
        neg_u[neg_eids[test_neg_size:(test_neg_size + val_neg_size)]],
        neg_v[neg_eids[test_neg_size:(test_neg_size + val_neg_size)]],
    )
    train_neg_u, train_neg_v = (
        neg_u[neg_eids[(test_neg_size + val_neg_size):]],
        neg_v[neg_eids[(test_neg_size + val_neg_size):]],
    )

    train_u = np.concatenate((train_pos_u, train_neg_u), axis = 0)
    train_v = np.concatenate((train_pos_v, train_neg_v), axis = 0)

    val_u = np.concatenate((val_pos_u, val_neg_u), axis = 0)
    val_v = np.concatenate((val_pos_v, val_neg_v), axis = 0)

    test_u = np.concatenate((test_pos_u, test_neg_u), axis = 0)
    test_v = np.concatenate((test_pos_v, test_neg_v), axis = 0)

    #Create train and test mask
    train_ids = np.stack((train_u, train_v), axis = 1)
    train_labels = np.concatenate([train_pos_d, np.zeros([train_neg_u.shape[0]])], axis = 0)

    val_ids = np.stack((val_u, val_v), axis = 1)
    val_labels = np.concatenate([val_pos_d, np.zeros([val_neg_u.shape[0]])], axis = 0)

    test_ids = np.stack((test_u, test_v), axis = 1)
    test_labels = np.concatenate([test_pos_d, np.zeros([test_neg_u.shape[0]])], axis = 0)

    
    train_g = dgl.remove_edges(graph, eids[:(test_size + val_size)])
    # val_g = dgl.graph((val_pos_u, val_pos_v), num_nodes = len(var_names))
    # test_g = dgl.remove_edges(graph, eids[test_size:])
    adj_train = train_g.adj(scipy_fmt = 'coo')
    ## gene expression
    features = data.T.values
    features = normalize_features(features)
    features = torch.FloatTensor(features)
    ## atac gene score
    features_atac = data_atac.T.values ##  gene * cell
    features_atac = normalize_features(features_atac) ## normalize per gene
    features_atac = torch.FloatTensor(features_atac) 

    train_eids = np.arange(len(train_labels))
    np.random.shuffle(train_eids)
    train_ids = train_ids[train_eids]
    train_labels = train_labels[train_eids]
    #Convert the matrix to torch sparse tensor
    adj_train = sparse_mx_to_torch_sparse_tensor(adj_train)
    train_ids = torch.LongTensor(train_ids)
    val_ids = torch.LongTensor(val_ids)
    test_ids = torch.LongTensor(test_ids)
    train_labels = torch.LongTensor(train_labels)
    val_labels = torch.LongTensor(val_labels)
    test_labels = torch.LongTensor(test_labels)

    return adj_train, features, features_atac, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels



def load_sc_causal_data_downsampling(data_path1, data_path2, label_path):
    if data_path1.split(".")[-1] == "h5":
        store = pd.HDFStore(data_path1)
        data = store['RPKMs']
        store.close()

        ##preprocess the raw expression data
        cellinfo = pd.DataFrame(data.index,index=data.index,columns=['sample_index'])
        geneinfo = pd.DataFrame(data.columns,index=data.columns,columns=['gene_index'])
        adata = sc.AnnData(data.values, obs = cellinfo, var = geneinfo)
        sc.pp.filter_cells(adata, min_genes = 200)
        sc.pp.filter_genes(adata, min_cells = 100)
        data = pd.DataFrame(adata.X, index = adata.obs.index, columns = adata.var.index)
        # 
        ##filter genes using sc_gene_list
        # gene_list = pd.read_csv("~/src_codes/CNNC-master/data/sc_gene_list.txt", sep = "\s+", header = None)
        # data = data[gene_list[1]][:]

        label_file = pd.read_csv(label_path, header = None, sep = "\t")
        print("read data complete!")
    else:
        label_file = pd.read_csv(label_path, header = 0, sep = ",")
        data = pd.read_csv(data_path1, header = 0, index_col = 0).T
        data_atac = pd.read_csv(data_path2, header = 0, index_col = 0).T ## cell * gene
        print("read data complete!")
    
    data = data.transform(lambda x: np.log(x + 1))
    print("log the expression data")

    u = []
    v = []
    d = []
    var_names = list(data.columns)
    for row_index, row in label_file.iterrows():
        u.append(var_names.index(row[0]))
        v.append(var_names.index(row[1]))
        d.append(1)

        u.append(var_names.index(row[1]))
        v.append(var_names.index(row[0]))
        d.append(2)

    u = np.array(u)
    u = torch.LongTensor(u)
    v = np.array(v)
    v = torch.LongTensor(v)
    d = np.array(d)
    d = torch.LongTensor(d)

    eids = np.arange(len(u))
    eids = np.random.permutation(eids)
    test_size = int(len(eids) * 0.1)
    val_size = test_size
    train_size = len(u) - test_size - val_size

    test_pos_u, test_pos_v, test_pos_d = u[eids[:test_size]], v[eids[:test_size]], d[eids[:test_size]]
    val_pos_u, val_pos_v, val_pos_d = u[eids[test_size:(test_size + val_size)]], v[eids[test_size:(test_size + val_size)]], d[eids[test_size:(test_size + val_size)]]
    train_pos_u, train_pos_v, train_pos_d = u[eids[(test_size + val_size):]], v[eids[(test_size + val_size):]], d[eids[(test_size + val_size):]]

    # Find all negative edges and split them for training and testing
    graph = dgl.graph((u, v), num_nodes=len(var_names))

    adj = sp.coo_matrix((np.ones(u.shape), (u, v)), shape=(len(var_names), len(var_names)), dtype=np.float32)
    adj_neg = 1 - adj.todense() - np.eye(len(var_names))
    neg_u, neg_v = np.where(adj_neg != 0)

    print("Find all negative edges and split them for training and testing!")

    # Split negative samples
    neg_eids = np.arange(len(neg_u))
    np.random.shuffle(neg_eids)

    test_neg_size = int(len(neg_eids) * 0.1)
    val_neg_size = test_neg_size
    train_neg_size = neg_eids.shape[0] - test_neg_size - val_neg_size

    test_neg_u, test_neg_v = (
        neg_u[neg_eids[:test_neg_size]],
        neg_v[neg_eids[:test_neg_size]],
    )
    val_neg_u, val_neg_v = (
        neg_u[neg_eids[test_neg_size:(test_neg_size + val_neg_size)]],
        neg_v[neg_eids[test_neg_size:(test_neg_size + val_neg_size)]],
    )
    train_neg_u, train_neg_v = (
        neg_u[neg_eids[(test_neg_size + val_neg_size):]],
        neg_v[neg_eids[(test_neg_size + val_neg_size):]],
    )

    # Downsample training negative samples to match the number of positive samples
    train_neg_indices = np.random.choice(len(train_neg_u), size=train_pos_u.shape[0], replace=False)
    train_neg_u = train_neg_u[train_neg_indices]
    train_neg_v = train_neg_v[train_neg_indices]

    # Combine positive and downsampled negative samples for training
    train_u = np.concatenate((train_pos_u, train_neg_u), axis=0)
    train_v = np.concatenate((train_pos_v, train_neg_v), axis=0)
    train_labels = np.concatenate([train_pos_d, np.zeros([train_neg_u.shape[0]])], axis=0)

    # Combine positive and negative samples for validation and testing (without downsampling)
    val_u = np.concatenate((val_pos_u, val_neg_u), axis=0)
    val_v = np.concatenate((val_pos_v, val_neg_v), axis=0)
    val_labels = np.concatenate([val_pos_d, np.zeros([val_neg_u.shape[0]])], axis=0)

    test_u = np.concatenate((test_pos_u, test_neg_u), axis=0)
    test_v = np.concatenate((test_pos_v, test_neg_v), axis=0)
    test_labels = np.concatenate([test_pos_d, np.zeros([test_neg_u.shape[0]])], axis=0)

    # Create train and test mask
    train_ids = np.stack((train_u, train_v), axis=1)
    val_ids = np.stack((val_u, val_v), axis=1)
    test_ids = np.stack((test_u, test_v), axis=1)

    train_g = dgl.remove_edges(graph, eids[:(test_size + val_size)])
    adj_train = train_g.adj(scipy_fmt='coo')

    ## gene expression
    features = data.T.values
    features = normalize_features(features)
    features = torch.FloatTensor(features)
    ## atac gene score
    features_atac = data_atac.T.values ##  gene * cell
    features_atac = normalize_features(features_atac) ## normalize per gene
    features_atac = torch.FloatTensor(features_atac) 

    train_eids = np.arange(len(train_labels))
    np.random.shuffle(train_eids)
    train_ids = train_ids[train_eids]
    train_labels = train_labels[train_eids]
    #Convert the matrix to torch sparse tensor
    adj_train = sparse_mx_to_torch_sparse_tensor(adj_train)
    train_ids = torch.LongTensor(train_ids)
    val_ids = torch.LongTensor(val_ids)
    test_ids = torch.LongTensor(test_ids)
    train_labels = torch.LongTensor(train_labels)
    val_labels = torch.LongTensor(val_labels)
    test_labels = torch.LongTensor(test_labels)

    return adj_train, features, features_atac, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels

def normalize_features(mx):
    """Row-normalize sparse matrix"""
    rowsum = np.array(mx.sum(1)).astype(float) 
    r_inv = np.power(rowsum, -1).flatten()
    r_inv[np.isinf(r_inv)] = 0.
    r_mat_inv = sp.diags(r_inv)
    mx = r_mat_inv.dot(mx)
    return mx



def normalize(mx):
    """Row-normalize sparse matrix"""
    rowsum = np.array(mx.sum(1))
    r_inv = np.power(rowsum, -1).flatten()
    r_inv[np.isinf(r_inv)] = 0.
    r_mat_inv = sp.diags(r_inv)
    mx = r_mat_inv.dot(mx)
    return mx

def encode_onehot(labels):
    classes = set(labels)
    classes_dict = {c: np.identity(len(classes))[i, :] for i, c in
                    enumerate(classes)}
    labels_onehot = np.array(list(map(classes_dict.get, labels)),
                             dtype=np.int32)
    return labels_onehot



def plot_confusion_matrix(probabilities, true_labels, class_names):
    """
    绘制混淆矩阵
    :param probabilities: 模型输出的概率，shape 为 (n_samples, n_classes)
    :param true_labels: 真实的标签，shape 为 (n_samples,)
    :param class_names: 类别名称，list of str
    """
    # 1. 从概率中提取预测标签
    # predicted_labels = torch.ge(torch.sigmoid(torch.from_numpy(probabilities)), 0.5).type(torch.IntTensor).numpy
    predicted_labels = np.argmax(probabilities, axis=1)

    # 2. 生成混淆矩阵
    cm = confusion_matrix(true_labels, predicted_labels)

    # 3. 绘制混淆矩阵
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=class_names, yticklabels=class_names)
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title('Confusion Matrix')
    plt.show()


def load_sc_data_clean(data_path1, data_path2, label_path):
    label_file = pd.read_csv(label_path, header = 0, sep = ",")
    data = pd.read_csv(data_path1, header = 0, index_col = 0).T ## transform to cell * gene 
    data = data.transform(lambda x: np.log(x + 1))
    data_atac = pd.read_csv(data_path2, header = 0, index_col = 0).T ## transform to cell * gene

    u = []
    v = []
    var_names = list(data.columns) # genes
    ## locate gene index for TF-gene 
    for row_index, row in label_file.iterrows(): 
        u.append(var_names.index(row[0]))
        v.append(var_names.index(row[1]))
    # gene1 index list
    u = np.array(u)
    u = torch.LongTensor(u)
    # gene2 index list
    v = np.array(v)
    v = torch.LongTensor(v)
    ## permutate
    eids = np.arange(label_file.shape[0])
    eids = np.random.permutation(eids)
    # all prior edges are used for training 
    test_size = 0
    val_size = 0
    train_size = label_file.shape[0] 
    
    test_pos_u, test_pos_v = u[eids[:test_size]], v[eids[:test_size]] # None
    train_pos_u, train_pos_v = u[eids[(test_size + val_size):]], v[eids[(test_size + val_size):]] ## all positive

    #find all negative edges and split them for training and testing 
    #the edges those not supported by ChIP-seq are negative
    graph = dgl.graph((u, v), num_nodes = len(var_names))
    
    adj = sp.coo_matrix((np.ones(u.shape), (u, v)),
                        shape=(len(var_names), len(var_names)),
                        dtype=np.float32)
    #adj = graph.adj(scipy_fmt = 'coo')
    adj_neg = 1 - adj.todense() - np.eye(len(var_names))
    neg_u, neg_v = np.where(adj_neg != 0)

    ##For 1:1 Pos-Neg
    # neg_eids = np.random.choice(len(neg_u), label_file.shape[0])
    # 
    ##For 1:All Pos-Neg
    ##  split negative ones
    neg_eids = np.arange(len(neg_u))
    np.random.shuffle(neg_eids)
    # using 0.1 as negative samples
    # train_neg_size = int(len(neg_eids) * 0.1)
    train_neg_size = train_size 

    test_neg_size = neg_eids.shape[0] - train_neg_size 

    test_neg_u, test_neg_v = (
        neg_u[neg_eids[:test_neg_size]],
        neg_v[neg_eids[:test_neg_size]],
    )

    train_neg_u, train_neg_v = (
        neg_u[neg_eids[test_neg_size:]],
        neg_v[neg_eids[test_neg_size:]],
    )

    train_u = np.concatenate((train_pos_u, train_neg_u), axis = 0)
    train_v = np.concatenate((train_pos_v, train_neg_v), axis = 0)


    test_u = np.concatenate((test_pos_u, test_neg_u), axis = 0)
    test_v = np.concatenate((test_pos_v, test_neg_v), axis = 0)

    #Create train and test mask
    train_ids = np.stack((train_u, train_v), axis = 1)
    train_labels = np.concatenate([np.ones(train_pos_u.shape[0]), np.zeros(train_neg_u.shape[0])], axis = 0)

    test_ids = np.stack((test_u, test_v), axis = 1)

    ## the prior graph is spliited corresponding to the splited sets
    train_g = dgl.remove_edges(graph, eids[:(test_size + val_size)])

    adj_train = train_g.adj(scipy_fmt = 'coo')


    features = data.T.values ##  gene * cell
    features = normalize_features(features) ## normalize per gene
    features = torch.FloatTensor(features) 

    features_atac = data_atac.T.values ##  gene * cell
    features_atac = normalize_features(features_atac) ## normalize per gene
    features_atac = torch.FloatTensor(features_atac) 

    #Convert the matrix to torch sparse tensor
    adj_train = sparse_mx_to_torch_sparse_tensor(adj_train)
    train_ids = torch.LongTensor(train_ids)
    test_ids = torch.LongTensor(test_ids)
    train_labels = torch.FloatTensor(train_labels)

    # return adj_train, features, features_atac, train_ids, val_ids, test_ids, train_labels, val_labels, test_labels
    return adj_train, features, features_atac, train_ids, test_ids, train_labels

def load_sc_causal_data_clean(data_path1, data_path2, label_path):
    if data_path1.split(".")[-1] == "h5":
        store = pd.HDFStore(data_path1)
        data = store['RPKMs']
        store.close()

        ##preprocess the raw expression data
        cellinfo = pd.DataFrame(data.index,index=data.index,columns=['sample_index'])
        geneinfo = pd.DataFrame(data.columns,index=data.columns,columns=['gene_index'])
        adata = sc.AnnData(data.values, obs = cellinfo, var = geneinfo)
        sc.pp.filter_cells(adata, min_genes = 200)
        sc.pp.filter_genes(adata, min_cells = 100)
        data = pd.DataFrame(adata.X, index = adata.obs.index, columns = adata.var.index)
        # 
        ##filter genes using sc_gene_list
        # gene_list = pd.read_csv("~/src_codes/CNNC-master/data/sc_gene_list.txt", sep = "\s+", header = None)
        # data = data[gene_list[1]][:]

        label_file = pd.read_csv(label_path, header = None, sep = "\t")
        print("read data complete!")
    else:
        label_file = pd.read_csv(label_path, header = 0, sep = ",")
        data = pd.read_csv(data_path1, header = 0, index_col = 0).T
        data_atac = pd.read_csv(data_path2, header = 0, index_col = 0).T ## cell * gene
        print("read data complete!")
    
    data = data.transform(lambda x: np.log(x + 1))
    print("log the expression data")

    u = []
    v = []
    d = []
    var_names = list(data.columns)
    for row_index, row in label_file.iterrows():
        u.append(var_names.index(row[0]))
        v.append(var_names.index(row[1]))
        d.append(1)

        u.append(var_names.index(row[1]))
        v.append(var_names.index(row[0]))
        d.append(2)

    print("process the ground truth!")

    u = np.array(u)
    u = torch.LongTensor(u)
    v = np.array(v)
    v = torch.LongTensor(v)
    d = np.array(d)
    d = torch.LongTensor(d)

    eids = np.arange(len(u))
    eids = np.random.permutation(eids)

    # all prior edges are used for training 
    test_size = 0
    val_size = 0
    train_size = label_file.shape[0] 

    test_pos_u, test_pos_v, test_pos_d = u[eids[:test_size]], v[eids[:test_size]], d[eids[:test_size]]
    train_pos_u, train_pos_v, train_pos_d = u[eids[(test_size + val_size):]], v[eids[(test_size + val_size):]], d[eids[(test_size + val_size):]]

    #find all negative edges and split them for training and testing
    graph = dgl.graph((u, v), num_nodes = len(var_names))
    
    adj = sp.coo_matrix((np.ones(u.shape), (u, v)),
                        shape=(len(var_names), len(var_names)),
                        dtype=np.float32)
    #adj = graph.adj(scipy_fmt = 'coo')
    adj_neg = 1 - adj.todense() - np.eye(len(var_names))
    neg_u, neg_v = np.where(adj_neg != 0)
    print("find all negative edges and split them for training and testing!")

    ##For 1:1 Pos-Neg
    # neg_eids = np.random.choice(len(neg_u), label_file.shape[0])
    # 
    ##For 1:All Pos-Neg
    neg_eids = np.arange(len(neg_u))
    np.random.shuffle(neg_eids)
    # balanced 
    train_neg_size = train_size
    test_neg_size = neg_eids.shape[0] - train_neg_size 

    test_neg_u, test_neg_v = (
        neg_u[neg_eids[:test_neg_size]],
        neg_v[neg_eids[:test_neg_size]],
    )

    train_neg_u, train_neg_v = (
        neg_u[neg_eids[test_neg_size:]],
        neg_v[neg_eids[test_neg_size:]],
    )

    train_u = np.concatenate((train_pos_u, train_neg_u), axis = 0)
    train_v = np.concatenate((train_pos_v, train_neg_v), axis = 0)


    test_u = np.concatenate((test_pos_u, test_neg_u), axis = 0)
    test_v = np.concatenate((test_pos_v, test_neg_v), axis = 0)

    #Create train and test mask
    train_ids = np.stack((train_u, train_v), axis = 1)
    train_labels = np.concatenate([train_pos_d, np.zeros([train_neg_u.shape[0]])], axis = 0)
    ## test are all unknown edges, without labels
    test_ids = np.stack((test_u, test_v), axis = 1)


    
    train_g = dgl.remove_edges(graph, eids[:(test_size + val_size)])
    # val_g = dgl.graph((val_pos_u, val_pos_v), num_nodes = len(var_names))
    # test_g = dgl.remove_edges(graph, eids[test_size:])
    adj_train = train_g.adj(scipy_fmt = 'coo')
    ## gene expression
    features = data.T.values
    features = normalize_features(features)
    features = torch.FloatTensor(features)
    ## atac gene score
    features_atac = data_atac.T.values ##  gene * cell
    features_atac = normalize_features(features_atac) ## normalize per gene
    features_atac = torch.FloatTensor(features_atac) 

    train_eids = np.arange(len(train_labels))
    np.random.shuffle(train_eids)
    train_ids = train_ids[train_eids]
    train_labels = train_labels[train_eids]
    #Convert the matrix to torch sparse tensor
    adj_train = sparse_mx_to_torch_sparse_tensor(adj_train)
    train_ids = torch.LongTensor(train_ids)
    val_ids = torch.LongTensor(val_ids)
    test_ids = torch.LongTensor(test_ids)
    train_labels = torch.LongTensor(train_labels)
    val_labels = torch.LongTensor(val_labels)
    test_labels = torch.LongTensor(test_labels)

    return adj_train, features, features_atac, train_ids, test_ids, train_labels

