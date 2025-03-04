import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import pyro
import pyro.distributions as dist
import torch.nn.functional as F
import torch
import torch.nn as nn
# import prism
from prism.layer import MLP, Exp, ExpM
from prism.layer import GraphConvolution, GraphAttentionLayer
torch.set_default_tensor_type(torch.FloatTensor)


# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
## Graph VAE encoder 
class Encoder(nn.Module):
    def __init__(self, nhid, nz) -> None:
        super().__init__()
        self.linear1 = nn.Linear(nhid, nz)
        self.linear2 = nn.Linear(nhid, nz)

    
    def forward(self, x):
        #x = torch.cat((x, y), dim=1)
        mu = F.relu(self.linear1(x))
        logvar = F.relu(self.linear2(x))
        return mu, logvar

## decode gene expression   
class Decoder(nn.Module):
    def __init__(self, ns, nhid) -> None:
        super().__init__()
        self.linear_x = nn.Linear(ns, nhid)
    
    def forward(self, s):
        x = F.sigmoid(self.linear_x(s))
        return x


## only return the embeddings of GAT   
class GATEncoder(nn.Module):
    def __init__(self, nfeat, nhid, nz, ns,dropout, alpha) -> None:
        super().__init__()
        self.dropout = dropout
        # concat: whether input elu layer
        # encoder

        self.attention_z = GraphAttentionLayer(nfeat, nhid, dropout=dropout, alpha=alpha, concat=True)
        self.encoder = Encoder(nhid, nz)  ## nhid -> nz

        # z->s 
        self.attention_s = GraphAttentionLayer(nfeat, nhid, dropout=dropout, alpha=alpha, concat=True)

        self.linear_s = nn.Linear(nhid + nz, ns)  ## nhid + nz -> ns
        
        self.nz = nz


    def forward(self, feat, adj, stage=None, z_mean=None):
        # x = Wx
        # gene * nhid
        h_z = self.attention_z(feat, adj) 
        h_s = self.attention_s(feat, adj)
        # x = F.dropout(h_z, self.dropout, training=self.training)
        x = h_z
        mu = [] 
        logvar = []
        z_sum = 0
        if stage == 'training':
        # encoder 
            mu, logvar = self.encoder(x)
            z = self.reparametrize(mu, logvar)
            z_sum = z.var(dim=1).sum() ## the var of each gene and summed up
            z_mean = torch.mean(z, dim=0).unsqueeze(0)  # 1 * nz the mean across all genes on each dimension 
        else: 
            z = z_mean.repeat(x.size()[0], 1) # gene * nz

        # z-> s
        x = h_s
        ## concat VAE & GAT
        x = torch.cat((z, x), dim=1)
        ## diffusion gene * ns 
        s = F.relu(self.linear_s(x))

        # predict GRN 
        # output = self.base_model(s, recon_x, adj, train_ids)
        ## only return the gene emeddings ~ Gaussian
        ## calculate the KL
        return z_mean, mu, logvar, z_sum, s
    
    def reparametrize(self, mu, logvar):

        eps = torch.randn_like(logvar)
        z = mu + eps * torch.exp(logvar/2)
        return z


class GCN(nn.Module):
    def __init__(self, nfeat, ns, nhid, dropout, flag = False):
        super(GCN, self).__init__()

        self.gc1 = GraphConvolution(nfeat, nhid)
        self.gc2 = GraphConvolution(nhid, nhid)
        self.linear1 = nn.Linear((nhid+ns) * 2, nhid)
        if flag:
            self.linear2 = nn.Linear(nhid, 3)
        else:
            self.linear2 = nn.Linear(nhid, 1)
        #self.dotprodcut = DotProductPredictor(dropout, act = lambda x: x)
        self.flag = flag
        self.dropout = dropout

    def forward(self, s, x, adj, train_ids):
        # tranditional gcn
        x = F.dropout(x, self.dropout, training=self.training)

        x = F.relu(self.gc1(x, adj))
        x = F.dropout(x, self.dropout, training=self.training)
        # x = F.relu(self.gc2(x, adj))
        x = self.gc2(x, adj)


        # add s
        # x = F.dropout(x, self.dropout, training=self.training)
        x = torch.cat((x, s), dim=1)
        x_src = x[train_ids[:, 0]]
        x_dst = x[train_ids[:, 1]]

        x_edges = torch.cat([x_src, x_dst], dim=1)
        output = F.relu(self.linear1(x_edges))
        output = self.linear2(output)
        # x = self.dotprodcut(x)
        if self.flag:
            return output
        else:
            return output.squeeze(1)

class GAT(nn.Module):
    def __init__(self, nfeat, nhid, dropout, alpha, nheads, ns, flag = False):
        """Dense version of GAT."""
        super(GAT, self).__init__()
        self.dropout = dropout

        self.attentions = [GraphAttentionLayer(nfeat, nhid, dropout=dropout, alpha=alpha, concat=True) for _ in range(nheads)]
        for i, attention in enumerate(self.attentions):
            self.add_module('attention_{}'.format(i), attention)

        self.out_att = GraphAttentionLayer(nhid * nheads, nhid * nheads, dropout=dropout, alpha=alpha, concat=False) 

        self.linear1 = nn.Linear((nhid * nheads + ns) * 2, nhid * nheads)
        if flag:
            self.linear2 = nn.Linear(nhid * nheads, 3)
        else:
            self.linear2 = nn.Linear(nhid * nheads, 1)
        self.flag = flag

    def forward(self, s, x, adj, train_ids):
        # tranditional GAT
        x = F.dropout(x, self.dropout, training=self.training)
        x = torch.cat([att(x, adj) for att in self.attentions], dim=1)
        # x = F.dropout(x, self.dropout, training=self.training)
        # x = F.elu(self.out_att(x, adj))
        x = self.out_att(x, adj)
        x = F.dropout(x, self.dropout, training=self.training)

        # # add s
        x = torch.cat((x, s), dim=1)
        x_src = x[train_ids[:, 0]]
        x_dst = x[train_ids[:, 1]]

        x_edges = torch.cat([x_src, x_dst], dim=1)
        output = F.relu(self.linear1(x_edges))
        if self.flag:
            output = self.linear2(output)
        else:
            output = self.linear2(output).squeeze(1)
        return output


class DotProductPredictor(nn.Module):
	def __init__(self, dropout = 0., act = F.sigmoid) -> None:
		super().__init__()
		self.dropout = dropout
		self.act = act

	def forward(self, h):
		h = F.dropout(h, self.dropout)
		x = torch.mm(h, h.T)
		return self.act(x)
    


##define the entity model
# gene * cell  
# nfeat, nhid, nz, ns,dropout, alpha, base_model='gcn', nheads=8, flag=False
class PRISM(nn.Module):
    def __init__(self,
                 nfeat = None,  ## cell num
                 nhid = None,   ## hidden size
                #  nTF = None,    ## gene embedding size -> TF num
                #  nz = None,     ## GCN gene embedding size
                 ns = None,     ## GCN -> VAE size
                #  hidden_layers = [50],
                ## default setting
                 dropout = 0.1,
                 alpha = 0.35,
                 flag = True, 
                 config_enum = "parallel",
                 use_cuda = True,
                 aux_loss_multiplier = 1, ## hyperparameter
    ):
        super().__init__()    
        ## input parameters
        self.nfeat = nfeat
        self.nhid = nhid  ## VAE hidden TF_number
        self.ns = ns    ## gene embeddung

        ##defalt parameters
        self.use_cuda = use_cuda
        ## if there is no gpu, device is cpu
        if self.use_cuda:
            self.device = torch.device("cuda:0")
        else:
            self.device = torch.device("cpu")
        self.allow_broadcast = config_enum 
        self.aux_loss_multiplier = aux_loss_multiplier
        ## predicting GRN edge based on gene embeddings
        if flag:
            self.ny = 3
        else:
            self.ny = 2
        self.edgeLinear = nn.Linear(self.nhid, self.ny)
        ## GAT for gene embedding 
        self.linear_s = nn.Linear(nhid + ns, ns)
        self.GeneGRNEncoder = GraphAttentionLayer(nfeat, nhid, dropout=dropout, alpha=alpha, concat=True)
        self.GeneGRNEncoder_VAE = GraphAttentionLayer(nfeat, nhid, dropout=dropout, alpha=alpha, concat=True) ## revised later 
        ## Recon GRN 
        #choice 1: directly based on the trained embeddings from GCN 
        #if this choice, deposit in the VAE networks 
        #choice 2：create another GCN/GAT based on the graph and the features(trained/initial) test later
        self.GRNDecoder = GCN(nfeat, ns, nhid, dropout, flag)
        # self.GRNDecoder = GAT(nfeat, nhid, dropout, alpha, nheads, ns, flag)


        ## VAE encoders/decoders
        self.setup_networks()
    

    ## VAE networks
    def setup_networks(self):
        ## encoders
        ## rp score -> hidden -> gene embedding
        self.encoder_ATAC = MLP(
            [self.nfeat] + [self.nhid] + [[self.ns, self.ns]],
            activation = nn.Softplus,
            output_activation = [None, Exp],
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )

        ## rna -> hidden -> gene embedding 
        self.encoder_RNA = MLP(
            [self.nfeat] + [self.nhid] + [[self.ns, self.ns]],
            activation = nn.Softplus,
            output_activation = [None, Exp],
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )
        ## GATgene encoder -> hidden -> GATVAE embedding
        self.encoder_GAT = MLP(
            [self.nhid] + [self.nhid] + [[self.ns, self.ns]],
            activation = nn.Softplus,
            output_activation = [None, Exp],
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )
        
        ## directly produce y based on the paired gene embeddings (2 * ns) VAE-based ONEHOT CATEGORIAL DIST 

        self.encoder_GRNedges = MLP(
            [self.ns * 2] + [self.nhid] + [self.ny],
            activation = nn.Softplus,
            output_activation = nn.Softmax,
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )
        
        ## decoders
        ## expression recon
        ## expression gene embedding + rp gene embedding + GCN embedding -> expression 
        self.decoder_RNA = MLP( 
            [self.ns + self.ns + self.ns ] + [self.nhid] + [self.nfeat],
            activation = nn.Softplus,
            output_activation = nn.Sigmoid,
            allow_broadcast = self.allow_broadcast,
            use_cuda = self.use_cuda,
        )

        self.cutoff = nn.Threshold(1.0e-9, 1.0e-9)
        if self.use_cuda:
            self.to(self.device)
    
   
    def model(self, XRNA, XATAC = None, adj = None, train_ids = None,train_y = None):
        """
        The model corresponds to the following generative process:		
	    the prior distribution
        θ= decoder_θ(zrna, zatac, zgrn)	
        X ~ Multinomial(θ)  
        :return: None
        """
        pyro.module('RNArecon', self)
        XRNA = XRNA.to(self.device)   
        batch_size = XRNA.size(0) ## all genes
        options = dict(dtype=XRNA.dtype, device=XRNA.device)
        XGATembedding = self.GeneGRNEncoder(XRNA, adj) ## revise later, test whether to delete this variable # nhid

        with pyro.plate('data'):
            # prior initialize
            ## zrna
            zrna_loc = torch.zeros(batch_size, self.ns, **options)
            zrna_scale = torch.ones(batch_size, self.ns, **options)
            zrna = pyro.sample('zrna', dist.Normal(zrna_loc , zrna_scale).to_event(1)).to(self.device)
            ## zatac 
            zatac_loc = torch.zeros(batch_size, self.ns, **options)
            zatac_scale = torch.ones(batch_size, self.ns, **options)
            zatac = pyro.sample('zatac', dist.Normal(zatac_loc, zatac_scale).to_event(1)).to(self.device)
            # ## zgrn based on GATGeneEncoder #nhid -> ns 
            zgrn_loc = torch.zeros(batch_size, self.ns, **options)
            zgrn_scale = torch.ones(batch_size, self.ns, **options)
            ## zgrn 
            zgrn = pyro.sample('zgrn', dist.Normal(zgrn_loc, zgrn_scale).to_event(1)).to(self.device)
            # if need the diffusion layer， uncomment the next two lines 
            zgrn = torch.cat((XGATembedding, zgrn), dim=1)
            zgrn = F.relu(self.linear_s(zgrn))
            ## zgrn for directly encode by GRNEncoder
            # zgrn is the embeddings after GAT and processed by VAE 
            ##rna decoder
            # thetas = self.decoder_RNA ([zrna, zatac, s]).to(device)
            thetas = self.decoder_RNA ([zrna, zatac, zgrn]).to(self.device)
            thetas = self.cutoff(thetas)
            max_count = torch.ceil(abs(XRNA).sum(1).sum()).int().item()
            pyro.sample('XRNA', dist.DirichletMultinomial(total_count = max_count, concentration=thetas), obs=XRNA)
            
        # ## GRN decoder
        # recon_GRN = self.GRNDecoder(zgrn, XRNA, adj, train_ids)
        # return z_mean, mu, logvar, z_sum, recon_GRN 

    def guide(self, XRNA, XATAC, adj = None,train_ids = None,train_y = None):
        ## post-distribution of VAE
        XRNA = XRNA.to(self.device)
        XATAC= XATAC.to(self.device)
        adj = adj.to(self.device)
        # XGATembedding = self.GeneGRNEncoder(XRNA, adj) # 共享 GAT 的 parameter
        XGATembedding = self.GeneGRNEncoder_VAE(XRNA, adj).to(self.device) #不共享 GAT 的 parameter
        with pyro.plate('data'):
            ## observed 
            ##rna -> gene_zran
            zrna_loc, zrna_scale = self.encoder_RNA(XRNA)
            zrna = pyro.sample('zrna', dist.Normal(zrna_loc, zrna_scale).to_event(1)).to(self.device)

            ##atac -> gene_zatac
            zatac_loc, zatac_scale = self.encoder_ATAC(XATAC)
            zatac = pyro.sample('zatac', dist.Normal(zatac_loc, zatac_scale).to_event(1)).to(self.device)

            ##GATenmbedding -> gene_zgrn
            ##atac -> zatac
            zgrn_loc, zgrn_scale = self.encoder_GAT(XGATembedding)
            zgrn = pyro.sample('zgrn', dist.Normal(zgrn_loc, zgrn_scale).to_event(1)).to(self.device)

    # def forward(self, XRNA, XATAC, adj, train_ids, stage=None, z_mean=None):
    #     z_mean, mu, logvar, z_sum, zgrn = self.GRNEncoder(XRNA, adj, stage)
    #     recon_GRN = self.GRNDecoder(zgrn, XRNA, adj, train_ids)
    #     return z_mean, mu, logvar, z_sum, recon_GRN

    def classifier(self, XRNA, adj, train_ids):
        """
        classify a cell (or a batch of cells)

        :param xs: a batch of vectors of gene counts from a cell
        :return: a batch of the corresponding class labels (as one-hots)
                 along with the class probabilities
        """
        # use the trained model q(y|x) = categorical(alpha(x))
        # compute all class probabilities for the cell(s)
        XRNA = XRNA.to(self.device)
        adj = adj.to(self.device)
        XGATembedding = self.GeneGRNEncoder(XRNA, adj).to(self.device) 
        XGATembedding2 = self.GeneGRNEncoder_VAE(XRNA, adj).to(self.device)
        GeneEmbedding_GAT,_ = self.encoder_GAT(XGATembedding2)
        GeneEmbedding_GVAE = torch.cat((XGATembedding, GeneEmbedding_GAT), dim=1)
        GeneEmbedding_GVAE = F.relu(self.linear_s(GeneEmbedding_GVAE))

        Gene1 = GeneEmbedding_GVAE[train_ids[:, 0]]
        Gene2 = GeneEmbedding_GVAE[train_ids[:, 1]]
        X_edges = torch.cat([Gene1, Gene2], dim=1)
        y_alpha = self.encoder_GRNedges(X_edges).to(self.device)
        res, ind = torch.topk(y_alpha, 1)
        edge_y = torch.zeros_like(y_alpha).scatter_(1, ind, 1.0).to(self.device)
        return edge_y, y_alpha
    
    ## GRN recon
    def model_GRNrecon(self, XRNA, XATAC = None, adj = None, train_ids = None, train_y = None):
        """
        this model is used to add auxiliary (supervised) loss as described in the
        Kingma et al., "Supervised Learning with Deep Generative Models".
        """
        # register all pytorch (sub)modules with pyro
        pyro.module('GRNrecon', self)
    
        # inform pyro that the variables in the batch of xs, ys are conditionally independent
        with pyro.plate('data'):
            # this here is the extra term to yield an auxiliary loss that we do gradient descent on
                XGATembedding = self.GeneGRNEncoder(XRNA, adj).to(self.device) ## GAT
                XGATembedding2 = self.GeneGRNEncoder_VAE(XRNA, adj).to(self.device)
                GeneEmbedding_GAT,_ = self.encoder_GAT(XGATembedding2) ## VAE_embedding based on GAT, using mean values instead of sampling
                GeneEmbedding_GVAE = torch.cat((XGATembedding, GeneEmbedding_GAT), dim=1)
                GeneEmbedding_GVAE = F.relu(self.linear_s(GeneEmbedding_GVAE)).to(self.device)
                Gene1 = GeneEmbedding_GVAE[train_ids[:, 0]]
                Gene2 = GeneEmbedding_GVAE[train_ids[:, 1]]
                X_edges = torch.cat([Gene1, Gene2], dim=1)
                alpha_y = self.encoder_GRNedges(X_edges).to(self.device)
                ## need the real label
                with pyro.poutine.scale(scale = self.aux_loss_multiplier):
                    pyro.sample('train_y', dist.OneHotCategorical(alpha_y), obs = train_y)
            

    def guide_GRNrecon(self, XRNA, XATAC = None, adj = None, train_ids = None, train_y = None):
        
        """
        Additional classification and none inference steps are included
        dummy guide function to accompany model_classify in inference
        """
        pass