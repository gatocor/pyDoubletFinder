import scanpy as scp
import numpy as np
import scipy as sp
import pandas as pd
from sklearn.neighbors import KernelDensity

def bimodality_coefficient(x):
    g = sp.stats.skew(x)
    k = sp.stats.kurtosis(x)
    n = np.size(x)
    b = ((g**2)+1)/(k+((3*((n-1)**2))/((n-2)*(n-3))))

    return b

def paramSweep_and_summarize(adata,PCs,
                metric = "euclidean",
                pNscreening=np.arange(0.05,0.35,0.05),
                pKscreening=np.append([0.0005,0.001,0.005],np.arange(0.01,0.31,0.01))
                ):
    
    # Extract relevant parameters
    n_real_cells = adata.shape[0]
    n_doublets = round(n_real_cells/(1-pNscreening[-1]) - n_real_cells)

    # Generate Doublets
    sel1 = np.random.randint(0,n_real_cells,size=n_doublets)
    sel2 = np.random.randint(0,n_real_cells,size=n_doublets)
    XDoublets = adata.X[sel1,:]+adata.X[sel2,:]

    # Screen over pN
    data = pd.DataFrame()
    for pN in pNscreening:
        #Compute doublets number
        n_doublets = round(n_real_cells/(1-pN) - n_real_cells)

        # Make Ann object
        X = adata.X.copy()
        X = sp.sparse.vstack((X,XDoublets[:n_doublets,:]))

        adataC = scp.AnnData(X)
        adataC.obs["Doublet"] = True
        adataC.obs.loc[:,"Doublet"].iloc[0:n_real_cells] = False

        # Normalize
        scp.pp.log1p(adataC)

        # Scale
        scp.pp.normalize_total(adataC,target_sum=100)

        # Find highly variable genes
        scp.pp.highly_variable_genes(adataC,flavor="seurat")

        # PCA
        scp.pp.pca(adataC,n_comps=PCs)

        # Screen over pK
        for pK in pKscreening:
            # Compute distance matrix
            k = max(round(adataC.shape[0]*pK),2)
            scp.pp.neighbors(adataC,n_neighbors=k,metric=metric)

            # Compute pANN
            l1,l2 = adataC.obsp["distances"].nonzero()
            m = pd.DataFrame()
            m["Doublet"] = adataC.obs.iloc[l2,-1]
            m["Index"] = l1
            m = m.groupby("Index").mean().values[:n_real_cells,0]

            # Smooth with gaussian kernel
            model = KernelDensity()
            model.fit(m.reshape(-1,1))
            x = model.score_samples(np.arange(0,1,0.01).reshape(-1,1))
            data.loc[pN,pK] = bimodality_coefficient(x)

    return data

def doubletFinder(adata, PCs, pK, nExp, pN = 0.25, metric="euclidean"):

    # Extract relevant parameters
    n_real_cells = adata.shape[0]
    n_doublets = round(n_real_cells/(1-pN) - n_real_cells)

    # Generate Doublets
    sel1 = np.random.randint(0,n_real_cells,size=n_doublets)
    sel2 = np.random.randint(0,n_real_cells,size=n_doublets)
    XDoublets = adata.X[sel1,:]+adata.X[sel2,:]

    # Make Ann object
    X = adata.X.copy()
    X = sp.sparse.vstack((X,XDoublets))

    adataC = scp.AnnData(X)
    adataC.obs["Doublet"] = True
    adataC.obs.loc[:,"Doublet"].iloc[0:n_real_cells] = False

    # Normalize
    scp.pp.log1p(adataC)

    # Scale
    scp.pp.normalize_total(adataC,target_sum=100)

    # Find highly variable genes
    scp.pp.highly_variable_genes(adataC,flavor="seurat")

    # PCA
    scp.pp.pca(adataC,n_comps=PCs)

    # Compute distance matrix
    k = round(adataC.shape[0]*pK)
    scp.pp.neighbors(adataC,n_neighbors=k,metric=metric)

    # Compute pANN
    l1,l2 = adataC.obsp["distances"].nonzero()
    m = pd.DataFrame()
    m["Doublet"] = adataC.obs.iloc[l2,-1]
    m["Index"] = l1
    m = m.groupby("Index").mean()
    adataC.obs["DoubletFinder_score"] = m.values

    # Classify doublets
    cut = np.sort(m.values[:,0])[-nExp]
    adataC.obs["DoubletFinder_imputed"] = m.values >= cut

    adata.obs

    return adataC