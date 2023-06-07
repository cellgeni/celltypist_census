#!/nfs/cellgeni/batu/conda-envs/scvi-env/bin/python

import numpy as np
import pandas as pd
import scanpy as sc
import scvi

def myRange(start,end,step):
    i = start
    while i < end:
        yield i
        i += step
    yield end

adata = sc.read('/lustre/scratch127/cellgen/cellgeni/shibla/cell-census/bone-marrow_cell-census.h5ad')

rangelist = [i for i in (myRange(0,adata.shape[0],10000))]

vae = scvi.model.SCVI.load('/lustre/scratch127/cellgen/cellgeni/cakirb/census/', prefix = 'bonemarrow_scvi', adata = adata)

for i in range(0,len(rangelist)+1):
    if i == 61:
        scales0.to_pickle(f'/lustre/scratch127/cellgen/cellgeni/cakirb/census/model/normalized_final.pkl')
        break
    scales1 = vae.get_normalized_expression(adata[rangelist[i]:rangelist[i+1],], library_size = 10000)
    if i == 0:
        scales0 = scales1.copy()
    else:
        scales0 = pd.concat([scales0, scales1])
