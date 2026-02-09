import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np
import dalmatian
import anndata as ad


wm = dalmatian.WorkspaceManager('getzlab-peddep-hcr-seq/Nagel-scHCRseq-prism-pilot')
S = wm.get_samples()

S[['pool','timepoint']] = S.reset_index()['sample_id'].str.split('_',expand=True).values
S = S[S['pool']=='P3']


adatas= list()
for ind,row in S.iterrows():
    adata = read_h5ad_gs(row['repair_annotated_h5ad'])
    df = pd.read_csv(row['demuxlet_output'] + '/' + ind +'.best',sep='\t')
    adata.obs = adata.obs.join(df.set_index('BARCODE'),how='left')
    adata.obs['timepoint'] = row['timepoint']
    adata.obs['pool'] = row['pool']
    adata.var_names_make_unique()
    adatas.append(adata)

adata = ad.concat(adatas)
adata.obs_names_make_unique()

sc.pp.filter_cells(adata, min_genes=200)  
sc.pp.filter_genes(adata, min_cells=3)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
adata.var_names_make_unique()


## FILTERING
adata = adata[adata.obs['pct_counts_mt']<7.5]
adata = adata[adata.obs['total_counts']>5000]
adata = adata[adata.obs['DROPLET.TYPE']=='SNG']
adata.obs['cell_line'] = adata.obs['SNG.BEST.GUESS']


adata.obs['transfected'] = adata.obs['GFP_BAR_010']>=50

## Perform scRNA processing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    flavor="seurat_v3",
)

adata.raw = adata

adata.layers["scaled"] = adata.X.toarray()
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"], layer="scaled")
sc.pp.scale(adata, max_value=10, layer="scaled")

sc.pp.pca(adata, layer="scaled", svd_solver="arpack")

sc.external.pp.harmony_integrate(adata,'timepoint')
sc.pp.neighbors(adata, use_rep = 'X_pca_harmony')
sc.tl.umap(adata)
sc.tl.leiden(adata)

adata.write("data/P3_clustered.v0.h5ad")

