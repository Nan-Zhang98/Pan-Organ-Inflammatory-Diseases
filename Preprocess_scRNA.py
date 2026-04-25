# Module
import scanpy as sc
import pandas as pd
import subprocess
import os
# pandas显示所有列
#pd.set_option('display.max_columns', None)

def adataQC(adata):
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
    '''
    sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
        )
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    '''
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs['pct_counts_mt'] < 20, :]

    sc.pp.scrublet(adata, batch_key="Sample_geo_accession")
    adata = adata[~adata.obs['predicted_doublet'],:]

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    #sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="Sample_geo_accession")
    #sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    #sc.pl.highly_variable_genes(adata)
    #adata = adata[:,adata.var['highly_variable']]

    sc.tl.pca(adata)
    #sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    '''
    sc.pl.pca(
        adata,
        color=["Sample_geo_accession", "Sample_geo_accession", "pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
    )
    '''
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    '''
        sc.pl.umap(
            adata,
            color="Sample_geo_accession",
            # Setting a smaller point size to get prevent overlap
            size=2,
        )
    '''
    return adata

def main():
    gse_id = "3.GSE149878"
    raw_folder = r"/storage/data/KAI/Pan_Inflammation/Lung/COVID-19/"
    new_folder = os.path.join(raw_folder, gse_id) + "/" + gse_id.split(".")[1] + "_py/"
    h5ad_file = os.listdir(new_folder)

    adata = sc.read_h5ad(os.path.join(new_folder, [file for file in h5ad_file if file.endswith('_scanpy.h5ad')][0]))
    adata.raw = adata
    adata = adataQC(adata)
    adata.write(new_folder + gse_id.split(".")[1] + "_preprocess.h5ad")
if __name__ == "__main__":
    main()

adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Joint/Inflammatory arthritis/1.E-MTAB-11791/E-MTAB-11791_py/E-MTAB-11791_preprocess.h5ad")


adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Lung/Chronic obstructive pulmonary disease/6.GSE270667/GSE270667_py/GSE270667_scanpy.h5ad")
adata.obs['Sample_Site'].unique()
lung_organoid_cells = adata.obs['Sample_Site'] == 'Lung_organoid'
adata = adata[~lung_organoid_cells, :]
adata.write("/storage/data/KAI/Pan_Inflammation/Lung/Chronic obstructive pulmonary disease/6.GSE270667/GSE270667_py/GSE270667_scanpy.h5ad")


