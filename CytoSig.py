import os
import sys
import pandas as pd
import numpy as np
import subprocess
import re
import scanpy as sc

import CytoSig

pan_Macro_filtered = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Macrophage/pan_Macro_filtered_scanvi.h5ad")

Y = pd.DataFrame(pan_Macro_filtered.X,
                 index=pan_Macro_filtered.obs_names,
                 columns=pan_Macro_filtered.var_names)

Y = pan_Macro_filtered.to_df().T

signature = CytoSig.find_signature_path()
signature = pd.read_csv(signature, sep='\t', index_col=0) #43 个高置信度细胞因子

beta, std, zscore, pvalue = CytoSig.ridge_significance_test(signature, Y, alpha=1E4, alternative="two-sided", nrand=1000, cnt_thres=10, flag_normalize=True, verbose = True)

beta.to_csv('/storage/data/KAI/Pan_Inflammation/Macrophage/Macrophage_CytoSig_beta.txt', sep = '\t', index = True)
std.to_csv('/storage/data/KAI/Pan_Inflammation/Macrophage/Macrophage_CytoSig_std.txt', sep = '\t', index = True)
zscore.to_csv('/storage/data/KAI/Pan_Inflammation/Macrophage/Macrophage_CytoSig_zscore.txt', sep = '\t', index = True)
pvalue.to_csv('/storage/data/KAI/Pan_Inflammation/Macrophage/Macrophage_CytoSig_pvalue.txt', sep = '\t', index = True)



