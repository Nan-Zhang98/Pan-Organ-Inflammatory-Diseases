#### Module ------------------------------------------
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import subprocess
import os
import re
import scvi
import seaborn as sns
from scib_metrics.benchmark import Benchmarker

#### load data ----------------------------------------
moMac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/moMac.h5ad")

moMac10 = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/sctour/data/Mac10_sctour.h5ad")

#### PAGA analysis -------------------------------------

sc.tl.paga(moMac, groups='cell_type_level3')
sc.pl.paga(moMac)

moMac.uns["iroot"] = np.flatnonzero(moMac.obs["cell_type_level3"] == "IL1B+NLRP3+_Macro")[0]
sc.tl.dpt(moMac)


sc.tl.draw_graph(moMac, init_pos="paga")
sc.pl.draw_graph(moMac, color=["cell_type_level3", "dpt_pseudotime"], legend_loc="on data")

#### PAGA analysis (downsampled 10%)-------------------------------------

sc.tl.paga(moMac10, groups='cell_type_level3')
sc.pl.paga(moMac10)
sc.pl.paga(moMac10, show=False)
plt.savefig('/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/PAGA/result/net_moMac10_PAGA.pdf')


moMac.uns["iroot"] = np.flatnonzero(moMac.obs["cell_type_level3"] == "IL1B+NLRP3+_Macro")[0]
sc.tl.dpt(moMac10)


sc.tl.draw_graph(moMac, init_pos="paga")
sc.pl.draw_graph(moMac, color=["cell_type_level3", "dpt_pseudotime"], legend_loc="on data")

