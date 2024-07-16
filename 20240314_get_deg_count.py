# %%
import glob
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append("/BiO/Access/kyungwhan1998/pyvenn")
import venn

# %%
dir_deg = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/5_deg"
list_file_deg = glob.glob(f"{dir_deg}/**/Visit*_Severe__Visit*_Mild_20231007.tsv", recursive=True)

dict_visit_deg = dict()
for file_deg in list_file_deg:
    df_deg = pd.read_csv(file_deg, sep="\t")
    df_deg["abslog2FoldChange"] = df_deg["log2FoldChange"].apply(abs)
    df_deg_sig = df_deg[np.logical_and(df_deg["abslog2FoldChange"] > 1.3, df_deg["padj"] < 0.05)]
    visit = os.path.basename(file_deg).split("_")[0]
    list_de_fc = df_deg_sig["log2FoldChange"].to_list()
    list_de_gene = df_deg_sig["ID"].to_list()
    dict_de_gene_fc = dict(zip(list_de_gene, list_de_fc))
    dict_visit_deg[visit] = dict_de_gene_fc

print(len(dict_visit_deg["Visit1"].items()))
list_up_visit1 = [k for k, v in dict_visit_deg["Visit1"].items() if v > 0]
print(len(list_up_visit1))
list_down_visit1 = [k for k, v in dict_visit_deg["Visit1"].items() if v < 0]
print(len(list_down_visit1))

print(len(dict_visit_deg["Visit2"].items()))
list_up_visit2 = [k for k, v in dict_visit_deg["Visit2"].items() if v > 0]
print(len(list_down_visit1))
list_down_visit2 = [k for k, v in dict_visit_deg["Visit2"].items() if v < 0]
print(len(list_down_visit2))

# %%
labels = venn.get_labels([set(dict_visit_deg["Visit1"]), set(dict_visit_deg["Visit2"])])
venn.venn2(labels, names=["Visit1", "Visit2"])
plt.show()
plt.close()

# %%
labels = venn.get_labels([set(list_up_visit1), set(list_up_visit2)])
venn.venn2(labels, names=["Visit1", "Visit2"])
plt.show()
plt.close()

# %%
labels = venn.get_labels([set(list_down_visit1), set(list_down_visit2)])
venn.venn2(labels, names=["Visit1", "Visit2"])
plt.show()
plt.close()

# %%
labels = venn.get_labels([set(list_up_visit1), set(list_down_visit2)])
venn.venn2(labels, names=["Visit1", "Visit2"])
plt.show()
plt.close()

# %%
labels = venn.get_labels([set(list_down_visit1), set(list_up_visit2)])
venn.venn2(labels, names=["Visit1", "Visit2"])
plt.show()
plt.close()


# %%
# %
# %%
