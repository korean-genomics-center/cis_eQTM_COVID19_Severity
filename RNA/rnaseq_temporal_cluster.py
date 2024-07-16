#%%
# from certifi import where
# from importlib_metadata import files
import numpy as np
import pandas
from matplotlib import pyplot as plt

files_to_include=["/BiO/Research/Project1/KogicInfectomics/Results/RNASeq_DEG/RepeatMeasure/Visit1_Visit2_DEG.tsv","/BiO/Research/Project1/KogicInfectomics/Results/RNASeq_DEG/RepeatMeasure/Visit1_Visit3_DEG.tsv","/BiO/Research/Project1/KogicInfectomics/Results/RNASeq_DEG/RepeatMeasure/Visit1_Visit4_DEG.tsv"];xlabels=['V2/V1','V3/V1','V4/V1']; directory_out="/BiO/Research/Project1/KogicInfectomics/Results/RNASeq_DEG/TemporalCluster/V2V1_V3V1_V4V1"
#files_to_include=["/BiO/Research/Project1/KogicInfectomics/Results/RNASeq_DEG/RepeatMeasure/Visit1_Visit2_DEG.tsv","/BiO/Research/Project1/KogicInfectomics/Results/RNASeq_DEG/RepeatMeasure/Visit2_Visit3_DEG.tsv","/BiO/Research/Project1/KogicInfectomics/Results/RNASeq_DEG/RepeatMeasure/Visit3_Visit4_DEG.tsv"];xlabels=['V2/V1','V3/V2','V4/V3'];
column_with_values="log2FoldChange"
genes_column="ID"
fin_genes_to_include="/BiO/Research/Project1/KogicInfectomics/Results/RNASeq_DEG/RepeatMeasure/V1234_DEG.all.tsv"
fin_ensembl_to_genesymbol_table="/BiO/Share/Resource/References/geneInfo.tab"

#%%
import os

if(not os.path.exists(directory_out)):
    os.makedirs(directory_out)

def get_ensembl_to_genesymbol_dict(fin=fin_ensembl_to_genesymbol_table):
    en_sym_dict={}
    with open(fin) as fh:
        fh.readline()
        aline=fh.readline()
        while(aline):
            (en, gs)=aline.rstrip().split("\t")
            en_sym_dict[en]=gs
    return en_sym_dict


#%%
genes=[]
with open(fin_genes_to_include) as fh:
    aline=fh.readline()
    while(aline):
        genes.append(aline.strip())
        aline=fh.readline()

#%%
import numpy as np

logfoldlist=np.ones((len(genes),len(files_to_include)))

for idx,afile in enumerate(files_to_include):
    with open(afile) as fh:
        aline=fh.readline()
        vs=aline.rstrip().split("\t")
        targetcol=vs.index(column_with_values)
        namescol=0
        aline=fh.readline()

        while(aline):
            parsed=aline.strip().split("\t")
            val=float(parsed[targetcol])
            genename=parsed[namescol]
            if(genename in genes):
                gid=genes.index(genename)
                logfoldlist[gid,idx]=val
            aline=fh.readline()


#%%
fileout=os.path.join(directory_out, "cluster.all.list.txt")
with open(fileout,'w') as fh:
    fh.write("ID\t"+"\t".join(xlabels)+"\n")
    for anid in range(len(logfoldlist)):
        fh.write(genes[anid]+"\t"+"\t".join(map(str,logfoldlist[anid,:]))+"\n")



#%%
from matplotlib import pylab as plt

ncluster=8
from scipy.cluster.vq import kmeans, kmeans2, vq, whiten

centroids, mean_value = kmeans2(logfoldlist, ncluster)
clusters,distances=vq(logfoldlist,centroids)

for idx in range(ncluster):
    curcluster=clusters==idx
    curids = np.where(curcluster)[0]
    fileout=os.path.join(directory_out, f"cluster.{idx}.list.txt")
    with open(fileout,'w') as fh:
        fh.write("ID\t"+"\t".join(xlabels)+"\n")
        for anid in curids:
            fh.write(genes[anid]+"\t"+"\t".join(map(str,logfoldlist[anid,:]))+"\n")

#%%

import matplotlib

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15 }

matplotlib.rc('font', **font)
plt.figure(figsize=(4,12))

cs=['tab:blue','tab:orange','tab:red','tab:purple','tab:pink','tab:gray','tab:olive','tab:cyan']
for idx in range(ncluster):
    curcluster=clusters==idx
    curvals=logfoldlist[curcluster,:]
#    plt.subplot(ncluster, 1, idx+1)
    plt.plot(curvals.T,'-',c=cs[idx],alpha=0.2)
    plt.plot(curvals.T,'o',c=cs[idx],alpha=0.7)
plt.xticks([0,1,2],labels=xlabels)
plt.ylim((-10,15))
plt.ylabel(r'$log_2$(Fold Change)')
plt.savefig(os.path.join(directory_out,'allcluster.pdf'))

for idx in range(ncluster):
    plt.figure(figsize=(4,12))
    curcluster=clusters==idx
    curvals=logfoldlist[curcluster,:]
#    plt.subplot(ncluster, 1, idx+1)
    plt.plot(curvals.T,'-',c=cs[idx],alpha=0.2)
    plt.plot(curvals.T,'o',c=cs[idx],alpha=0.7)
    plt.xticks([0,1,2],labels=xlabels)
    plt.ylim((-10,15))
    plt.ylabel(r'$log_2$(Fold Change)')
    plt.savefig(os.path.join(directory_out,f'cluster_{idx}.pdf'))

       
# %%
# %%