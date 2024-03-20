#assign logical group identification to maize_ortholog_degTable.tsv

#import necessary packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#read in ortholog DEG table
df = pd.read_csv("maize_ortholog_degTable.tsv", sep='\t', header=0)


#for each genotype, specify each of 3 possible DEG outcomes (up, down, notsignificant) and output a boolean filter for df

b73up = df['log2FoldChange_B73'] > 0

b73down = df['log2FoldChange_B73'] < 0

b73ns = df['log2FoldChange_B73'] == 0

ky21up = df['log2FoldChange_Ky21'] > 0

ky21down = df['log2FoldChange_Ky21'] < 0

ky21ns = df['log2FoldChange_Ky21'] == 0

oh7bup = df['log2FoldChange_Oh7b'] > 0

oh7bdown = df['log2FoldChange_Oh7b'] < 0

oh7bns = df['log2FoldChange_Oh7b'] == 0


#construct a list of each genotype's possible outcomes
b=[b73up, b73down, b73ns]
bn=["b73up", "b73down", "b73ns"]

k=[ky21up, ky21down, ky21ns]
kn=["ky21up", "ky21down", "ky21ns"]

o=[oh7bup, oh7bdown, oh7bns]
on=["oh7bup", "oh7bdown", "oh7bns"]

#initialize an empty list named 'categories' and iterate through the genotype DEG outcome lists to produce a list of all possible categories

categories=[]

for i in bn:
        for x in kn:
                for y in on:
                    categories.append(i + ", " + x + ", " + y)

#initialize a place-holder column in df named 'logicCat' where all records are initially assigned values of 'none'


df['logicCat'] = 'none'


#iterate through the lists of categories and all possible outcomes, using .loc[] to assign the appropriate value of 'categories' to the 'logicCat' column.

a=0
while a<len(categories):
    for i in b:
        for x in k:
            for y in o:
                df.loc[(i) & (x) & (y), 'logicCat'] = categories[a]
                a += 1


#sort df by logicCat, log2FoldChange_B73, log2FoldChange_Ky21, log2FoldChange_Oh7b, and plot a heatmap of the log2foldchange data

dfs = df.sort_values(by=['logicCat', 'log2FoldChange_B73', 'log2FoldChange_Ky21', 'log2FoldChange_Oh7b'], ascending=[False, False, False, False])

dfs.set_index('B73SyntelogId', inplace=True)

dfs.to_csv('logicClusteredSyntelogDEGs.tsv', sep='\t')

#get index number of first row in each group

dividers=list(dfs.reset_index(drop=True).groupby('logicCat').head(1).index)


clusterdf = dfs[['log2FoldChange_B73', 'log2FoldChange_Ky21', 'log2FoldChange_Oh7b', 'logicCat']]

clusterdf.set_index('logicCat', inplace=True)


plt.figure(figsize=(8,6))
col=sns.diverging_palette(0, 180, s=100, as_cmap=True)
ax=sns.heatmap(clusterdf, cmap=col, center=0, yticklabels=False, robust=True)
ax.set_yticks(dividers)
ax.set_yticklabels(usedCat)
ax.set_title('Heatmap of log2(flood/control) for Logically Grouped Flood DEGs')
ax.set_ylabel('')
plt.savefig(fname='./logic_clustered_maize_DEG_heatmap.pdf', format='pdf')




#split dfs by logicCat to genelists for each used category
#build lists of categories and output filenames 
categories=list(dfs['logicCat'].unique())
names=[]
for a in categories:
    names.append(a.replace(", ", "_") + ".tsv")

splits=[]
for z in categories:
        splits.append(dfs[dfs['logicCat'] == z])

for num, frame in enumerate(splits):
   frame['B73SyntelogId'].to_csv("./logicClusterLists/" + names[num], sep='\t', header=False, index=False)


###Annotate known DEG syntelog TFs (grassius) in logicClusteredGeneSet

TFs=pd.read_csv("/psi/basslab/zturpin/genomes/Zea_mays/B73_v5/B73v5_GrassiusTF_Mar15_2024_zmtCommaFix.csv", header=0, sep=',', usecols=['protein name', 'family', 'gene ID'])

##Extract records from logicClusteredGeneSet that correspond to a maize Grassius TF

lCGS = pd.read_csv('logicClusteredSyntelogDEGs.tsv', sep='\t')

lCTF = lCGS.merge(TFs, how='inner', left_on='B73SyntelogId', right_on='gene ID')

lCTF.to_csv("logicClusteredSyntelogDEG_TFs_Grassius_mar152024.tsv", sep='\t')

