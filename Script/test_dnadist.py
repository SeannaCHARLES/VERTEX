import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import AlignIO
from Bio import Phylo


df = pd.read_csv(f"/home/seannacharles/Retrovirus_EL_Rtarticle/Retrovirus_EL_Rtarticle_dnadist.txt", delimiter="\t",index_col=0,header=0)
df.iloc[:, 1:]=df.iloc[:, 1:].astype(float)
# df = df.set_index('Unnamed: 0')
rename=df.index.tolist()
rename.append("empty")
df.columns =rename
df=df.drop(columns="empty")
end_line=df.shape[0]
end_col=df.shape[1]


for line in range(0,end_line) : 
    for col in range (0,end_col) : 
        #print(f"name line : {df.index[line]} name col : {df.columns[col]} Dist : {df.iloc[line,col]}")
        if df.iloc[line,col] <= 0.08 and df.iloc[line,col]>0 : 
            print(f"need to delete : name line : {df.index[line]} name col : {df.columns[col]} Dist : {df.iloc[line,col]}")