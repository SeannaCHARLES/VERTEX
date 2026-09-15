import pandas as pd

Diff_name={}

df=pd.read_csv(f"/home/seannacharles/Protein_equivalent.txt",comment="#",sep="\t",header=None)

end_line=df.shape[0]

for line in range(0,end_line-1,1) :
    Diff_name[df.iloc[line,0]]=df.iloc[line,1]
    #print(df.iloc[line,0],df.iloc[line,1])

print(Diff_name)
for prot,name in Diff_name.items() :
    print(prot,name)