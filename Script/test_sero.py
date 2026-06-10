import numpy as np

import pandas as pd

def find_serogroup(organism,path) : 
    df = pd.read_csv(f"{path}Serogroup.txt", delimiter="|",header=0,index_col=0)
    df.iloc[:, 1:]=df.iloc[:, 1:].astype(str)

    if organism in df.index: 
        return(df.loc[organism].iloc[0])
    else: 
        name=f"{path}Serogroup.txt"

        file = open(name,"a")
        file.write(f"{organism}|nan\n")       
        file.close()
        return "nan"

sero=find_serogroup("virus1","/home/seannacharles/")
print(sero)