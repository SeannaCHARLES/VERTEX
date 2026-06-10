##########################Description##########################
""" This script delete all the sequences that are to close from each other.
According to the threshold choosen in information_file.txt.
This script need the DNA_dist matrix as input.

"""
##########################################################

##########################Import##########################
import sys
import subprocess
import pandas as pd

from Bio import SeqIO

from All_Functions import select_doublon,delete_dbl,Modif_ID
##########################################################

###########################Input##########################
path =str(sys.argv[1])
studies=list(sys.argv[2].split(" "))
keep_list=list(sys.argv[3].replace("\n","").split(";"))
ths_list=list(sys.argv[4].split())
target_list=list(sys.argv[5].split())

keep=[]
##########################################################


##########################Script##########################



for element in keep_list : 
    keep.append(Modif_ID(element))

for study,ths,target in zip(studies,ths_list,target_list) : 
    Doublon=[]

    Doublon=select_doublon(path,study,float(ths),keep,Modif_ID(target))
    delete_dbl(path,Doublon,study)

    doublon_file=f"{path}{study}/Doublon_file_{study}.txt"

    subprocess.run(f"rm {doublon_file}",shell=True)
    for db in Doublon : 
        dbl_file=open(doublon_file,'a')
        dbl_file.write(f"{db}\n")
        dbl_file.close()
        
    change_name={}
    
    label=pd.read_csv(f"{path}{study}/legend/Label_{study}_itol.txt",sep="\t",skiprows=4,comment="#")
        
    end_line=label.shape[0]
    end_col=label.shape[1]

    for line in (range(0,end_line)) : 
        change_name[label.iloc[line,0]]=label.iloc[line,1]

    out=pd.read_csv(f"{path}{study}/{study}_dnadist.txt",sep="\t",skiprows=0,index_col=0)
  
    out.rename(columns=change_name, inplace=True,index=change_name)
  
    out.to_csv(f"{path}{study}/{study}_dnadist_filter.txt",sep="\t")

##########################################################

