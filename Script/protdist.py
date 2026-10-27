import subprocess
import sys
import os
import shutil
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import AlignIO
from Bio import Phylo
from multiprocessing import Pool



path = sys.argv[1]
studies=list(sys.argv[2].split(" "))
spe_prot_list=list(sys.argv[3].replace(" ","").split("\n")) 
names=list(sys.argv[4].split(" "))
path_script=sys.argv[5]
resfile=str(sys.argv[6])
script_path=str(sys.argv[7])
list_prot = [x.split(';') for x in sys.argv[8].split()]
studiesspeprot=list(sys.argv[9].split(" "))


def Lauch_MSP(list_record) : 

    global script_path
    record_1=list_record[0]
    record_2=list_record[1]
    name=list_record[2]
    name_prot=list_record[3]
    Query = record_1.id
    Qseq = record_1.seq
    subject = record_2.id
    Sseq = record_2.seq

    new_path=f"{path}{name}/"

    subprocess.run(f"./{script_path}PROT_dist.sh {new_path} {Query} {Qseq} {subject} {Sseq} {name_prot}",shell=True)


if studiesspeprot : 
    for name,spe_prot in zip(names,spe_prot_list) : 
    
        if name in studies and spe_prot!='""' and spe_prot!='' and "+" in  spe_prot:
   
            subprocess.run(f"rm -fr {path}{name}/Align_all/*",shell=True)
            subprocess.run(f"rm -fr {path}{name}/All_table/*",shell=True)
            change_name={}
            couple_list=[]
            
            proteins=spe_prot.split("=")[1].replace('"',"").split("+") 
    
        
            label=pd.read_csv(f"{path}{name}/legend/label_{name}{proteins[0].lower()}{proteins[1].lower()}_itol.txt",sep="\t",skiprows=4,comment="#")
       
            end_line=label.shape[0]
            end_col=label.shape[1]

            for line in (range(0,end_line)) : 
                change_name[label.iloc[line,0][1:10]]=label.iloc[line,1][1:10]
            name_prot=str(proteins[0].lower())+str(proteins[1].lower())
        
            subprocess.run(f"rm -f {path}{name}/temp_protdist_filter{name_prot}.txt",shell=True)
            cmd=f"sed 's/>/>{proteins[1].upper()[0]}/g' {path}{name}{resfile}{proteins[1].lower()}_{name}_prot.fasta > {path}{name}{resfile}{proteins[1].lower()}_{name}_prot_.fasta"
            subprocess.run(cmd,shell=True)

            cmd=f"sed 's/>/>{proteins[0].upper()[0]}/g' {path}{name}{resfile}{proteins[0].lower()}_{name}_prot.fasta > {path}{name}{resfile}{proteins[0].lower()}_{name}_prot_.fasta"
            subprocess.run(cmd,shell=True)


            name_prot=str(proteins[0].lower())+str(proteins[1].lower())
            if os.path.exists(f"{path}{name}{resfile}{proteins[1].lower()}_{name}_prot.fasta") and os.path.exists(f"{path}{name}{resfile}{proteins[0].lower()}_{name}_prot.fasta") :
                for record_1 in SeqIO.parse(f"{path}{name}{resfile}{proteins[1].lower()}_{name}_prot_.fasta","fasta") : 

                    for record_2 in SeqIO.parse(f"{path}{name}{resfile}{proteins[0].lower()}_{name}_prot_.fasta","fasta") : 
                        if record_1.seq!="	" and record_2.seq!="	":
                            
                            Lauch_MSP([record_1,record_2,name,name_prot])
                            try :
                                table=pd.read_csv(f"{path}{name}/All_table/{name_prot}{record_1.id}{record_2.id}_table.txt",sep="\t",skiprows=0,index_col=0,header=0)
            
                                if not os.path.exists(f"{path}{name}/temp_protdist_filter{name_prot}.txt"):
             
                                    table.to_csv(f"{path}{name}/temp_protdist_filter{name_prot}.txt",sep="\t")
                                else:
                                    final_table=pd.read_csv(f"{path}{name}/temp_protdist_filter{name_prot}.txt",sep="\t",index_col=0,header=0)
                                    Final_table_v=final_table.combine_first(table)
                                    Final_table_v.to_csv(f"{path}{name}/temp_protdist_filter{name_prot}.txt",sep="\t")

                            except Exception as e :
                                Lauch_MSP([record_1,record_2,name,name_prot])
                                table=pd.read_csv(f"{path}{name}/All_table/{name_prot}{record_1.id}{record_2.id}_table.txt",sep="\t",skiprows=0,index_col=0,header=0)
                                final_table=pd.read_csv(f"{path}{name}/temp_protdist_filter{name_prot}.txt",sep="\t",index_col=0,header=0)
                                Final_table_v=final_table.combine_first(table)
                                Final_table_v.to_csv(f"{path}{name}/temp_protdist_filter{name_prot}.txt",sep="\t")

                if os.path.exists(f"{path}{name}/temp_protdist_filter{name_prot}.txt") :
                    print(f"Creating the protdist table for the protein : {name_prot} study :{name}")           
                    final_table=pd.read_csv(f"{path}{name}/temp_protdist_filter{name_prot}.txt",sep="\t",index_col=0,header=0)    
                    final_table.drop(final_table.index[final_table.index.str.startswith(f'{proteins[0].upper()[0]}')],axis=0,inplace=True)
                    final_table.drop(final_table.columns[final_table.columns.str.startswith(f'{proteins[1].upper()[0]}')],axis=1,inplace=True) 
                    final_table.rename(columns=change_name, inplace=True,index=change_name)   
                    subprocess.run(f"rm -f {path}{name}/{name}_protdist_filter{name_prot}.txt",shell=True)   
                    final_table.to_csv(f"{path}{name}/{name}_protdist_filter{name_prot}.txt",sep="\t")
                    subprocess.run(f"rm -f {path}{name}/temp_protdist_filter{name_prot}.txt",shell=True)

                            


for name,prots in zip(names,list_prot) : 
    
    subprocess.run(f"rm -f {path}{name}/Maffts_protdist.log",shell=True)
    if name in studies :
    
        subprocess.run(f"rm -fr {path}{name}/Align_all/*",shell=True)
        subprocess.run(f"rm -fr {path}{name}/All_table/*",shell=True)
        change_name={}
        couple_list=[]
        
        label=pd.read_csv(f"{path}{name}/legend/label_{name}_itol.txt",sep="\t",skiprows=4,comment="#")
        end_line=label.shape[0]
        end_col=label.shape[1]

        for line in (range(0,end_line)) : 
            change_name[label.iloc[line,0][1:10]]=label.iloc[line,1][1:10]


        for prot in prots:
  
            couple_list=[]
            open(f"{path}{name}/temp_protdist_filter{prot}.txt", "w").close()
            if os.path.exists(f"{path}{name}{resfile}{prot}_{name}_prot.fasta") :
     
                for record_1 in SeqIO.parse(f"{path}{name}{resfile}{prot}_{name}_prot.fasta","fasta") : 
                    for record_2 in SeqIO.parse(f"{path}{name}{resfile}{prot}_{name}_prot.fasta","fasta") : 
                        if record_1.id != record_2.id and (record_1.seq!="	" and record_2.seq!="	"):
                       
                            Lauch_MSP([record_1,record_2,name,prot])
                            try :
                                table=pd.read_csv(f"{path}{name}/All_table/{prot}{record_1.id}{record_2.id}_table.txt",sep="\t",skiprows=0,index_col=0,header=0)
                    
                                if not os.path.exists(f"{path}{name}/temptable_protdist_filter{prot}.txt"):
             
                                    table.to_csv(f"{path}{name}/temptable_protdist_filter{prot}.txt",sep="\t")
                                else:
                                    final_table=pd.read_csv(f"{path}{name}/temptable_protdist_filter{prot}.txt",sep="\t",index_col=0,header=0)
                        
                                    Final_table_v=final_table.combine_first(table)
                        
                                    Final_table_v.to_csv(f"{path}{name}/temptable_protdist_filter{prot}.txt",sep="\t")

                            except Exception as e :
                                print(f"{e}")
            if os.path.exists(f"{path}{name}/temptable_protdist_filter{prot}.txt") :
                print(f"Creating the protdist table for the protein : {prot} study :{name}")
                final_table=pd.read_csv(f"{path}{name}/temptable_protdist_filter{prot}.txt",sep="\t",index_col=0,header=0)    
                final_table.rename(columns=change_name, inplace=True,index=change_name)      
                subprocess.run(f"rm -f {path}{name}/{name}_protdist_filter{prot}.txt",shell=True)  
                final_table.to_csv(f"{path}{name}/{name}_protdist_filter{prot}.txt",sep="\t")
                subprocess.run(f"rm -f {path}{name}/temp_protdist_filter{prot}.txt",shell=True)
                          