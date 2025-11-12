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
spe_prot=list(sys.argv[3].replace(" ","").split("\n")) 
names=list(sys.argv[4].split(" "))
path_script=sys.argv[5]
resfile=str(sys.argv[6])
script_path=str(sys.argv[7])
#print(spe_prot)

def Lauch_MSP(list_record) : 
    global script_path
    record_1=list_record[0]
    record_2=list_record[1]
    name=list_record[2]
    Query = record_1.id
    Qseq = record_1.seq
    subject = record_2.id
    Sseq = record_2.seq

    new_path=f"{path}{name}/Align_all/{Query}{subject}/"

    
    subprocess.run(f"./{script_path}PROT_dist.sh {new_path} {Query} {Qseq} {subject} {Sseq} ",shell=True)
    # subprocess.run(f"cd {new_path}",shell=True)
    # file_name=f"{new_path}temp"
    # file = open(f"{file_name}.fasta","w")
    # file.write(f">{Query}\n{Qseq}\n>{subject}\n{Sseq}\n")       
    # file.close()

    # subprocess.run(f"mafft --auto {file_name}.fasta > {file_name}.aln-fasta",shell=True)
    # subprocess.run(f"seqret -sequence {file_name}.aln-fasta -outseq {file_name}.phy -osformat phylip",shell=True)
    # subprocess.run(f"cp {file_name}.phy {new_path}infile",shell=True)
    # #shutil.copy(f"{file_name}.phy", f"{new_path}infile")
    # subprocess.run(f"phylip protdist  {new_path}infile > {new_path}outfile", input="Y\n", text=True,shell=True)

    # subprocess.run(f"sed 's/   /	/g' {new_path}outfile > {new_path}outfile_prot.txt",shell=True)
    # subprocess.run(f"sed 's/  /	/g' {new_path}outfile > {new_path}outfile_prot.txt",shell=True)


    out=pd.read_csv(f"{new_path}outfile_prot.txt",sep="\t",skiprows=0,index_col=0)
    out.columns=out.index
    out.to_csv(f"{path}{name}/Align_all/{Query}{subject}_table.txt",sep="\t")



for name,spe_prot in zip(names,spe_prot) : 
    if name in studies :

        subprocess.run(f"rm -r {path}{name}/Align_all/*",shell=True)
        change_name={}
        couple_list=[]
        proteins=spe_prot.split("=")[1].replace('"',"").split("+") 
        
        label=pd.read_csv(f"{path}{name}/legend/Label_{name}{proteins[0]}{proteins[1]}_itol.txt",sep="\t",skiprows=4,comment="#")
        
        end_line=label.shape[0]
        end_col=label.shape[1]

        for line in (range(0,end_line)) : 
            change_name[label.iloc[line,0]]=label.iloc[line,1]
        

        open(f"{path}{name}/temp_protdist_filter.txt", "w").close()
        
        

        cmd=f"sed 's/>/>{proteins[1][0]}/g' {path}{name}{resfile}{proteins[1]}_{name}_prot.fasta > {path}{name}{resfile}{proteins[1]}_{name}_prot_.fasta"
        subprocess.run(cmd,shell=True)

        cmd=f"sed 's/>/>{proteins[0][0]}/g' {path}{name}{resfile}{proteins[0]}_{name}_prot.fasta > {path}{name}{resfile}{proteins[0]}_{name}_prot_.fasta"
        subprocess.run(cmd,shell=True)


        for record_1 in SeqIO.parse(f"{path}{name}{resfile}{proteins[1]}_{name}_prot_.fasta","fasta") : 

            for record_2 in SeqIO.parse(f"{path}{name}{resfile}{proteins[0]}_{name}_prot_.fasta","fasta") : 

                couple_list.append([record_1,record_2,name])

        if __name__ == '__main__' : 
            with Pool(10) as p : 
                print(p.map(Lauch_MSP,couple_list))


        for elements in couple_list   :
            record_1=elements[0]   
            record_2=elements[1] 

            table=pd.read_csv(f"{path}{name}/Align_all/{record_1.id}{record_2.id}_table.txt",sep="\t",skiprows=0,index_col=0,header=0)

            if os.path.getsize(f"{path}{name}/temp_protdist_filter.txt")==0 :
                table.to_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t")
                print(f"empty {os.path.getsize(f"{path}{name}/temp_protdist_filter.txt")}")
            else : 
                final_table=pd.read_csv(f"{path}{name}/temp_protdist_filter.txt")
                Final_table_v=final_table.combine_first(table)
                Final_table_v.to_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t")
                        
                


            #Final_table_v=final_table.combine_first(table)
            #Final_table_v.to_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t")

        #final_table=pd.read_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t",index_col=0,header=0)
        final_table.drop(final_table.index[final_table.index.str.startswith(f'{proteins[0][0]}')],axis=0,inplace=True)
        final_table.drop(final_table.columns[final_table.columns.str.startswith(f'{proteins[1][0]}')],axis=1,inplace=True)
        print(f"drop {os.path.getsize(f"{path}{name}/temp_protdist_filter.txt")}")    
        final_table.rename(columns=change_name, inplace=True,index=change_name)
        print(f"rename {os.path.getsize(f"{path}{name}/temp_protdist_filter.txt")}")
        final_table.to_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t")
        print(f"end {os.path.getsize(f"{path}{name}/temp_protdist_filter.txt")}")
        # for i in range(0,len(couple_list)-1,2) in couple_list   :
        #     couple1=couple_list[i]
        #     couple2=couple_list[i+1]
            
        #     record_1_C1=couple1[0]   
        #     record_2_C1=couple1[1] 

        #     record_1_C2=couple2[0]   
        #     record_2_C2=couple2[1] 

        #     table_C1=pd.readcsv(f"{path}{study}/Align_all/{record_1_C1.id}{record_2_C1.id}_table.txt",sep="\t",skiprows=0,index_col=0,header=0)
        #     table_C2=pd.readcsv(f"{path}{study}/Align_all/{record_1_C2.id}{record_2_C2.id}_table.txt",sep="\t",skiprows=0,index_col=0,header=0)

        #     Final_table_v=table_C1.combine_first(table_C2)