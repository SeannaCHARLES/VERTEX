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

path = sys.argv[1]
studies=list(sys.argv[2].split(" "))
spe_prot=list(sys.argv[3].replace(" ","").split("\n")) 
names=list(sys.argv[4].split(" "))
path_script=sys.argv[5]
resfile=str(sys.argv[6])
#print(spe_prot)





for name,spe_prot in zip(names,spe_prot) : 
    if name in studies :
        change_name={}
        proteins=spe_prot.split("=")[1].replace('"',"").split("+") 
        
        label=pd.read_csv(f"{path}{name}/legend/Label_{name}{proteins[0]}{proteins[1]}_itol.txt",sep="\t",skiprows=4,comment="#")
        end_line=label.shape[0]
        end_col=label.shape[1]
        #change_name = label.set_index(0)[1].to_dict()
        for line in (range(0,end_line)) : 
            change_name[label.iloc[line,0]]=label.iloc[line,1]
        

        open(f"{path}{name}/temp_protdist_filter.txt", "w").close()
        
        

        cmd=f"sed 's/>/>{proteins[1][0]}/g' {path}{name}{resfile}{proteins[1]}_{name}_prot.fasta > {path}{name}{resfile}{proteins[1]}_{name}_prot_.fasta"
        subprocess.run(cmd,shell=True)

        cmd=f"sed 's/>/>{proteins[0][0]}/g' {path}{name}{resfile}{proteins[0]}_{name}_prot.fasta > {path}{name}{resfile}{proteins[0]}_{name}_prot_.fasta"
        subprocess.run(cmd,shell=True)


        for record_1 in SeqIO.parse(f"{path}{name}{resfile}{proteins[1]}_{name}_prot_.fasta","fasta") : 
            Query = record_1.id
            Qseq = record_1.seq
            for record_2 in SeqIO.parse(f"{path}{name}{resfile}{proteins[0]}_{name}_prot_.fasta","fasta") : 
                subject = record_2.id
                Sseq = record_2.seq
        
                file_name=f"{name}/temp"
                file = open(f"{path}{file_name}.fasta","w")
                file.write(f">{Query}\n{Qseq}\n>{subject}\n{Sseq}\n")       
                file.close()

                subprocess.run(f"mafft --auto {path}{file_name}.fasta > {path}{file_name}.aln-fasta",shell=True)
                subprocess.run(f"seqret -sequence {path}{file_name}.aln-fasta -outseq {path}{file_name}.phy -osformat phylip",shell=True)

                shutil.copy(f"{path}{file_name}.phy", f"infile")
                subprocess.run("phylip protdist", input="Y\n", text=True,shell=True)

                subprocess.run(f"sed 's/   /	/g' outfile > outfile_prot.txt",shell=True)
                subprocess.run(f"sed 's/  /	/g' outfile > outfile_prot.txt",shell=True)

                subprocess.run(f"rm infile outfile",shell=True)

                out=pd.read_csv(f"outfile_prot.txt",sep="\t",skiprows=0,index_col=0)
                out.columns=out.index
                try:
                  
                    final_table=pd.read_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t",index_col=0,header=0)
                    if final_table.empty:
                        
                        out.to_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t")
                    else:
                        
                        Final_table_v=out.combine_first(final_table)
                        Final_table_v.to_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t")
                except pd.errors.EmptyDataError:
                    out.to_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t")
                    
        final_table=pd.read_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t",index_col=0,header=0)
        final_table.drop(final_table.index[final_table.index.str.startswith(f'{proteins[0][0]}')],axis=0,inplace=True)
        final_table.drop(final_table.columns[final_table.columns.str.startswith(f'{proteins[1][0]}')],axis=1,inplace=True)    
        final_table.rename(columns=change_name, inplace=True,index=change_name)
        final_table.to_csv(f"{path}{name}/temp_protdist_filter.txt",sep="\t")
                
