import sys
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq

def clean_filename(s):
    return s.replace(" ", "").replace('"', "")

path = sys.argv[1]
names=list(sys.argv[2].split(" "))
spe_prot=list(sys.argv[3].replace(" ","").split("\n")) 
types_seq=list(sys.argv[4].split(";"))
res_file=str(sys.argv[5])
studies=list(sys.argv[6].split(" "))
Allstudies=list(sys.argv[7].split(" "))
proteins=list(sys.argv[8].split(" "))





dict_test={}
mafftfiles={}

onlyOG= [x for x in names if x not in studies]

for (study,name,prot,pp)in zip(Allstudies,names,spe_prot,proteins):
    prot=prot.replace('"',"")
    pp=pp.split(";")
    
    if prot!='':
        compo_prot=prot.split("=")[1].split("+")
        name_spe_prot=prot.split("=")[0].replace('"',"")


        for type_seq in types_seq :
            for OGprot in compo_prot : 
                need_prot_namefile=clean_filename(f"{path}{study}{res_file}{OGprot.lower()}_{name}_{type_seq}.fasta")
                OG_filename=clean_filename(f"{path}{study}{res_file}{name_spe_prot}_{study}_{type_seq}_OG.fasta")
                mafft_outfile=clean_filename(f"{path}{study}{res_file}{name_spe_prot}_{study}_{type_seq}_OG.aln-fasta")
                print(need_prot_namefile)
                mafftfiles[OG_filename]=mafft_outfile
                    
                file=open(OG_filename,"a")
                for record in SeqIO.parse(f"{need_prot_namefile}","fasta") : 
                    file.write(f">{OGprot.lower()}{record.id}\n{record.seq}\n")
                file.close()
            if len(pp)>1 :
                for prot in pp : 

                    need_prot_namefile=clean_filename(f"{path}{study}{res_file}{prot.lower()}_{name}_{type_seq}.fasta")
                    OG_filename=clean_filename(f"{path}{study}{res_file}{prot.lower()}_{study}_{type_seq}_OG.fasta")
                    mafft_outfile=clean_filename(f"{path}{study}{res_file}{prot.lower()}_{study}_{type_seq}_OG.aln-fasta")
                
                    mafftfiles[OG_filename]=mafft_outfile
                    
                    file=open(OG_filename,"a")
                    for record in SeqIO.parse(f"{need_prot_namefile}","fasta") : 
                        file.write(f">{record.id}\n{record.seq}\n")
                    file.close()

    else : 
        

        for type_seq in types_seq :
                for prot in pp : 

                    need_prot_namefile=clean_filename(f"{path}{study}{res_file}{prot.lower()}_{name}_{type_seq}.fasta")
                    OG_filename=clean_filename(f"{path}{study}{res_file}{prot.lower()}_{study}_{type_seq}_OG.fasta")
                    mafft_outfile=clean_filename(f"{path}{study}{res_file}{prot.lower()}_{study}_{type_seq}_OG.aln-fasta")
                    
                    mafftfiles[OG_filename]=mafft_outfile
                    
                    file=open(OG_filename,"a")
                    for record in SeqIO.parse(f"{need_prot_namefile}","fasta") : 
                        file.write(f">{record.id}\n{record.seq}\n")
                    file.close()


for filename,outfile in mafftfiles.items() : 
    subprocess.run(f"mafft --auto {filename} > {outfile}",shell=True)
