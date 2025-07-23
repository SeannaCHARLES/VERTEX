import sys
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
def clean_filename(s):
    return s.replace(" ", "").replace('"', "")

path = sys.argv[1]
names=list(sys.argv[2].split(" "))
spe_prot=list(sys.argv[3].replace(" ","").split("\n")) #.replace(" ", "\n").splitlines())
types_seq=list(sys.argv[4].split(";"))
res_file=str(sys.argv[5])
studies=list(sys.argv[6].split(" "))
Allstudies=list(sys.argv[7].split(" "))
proteins=list(sys.argv[8].split(" "))


# print(f" path: {path}")

# print(f" types_seq: {types_seq}")

# print(f" res_file: {res_file}")

# print(f" studies: {studies}")

# print(f" Allstudies: {Allstudies}")
# print(f" names: {names}")
# print(f" proteins: {proteins}")
# print(f" spe_prot: {spe_prot}")

dict_test={}
mafftfiles={}

onlyOG= [x for x in names if x not in studies]

for (study,name,prot,pp)in zip(Allstudies,names,spe_prot,proteins):
    prot=prot.replace('"',"")
    pp=pp.split(";")
    if prot!='' :
        compo_prot=prot.split("=")[1].split("+")
        name_spe_prot=prot.split("=")[0].replace('"',"")
        print(f"cette étude {name} a besoin de cette proteine {name_spe_prot} qui est composé de : {compo_prot}")

        for type_seq in types_seq :
            for OGprot in compo_prot : 
                need_prot_namefile=clean_filename(f"{path}{study}{res_file}{OGprot}_{name}_{type_seq}.fasta")
                OG_filename=clean_filename(f"{path}{study}{res_file}{name_spe_prot}_{study}_{type_seq}_OG.fasta")
                mafft_outfile=clean_filename(f"{path}{study}{res_file}{name_spe_prot}_{study}_{type_seq}_OG.aln-fasta")
                    
                mafftfiles[OG_filename]=mafft_outfile
                    
                file=open(OG_filename,"a")
                for record in SeqIO.parse(f"{need_prot_namefile}","fasta") : 
                    file.write(f">{OGprot[0]}{record.id}\n{record.seq}\n")
                file.close()
            if len(pp)>1 :
                for prot in pp : 
                    print(f"cette étude {study} a besoin de cette proteine {prot} de {name}")
                    need_prot_namefile=clean_filename(f"{path}{study}{res_file}{prot}_{name}_{type_seq}.fasta")
                    OG_filename=clean_filename(f"{path}{study}{res_file}{prot}_{study}_{type_seq}_OG.fasta")
                    mafft_outfile=clean_filename(f"{path}{study}{res_file}{prot}_{study}_{type_seq}_OG.aln-fasta")
                
                    mafftfiles[OG_filename]=mafft_outfile
                    
                    file=open(OG_filename,"a")
                    for record in SeqIO.parse(f"{need_prot_namefile}","fasta") : 
                        file.write(f">{record.id}\n{record.seq}\n")
                    file.close()

    else : 
        

        for type_seq in types_seq :
                for prot in pp : 
                    print(f"cette étude {study} a besoin de cette proteine {prot} de {name}")
                    need_prot_namefile=clean_filename(f"{path}{study}{res_file}{prot}_{name}_{type_seq}.fasta")
                    OG_filename=clean_filename(f"{path}{study}{res_file}{prot}_{study}_{type_seq}_OG.fasta")
                    mafft_outfile=clean_filename(f"{path}{study}{res_file}{prot}_{study}_{type_seq}_OG.aln-fasta")
                    
                    mafftfiles[OG_filename]=mafft_outfile
                    
                    file=open(OG_filename,"a")
                    for record in SeqIO.parse(f"{need_prot_namefile}","fasta") : 
                        file.write(f">{record.id}\n{record.seq}\n")
                    file.close()
# for (study,name,prot,pp)in zip(Allstudies,names,spe_prot,proteins): 
#     print("debut",name,prot)
#     if prot !="''" :            
#         comprot=list(prot.split("=")[1].split("+"))
#         for element_prot in prot.split(" ") : 
#             name_prot=element_prot.split("=")[0]

#             if name not in onlyOG :             

#                 for type_seq in types_seq :


#                     for prot in comprot :
#                         filename=clean_filename(f"{path}{study}{res_file}{name_prot}_{study}_{type_seq}_OG.fasta")
#                         namefile=clean_filename(f"{path}{study}{res_file}{prot}_{study}_{type_seq}.fasta")
#                         outfile=clean_filename(f"{path}{study}{res_file}{name_prot}_{study}_{type_seq}_OG.aln-fasta")

#                         print(f"For the specific protein{name_prot}{type_seq}")
#                         print(f"add the protein {prot}{type_seq}")
#                         for record in SeqIO.parse(f"{namefile}","fasta") :   
#                             file = open(filename,"a")
#                             file.write(f">{prot[0]}{record.id}\n{record.seq}\n")       
#                             file.close()
#                     mafftfile[filename]=[outfile]
#             else : 
#                 for type_seq in types_seq :
#                     for prot in comprot : 
#                         filename=clean_filename(f"{path}{study}{res_file}{name_prot}_{study}_{type_seq}_OG.fasta")
#                         namefile=clean_filename(f"{path}{study}{res_file}{prot}_{name}_{type_seq}.fasta")
#                         outfile=clean_filename(f"{path}{study}{res_file}{name_prot}_{study}_{type_seq}_OG.aln-fasta")

#                         print(f"For the specific protein{name_prot}{type_seq}")
#                         print(f"add the protein {prot}{type_seq}")
#                         for record in SeqIO.parse(f"{namefile}","fasta") :
#                             file = open(filename,"a")
#                             file.write(f">{prot[0]}{record.id}\n{record.seq}\n")       
#                             file.close()
#                     mafftfile[filename]=[outfile]
        
            
#     else : 
#         for element in pp.split(";") : 
#             filename=clean_filename(f"{path}{study}{res_file}{element}_{study}_{type_seq}_OG.fasta")
#             namefile=clean_filename(f"{path}{study}{res_file}{element}_{name}_{type_seq}.fasta")
#             outfile=clean_filename(f"{path}{study}{res_file}{name_prot}_{study}_{type_seq}_OG.aln-fasta")

#             print(f"For the specific protein{name_prot}{type_seq}")
#             print(f"add the protein {prot}{type_seq}")
#             for record in SeqIO.parse(f"{namefile}","fasta") :
#                 file = open(filename,"a")
#                 file.write(f">{prot[0]}{record.id}\n{record.seq}\n")       
#                 file.close()    
#             mafftfile[filename]=[outfile]
     

for filename,outfile in mafftfiles.items() : 
    subprocess.run(f"mafft --auto {filename} > {outfile}",shell=True)
