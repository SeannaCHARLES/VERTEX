
import sys

from All_Functions import select_doublon,delete_dbl,col_name

path ="/home/seannacharles/Documents/"


study=str(sys.argv[1])
target=str(sys.argv[2])
ths=float(sys.argv[3])
Prot_list=list(sys.argv[4])
print(Prot_list) 

final_dbl=[]
old_col_names=[]

for protein in Prot_list : 

    Doublon,name_colonnes=select_doublon(protein,study,target,ths)
    old_col_names=col_name(old_col_names,name_colonnes,protein)

    if len(final_dbl)> len(Doublon) or len(final_dbl)==0 : 
        final_dbl=Doublon

for protein in Prot_list : 
        delete_dbl(final_dbl,study,protein,old_col_names)

doublon_file=f"{path}{study}/Doublon_file_{study}.txt"
dbl_file=open(doublon_file,'a')

for db in final_dbl : 
    dbl_file.write(f"{db}\n")
dbl_file.close()


NG_old_col_names=[]
NG_final_dbl=[]

for name in old_col_names:
    NG_old_col_names.append(str("N"+name))
    NG_old_col_names.append(str("G"+name))

NG_old_col_names.append(str("Gkhurdunyp"))#Khurdun is not present in the sequence file for the other proteins because it will only be used for NSmGn as an outgroup.

for name in final_dbl:
    NG_final_dbl.append(str("N"+name))
    NG_final_dbl.append(str("G"+name))


delete_dbl(NG_final_dbl,study,"NSmGn",NG_old_col_names)

