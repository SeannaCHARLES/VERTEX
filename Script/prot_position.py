import subprocess
import sys 
import pandas as pd

path = sys.argv[1]
studies=list(sys.argv[2].split(" "))
studiesspeOG=list(sys.argv[3].split(" "))
orgaWIWO=list(sys.argv[4].replace(" ", "\n").splitlines())
INOUT=list(sys.argv[5].replace(" ", "\n").splitlines())
protINOUT=list(sys.argv[6].replace(" ", "\n").splitlines())
target_list=list(sys.argv[7].split())
names=list(sys.argv[8].split(" "))
Allstudies=list(sys.argv[9].split(" "))
OG_protOUT={}
OG_protIN={}

def Blastp_run(DB,File,output) :
    subprocess.run(f"rm -f {output}.txt",shell=True)
    subprocess.run(f"blastp -db {DB}_DB -query {File} -outfmt 6 -out {output}.txt",shell=True)
    # cmd=str("bioawk -c fastx '{ print $name, length($seq) }' ")
    # subprocess.run(f"{cmd}{DB}.fasta",shell=True)
    #print(f"{File}"
    # subprocess.run(f"water -asequence {DB}.fasta -bsequence {File} -outfile {output}_water.txt -auto",shell=True)
    # subprocess.run(f"needle -asequence {DB}.fasta -bsequence {File} -outfile {output}_needle.txt -auto",shell=True)
    
def add_new_pos(protein,output) :
    prot=protein.replace("Start","").replace("End","")
    if "Start" in protein :
        missing_pos=7
    else : 
        missing_pos=8

    #cmd=str("awk '/"+{prot}+"/ && (NR==1 || $3+0 > max) { max = ${missing_pos}+0; line = $0 } END { print line }'"+ {output}+".txt")
    #cmd = f"awk -F'\t' '$2 ~ /Gn/{ if ($4 > max) { max = $4; line = ${missing_pos}} } END { print $line } ' {output}.txt"
    cmd=f"awk '/{prot}/ && (NR==1 || $4 > max) {{ max = $4; line = ${missing_pos} }} END {{ print line }}' {output}.txt"
   
    pos=subprocess.run(f"{cmd}",capture_output=True,shell=True, text=True)
   
    if pos.stdout=="\n" or pos.stdout=='' : 
        return int(-1)
    else  :
        return int(pos.stdout)

def change_pos(line,end_col,Pos_prot_final) : 
    for col in range (1,end_col) :
        if "Start" in str(Pos_prot_final.columns[col]) : 
            #print(type(Pos_prot_final.iloc[line,col]))
            if "str" not in str(type(Pos_prot_final.iloc[line,col]))   :

                if Pos_prot_final.iloc[line,col]==-1 :
                    Pos_prot_final.iloc[line,col]="Na"
                    Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                elif "str" not in str(type(Pos_prot_final.iloc[line,col+1])) and  Pos_prot_final.iloc[line,col+1]!=-1:
                    if Pos_prot_final.iloc[line,col] > Pos_prot_final.iloc[line,col+1]  : 
                        if col ==1 : 
                            Pos_prot_final.iloc[line,col]=1
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                            #print(f"Start- {Pos_prot_final.iloc[line,col]}={Pos_prot_final.index[line]}")
                        else :
                            #print(f"Start-else {Pos_prot_final.iloc[line,col]}={Pos_prot_final.iloc[line,col-1]} : {Pos_prot_final.iloc[line,0]}")
                            Pos_prot_final.iloc[line,col]=int(Pos_prot_final.iloc[line,col-1])+1
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
            elif "str" not in str(type(Pos_prot_final.iloc[line,col+1])) and "str" not in str(type(Pos_prot_final.iloc[line,col-1])) : 
                Pos_prot_final.iloc[line,col]=int(Pos_prot_final.iloc[line,col-1])+1
                Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                        
        else : 
            if "str" not in str(type(Pos_prot_final.iloc[line,col]))  :

                if Pos_prot_final.iloc[line,col]==-1 :
                    Pos_prot_final.iloc[line,col]="Na"
                    Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                elif "str" not in str(type(Pos_prot_final.iloc[line,col-1])) and  Pos_prot_final.iloc[line,col-1]!=-1:
                    if  Pos_prot_final.iloc[line,col-1] > Pos_prot_final.iloc[line,col] :
                        if col !=end_col : 
                            Pos_prot_final.iloc[line,col]=int(Pos_prot_final.iloc[line,col+1])-1
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                        else : 
                            Pos_prot_final.iloc[line,col]="END"
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
            elif "str" not in str(type(Pos_prot_final.iloc[line,col+1])) and "str" not in str(type(Pos_prot_final.iloc[line,col-1])) : 
                Pos_prot_final.iloc[line,col]=int(Pos_prot_final.iloc[line,col+1])-1
                Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                    #print(f"end {Pos_prot_final.iloc[0,col]}={Pos_prot_final.columns[col]}")
        #print(Pos_prot_final.iloc[line,col])
    
#print(f"{INOUT}\n{protINOUT}")
dicy={}
#print(Allstudies,INOUT,protINOUT,orgaWIWO)
print(Allstudies)
print(INOUT)
print(protINOUT)
print(orgaWIWO)
for all_s, IO,prots,virus in zip(Allstudies,INOUT,protINOUT,orgaWIWO) :
    metadata = pd.read_csv(f"{path}{all_s}/Metadata_{all_s}.txt", delimiter="\t",header=0,index_col=0)
    print(all_s)
    if IO!='': #
        dicy[all_s]={}
        if IO=="OUT" : 
            for v in virus.split(";") : 
                dicy[all_s][v]=[prots]
        else : 
            list_ID=metadata["Nucleic accession number"].tolist()
            for  ID in list_ID :
                if ID not in virus.split(";") :
                    if ID not in dicy[all_s] :
                        dicy[all_s][ID]=[prots]
                    elif prots not in  dicy[all_s][ID] : 
                        dicy[all_s][ID].append(prots)
    elif IO=='none' : 
        dicy[all_s]={}

if not dicy : 
    for all_s in Allstudies : 
        dicy[all_s]=["-","-","-"]
#print(dicy)
#print(dicy.items())
#print(studies,target_list,dicy.items())
for study,target in zip(studies,target_list) : 

    pos_prot = pd.read_csv(f"{path}{study}/Pos_prot_{study}.txt", delimiter="\t",header=0,index_col=0)
    metadata = pd.read_csv(f"{path}{study}/Metadata_{study}.txt", delimiter="\t",header=0,index_col=0)
    pos_prot.iloc[:, 1:]=pos_prot.iloc[:, 1:].astype(str)
    pos_prot.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt",sep="\t")
    Pos_prot_final=pd.read_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", delimiter="\t",header=0,index_col=0)
    #OG_protOUT={}
    #print(pos_prot)
    # for element in dicy[study] : 
    #     print(element)
    #     print(element[1][0],element[1][1],element[1][2],element)
    #     if element[1][0]=="OUT" : 
    #         for virus in element[1][2].split(";") :
    #             OG_protOUT[virus]=element[1][1]
    #     elif element[1][0]=="IN":
    #         list_ID=metadata["Nucleic accession number"].tolist()
       
    #         for  ID in list_ID :
    #             if ID not in element[1][2].split(";") :
    #                 if ID not in OG_protOUT :
    #                     OG_protOUT[ID]=[element[1][1]]
    #                 elif element[1][1] not in  OG_protOUT[ID] : 
    #                     OG_protOUT[ID].append(element[1][1])
    #print(f"OG_protOUT : {OG_protOUT}")
    end_line=metadata.shape[0]
    end_col=metadata.shape[1]
    #print(end_line)
    match_ID={}
    for line in range(0,end_line) : 
        ID_itol=metadata.columns.get_loc('ID itol tree')
        ID_Nu=metadata.columns.get_loc('Nucleic accession number')
        match_ID[metadata.iloc[line,ID_itol]]=metadata.iloc[line,ID_Nu]

    end_line=pos_prot.shape[0]
    end_col=pos_prot.shape[1]  
    missing_val={}
    for line in range(0,end_line) : 
        for col in range(1,end_col-1) :  
            if "nan" in str(pos_prot.iloc[line,col]) : 
                DB = str(f"{path}{study}/{target}_ref")
                File=str(f"{path}{study}/prot/{match_ID[pos_prot.iloc[line,0]]}_prot.fasta")
                output=str(f"{path}{study}/Align/{target}_{match_ID[pos_prot.iloc[line,0]]}")
                #print(output)
                protein=(str(f"{pos_prot.columns[col]}"))
                #print(f"{OG_protOUT}")
                if match_ID[pos_prot.iloc[line,0]] in dicy[study] :
                    #print(f"element {str(pos_prot.columns[col]).lower()}")
                    if not any(str(elem).lower() in str(pos_prot.columns[col]).lower()  for elem in dicy[study][match_ID[pos_prot.iloc[line,0]]]) : 
                        
                        
                        if pos_prot.iloc[line,0] not in missing_val : 
                            missing_val[match_ID[pos_prot.iloc[line,0]]]=[pos_prot.columns[col]]
                            Blastp_run(DB,File,output)
                            
                            new_pos=add_new_pos(protein,output)
                            # if new_pos =="Aucune correspondance" :
                            #     Pos_prot_final.iloc[line,col].astype(str)
                            Pos_prot_final.iloc[line,col]=new_pos
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")

                        elif pos_prot.columns[col] not in missing_val[match_ID[pos_prot.iloc[line,0]]]: 
                            missing_val[match_ID[pos_prot.iloc[line,0]]].append(pos_prot.columns[col])
                            
                            new_pos=add_new_pos(protein,output)
                            # if new_pos =="Aucune correspondance" :
                            #     Pos_prot_final.iloc[line,col].astype(str)
                            Pos_prot_final.iloc[line,col]=new_pos
                            Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                    else : 
                        #print("else")
                        
                        #Pos_prot_final.iloc[line,col].astype(str)
                        Pos_prot_final.iloc[line,col]="ABS"
                        Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                else :  
                    
                    #Blastp_run(DB,File,output)
                    if pos_prot.iloc[line,0] not in missing_val : 
                        missing_val[match_ID[pos_prot.iloc[line,0]]]=[pos_prot.columns[col]]
                        Blastp_run(DB,File,output)

                        new_pos=add_new_pos(protein,output)
                        
                        # if new_pos =="Aucune correspondance" :
                        #     Pos_prot_final.iloc[line,col].astype(str)
                        Pos_prot_final.iloc[line,col]=new_pos
                        Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")

                    elif pos_prot.columns[col] not in missing_val[match_ID[pos_prot.iloc[line,0]]]: 
                        missing_val[match_ID[pos_prot.iloc[line,0]]].append(pos_prot.columns[col])

                        new_pos=str(add_new_pos(protein,output))
                        
                        # if new_pos =="Aucune correspondance" :
                        #     Pos_prot_final.iloc[line,col].astype(str)
                        Pos_prot_final.iloc[line,col]=new_pos
                        Pos_prot_final.to_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", sep="\t")
                        Pos_prot_final=pd.read_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", delimiter="\t",header=0,index_col=0)
        change_pos(line,end_col-1,Pos_prot_final)
        Pos_prot_final=pd.read_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", delimiter="\t",header=0,index_col=0)
        change_pos(line,end_col-1,Pos_prot_final)
    Pos_prot_final=pd.read_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", delimiter="\t",header=0,index_col=0)

    #print(f"{missing_val}")
                
                    
