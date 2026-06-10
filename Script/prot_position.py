import subprocess
import sys 
import pandas as pd


from All_Functions import Blastp_run,add_new_pos,change_pos
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


dicy={}

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
    pos_prot.convert_dtypes()
    Pos_prot_final.convert_dtypes()
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
        change_pos(line,end_col-1,Pos_prot_final,path,study)
        Pos_prot_final=pd.read_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", delimiter="\t",header=0,index_col=0)
        change_pos(line,end_col-1,Pos_prot_final,path,study)
    Pos_prot_final=pd.read_csv(f"{path}{study}/Pos_prot_Final_{study}.txt", delimiter="\t",header=0,index_col=0)

    #print(f"{missing_val}")
                
                    
