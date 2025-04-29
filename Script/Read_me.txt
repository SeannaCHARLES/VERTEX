





This script  : 
It is currently working for the peribunyaviridae the nairoviridae and CCHFV strains 
    1. Download all the sequences in the Genbank format
        This step download all the sequences in a Genbank format concatenated all the Genbank file into one and delete the single sequences Genbank file.
    2. Seperate all the nucleic sequences and the proteic sequences into single files
        From the Genbank files the nucleic and proteic sequences are being extract and save in 2 folder M_nucl and M_prot (M that represent the segment)
    3. Create a Metadata_files.txt that contains different informations such as ID nucleic, ID proteic , taxonomy, and the position of the differents proteins.
        This step download all the sequences in the gp format and delete it after
    4.Mark the position of the different proteins in the proteic sequences and the nucleic sequences
        For this step it is important to change the different position of each sequences and to change the target ID.
        You also need to change the startseq 
            The startseq was found by using blast and align the proteic sequences to the nucleic sequences 
    5.Seperate the different proteins and save them 
        one file for each proteins.
    6.
    
        From these files multiple informations such as the taxonomy, , ect are being saved in a Metadata_files.txt
If you you are not using the serogroup : 
If you want to use these script for other species
If you only want to create the itol legend
If you only want to create the Metadata files 




Output files 
Metadata_files.txt
