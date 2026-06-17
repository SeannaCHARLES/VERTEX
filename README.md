# VERTEX

![Workflow of the vertex](images/Workflow_v5.pdf)

## Project description

The purpose of this pipeline is to generate phylogenetic trees so that you can analyse and interpret them. This script allows you to generate phylogenetic trees, captions of these phylogenetic trees.

### Execution.sh

- Download all the sequences in the Genbank format.

This step download all the sequences in a Genbank format concatenated all the Genbank file into one and delete the single sequences Genbank file.

+---------------------+------------------+
| Input               | Output           |
+:===================:+:================:+
| All_information.txt | All_sequences.gb |
|                     |                  |
| (From user)         |                  |
+---------------------+------------------+

: Input-Output Execution.sh

### Genbank_extract.py

- Seperate all the nucleic sequences and the proteic sequences into single files

From the Genbank file the nucleic and proteic sequences are being extract and save in 2 folders nucl and prot.

Create a **Metadata_files.txt** that contains different informations.

This step create a file that contains all the informations available in the genbank file format and the Gp file format (the gp file format are being download and concatenate after in the ***prot_all_Peribunyaviridae_Rpacu.gp*** file). This file contains the following information : Study–Virus name–ID itol tree–Taxonomy–Serogroup–Protein accesion number–Size Proteic sequence–Startproteins–Endproteines–Nucleic accession number–Size Nucleic sequence–Strand–Strain–Geographical region–Host or vector References

- Create the caption for the Serogroup, Geographical region, and Host.

The caption are created based on the Itol format.

One color will be associated with the same values for each studies.

The Serogroup captions are created using the Serogroup.txt file. You have to provide you own Serogroup.txt file since this information cannot be found in the genbank files format

+---------------------+---------------------------------------+
| Input               | Output                                |
+:===================:+:=====================================:+
| All_sequences.gb    | proteic_sequence.fasta                |
|                     |                                       |
| (From Execution.sh) | nucleic_sequence.fasta                |
|                     |                                       |
|                     | All_proteic_sequences.fasta           |
|                     |                                       |
|                     | All_nucleic_sequences.fasta           |
|                     |                                       |
|                     | All_proteic_sequences_align.aln-fasta |
|                     |                                       |
|                     | All_nucleic_sequences_align.aln-fasta |
|                     |                                       |
|                     | All_proteic_sequences.gp              |
|                     |                                       |
|                     | **Serogroup_caption.txt**             |
|                     |                                       |
|                     | **Geographical_caption.txt**          |
|                     |                                       |
|                     | **Host_caption.txt**                  |
|                     |                                       |
|                     | **Metadata.txt**                      |
|                     |                                       |
|                     | Pos_prot.txt                          |
+---------------------+---------------------------------------+

: Input-Output Genbank.py

### Prot_position.py

- This file allows the reference proteins to be aligned(using blastp) with proteins whose positions are missing.

This file ensures that as many protein positions as possible are retrieved in cases where the information is not available in the Gp file format.

+-------------------+------------------------+
| Input             | Output                 |
+:=================:+:======================:+
| Pos_prot.txt      | **Final_pos_prot.txt** |
|                   |                        |
| (From Genbank.py) |                        |
+-------------------+------------------------+

: Input-Output Prot_position.py

### Execution.sh

- Create a Blast database of all the reference virus/organisms proteins

+---------------------------+----------------------------+
| Input                     | Output                     |
+:=========================:+:==========================:+
| Ref_prot.fasta            | Proteic_prot_1-n_ref.fasta |
|                           |                            |
| (From Genbank_extract.py) |                            |
+---------------------------+----------------------------+

: Input-Output Execution.sh

### DNA_dist.sh

- This file lauch DNA dist from phylip and create a divergence sequence table.

This file used

+-----------------------------------------+----------------------------------+
| Input                                   | Output                           |
+:=======================================:+:================================:+
| All \_nucleic_sequences_align.aln-fasta | **Divergence_nucleic_table.txt** |
|                                         |                                  |
| (From Execution.sh)                     |                                  |
+-----------------------------------------+----------------------------------+

: Input-Output DNA_dist.sh

### Delete_dbl.py

- This file delete de sequences that are too close using the **DNA_dist.sh** output table

+------------------------------+----------------------------------+
| Input                        | Output                           |
+:============================:+:================================:+
| Divergence_nucleic_table.txt | All_nucleic_sequences_sdbl.fasta |
|                              |                                  |
| (From DNA_dist.sh)           | All_proteic_sequences_sdbl.fasta |
+------------------------------+----------------------------------+

: Input-Output Delete_dbl.py

### Test_script_mark_spe.py

- Seperate the different proteins and save them

One file for each proteins.

+----------------------------------+------------------------+
| Input                            | Output                 |
+:================================:+:======================:+
| All_nucleic_sequences_sdbl.fasta | Nucleic_prot_1-n.fasta |
|                                  |                        |
| (From Delete_dbl.py)             | Proteic_prot_1-n.fasta |
|                                  |                        |
| All_proteic_sequences_sdbl.fasta | Nucleic_prot_OG.fasta  |
|                                  |                        |
| (From Delete_dbl.py)             | Proteic_prot_OG.fasta  |
+----------------------------------+------------------------+

: Input-Output Test_sript_mark_spe.py

### Merge_spe_prot.py

- This file merge the specific OG species, with the specific study.

+-------------------------------+------------------------+
| Input                         | Output                 |
+:=============================:+:======================:+
| Nucleic_prot_1-n.fasta        | Nucleic_prot_spe.fasta |
|                               |                        |
| (From Test_sript_mark_spe.py) | Proteic_prot_spe.fasta |
|                               |                        |
| Proteic_prot_1-n.fasta        |                        |
|                               |                        |
| (From Test_sript_mark_spe.py) |                        |
|                               |                        |
| Nucleic_prot_OG.fasta         |                        |
|                               |                        |
| (From Test_sript_mark_spe.py) |                        |
|                               |                        |
| Proteic_prot_OG.fasta         |                        |
|                               |                        |
| (From Test_sript_mark_spe.py) |                        |
+-------------------------------+------------------------+

: Input-Output Merge_spe_prot.py

### Protdist.py

- This file lauch Protdist from phyllip and create a sequences divergence table

+-------------------------------+---------------------------------+
| Input                         | Output                          |
+:=============================:+:===============================:+
| Nucleic_prot_1-n.fasta        | **Divergence_proteic_table.py** |
|                               |                                 |
| (From Test_sript_mark_spe.py) |                                 |
|                               |                                 |
| Proteic_prot_1-n.fasta        |                                 |
|                               |                                 |
| (From Test_sript_mark_spe.py) |                                 |
|                               |                                 |
| Nucleic_prot_spe.fasta        |                                 |
|                               |                                 |
| (From Merge_spe_prot.py)      |                                 |
|                               |                                 |
| Proteic_prot_spe.fasta        |                                 |
|                               |                                 |
| (From Merge_spe_prot.py)      |                                 |
+-------------------------------+---------------------------------+

: Input-Output Protdist.py

### Execution.sh

- Trimal

- IQtree2

Identified the best model (using iqtree2)

iqtree2 was used because it can provide a test to identify the best model for the nucleic sequences and the proteic sequences. If you only want to study the proteic sequences you may used prottest instead.

- Create the phylogenetic tree using the best model identified (using iqtree2)

+-------------------------------+-------------------------------+
| Input                         | Output                        |
+:=============================:+:=============================:+
| Proteic_prot_1-n.fasta        | **Proteic_prot_1-n.treefile** |
|                               |                               |
| (From Test_sript_mark_spe.py) | **Proteic_prot_spe.treefile** |
|                               |                               |
| Proteic_prot_spe.fasta        |                               |
|                               |                               |
| (From Merge_spe_prot.py)      |                               |
+-------------------------------+-------------------------------+

: Input-Output Protdist.py

## Instructions to use

Informations_file.txt

This file is the Execution input file. You need to fill up all the information in this file in order to use the pipeline.

If an information is missing it will not work!

**DO NOT USE "," AS A SEPARATOR INSIDE EACH COLONE USE ";" INSTEAD**

Example for the study Orthobunyaviridae with a specific outgroup

A : Execution Project name

Ce nom sera utilisé comme nom de dossier. Il vous permet également d'indiquer à quelle projet sont rattaché les outgroup specific le cas écheant.

Si vous voulez étudier plusieurs ségments ou polyprotéines d'une même famille ou espèce vous devez saisir deux noms de projet différents, sans quoi tous les fichiers produit par le pipeline seront pour les deux projets seront dans le même dossier.

+-------------------+
| Project name      |
+===================+
| Orthobunyaviridae |
+-------------------+
| Orthobunyaviridae |
+-------------------+
| Orthobunyaviridae |
+-------------------+

: Project name

B : Study name (if Specific OG)

Ce nom est utilisé pour différencier les outgroups specifique.

Si la ligne ne concerne pas un OG spécifique saisissez le nom de projet. Dans le cas contraire saisir une dénomination pour l'outgroup specifique( son nom par exemple)

Si vous voulez utilisez different outgroup spécifique pour le même projet.

Vous devez avoir une ligne par outgroup.

+-------------------+-------------------+
| A                 | Study name        |
+===================+===================+
| Orthobunyaviridae | Khurdun           |
+-------------------+-------------------+
| Orthobunyaviridae | Pacuvirus         |
+-------------------+-------------------+
| Orthobunyaviridae | Orthobunyaviridae |
+-------------------+-------------------+

: Study name : Where Khurdun and Pacuvirus will be used as outgroups

C : List of ID

La liste des ID du projet pour la ligne projet et les outgroup. Il est possible d'utiliser plusieurs organisme comme outgroup spécifique.

+-------------------+-------------------+-----------------------------------------+
| A                 | B                 | List of ID                              |
+===================+===================+=========================================+
| Orthobunyaviridae | Khurdun           | NC_078560.1                             |
+-------------------+-------------------+-----------------------------------------+
| Orthobunyaviridae | Pacuvirus         | NC_077871.1;NC_055187.1;NC_055410.1;... |
+-------------------+-------------------+-----------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | ID1;IDn                                 |
+-------------------+-------------------+-----------------------------------------+

: List of ID

D : Outgroup Spe (Y/N)

Y : if your project have a specific outgroup for one or multiple proteines.

N : If your project doesn't have a specific outgroup.

+-------------------+-------------------+-------------+--------------+
| A                 | B                 | C           | Outgroup Spe |
+===================+===================+=============+==============+
| Orthobunyaviridae | Khurdun           | IDk         | Y            |
+-------------------+-------------------+-------------+--------------+
| Orthobunyaviridae | Pacuvirus         | IDp         | Y            |
+-------------------+-------------------+-------------+--------------+
| Orthobunyaviridae | Orthobunyaviridae | IDn         | N            |
+-------------------+-------------------+-------------+--------------+

: Outgroup Specific : Is it a specifc outgroup?

E : For wich prot

Quels sont les protéines concernées par cette étude. Si cette ligne concerne la ligne du projet saisissez le nom des protéines des espèces étudiés.

Si cette ligne concerne les outgroups spécifiques indiquer pour quelles protéines cet outgroup sera utilisé

+-------------------+-------------------+---------------+-----------------------+
| A                 | B                 | D             | For wich prot         |
+===================+===================+===============+=======================+
| Orthobunyaviridae | Khurdun           | Y             | nsmgn                 |
+-------------------+-------------------+---------------+-----------------------+
| Orthobunyaviridae | Pacuvirus         | Y             | polyprotein;gn;nsm;gc |
+-------------------+-------------------+---------------+-----------------------+
| Orthobunyaviridae | Orthobunyaviridae | N             | polyprotein;gn;nsm;gc |
+-------------------+-------------------+---------------+-----------------------+

: For wich prot

F : Proteins name (in order)

Proteins name of the study name

+-------------------+-------------------+---------------+--------------------------+
| A                 | B                 | E             | Proteins name (in order) |
+===================+===================+===============+==========================+
| Orthobunyaviridae | Khurdun           | nsmgn         | polyprotein;gn;gc        |
+-------------------+-------------------+---------------+--------------------------+
| Orthobunyaviridae | Pacuvirus         | prots         | polyprotein;gn;nsm;gc    |
+-------------------+-------------------+---------------+--------------------------+
| Orthobunyaviridae | Orthobunyaviridae | prots         | polyprotein;gn;nsm;gc    |
+-------------------+-------------------+---------------+--------------------------+

: Proteins name (in order)

G : spe prot

Saisir les protéines qui seront impliqués à la création de la protéine spécifique

+-------------------+-------------------+---------------+---------------+
| A                 | B                 | F             | spe prot      |
+===================+===================+===============+===============+
| Orthobunyaviridae | Khurdun           | prots         | nsmgn=gn      |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Pacuvirus         | prots         |               |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Orthobunyaviridae | prots         | nsmgn=nsm+gn  |
+-------------------+-------------------+---------------+---------------+

: spe prot

H : TargetID

Saisir un ID d'une organisme qui servira de référence pour l'alignement blast.

+-------------------+-------------------+---------------+---------------+
| A                 | B                 | G             | RefID         |
+===================+===================+===============+===============+
| Orthobunyaviridae | Khurdun           | nsmgn=gn      |               |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Pacuvirus         |               |               |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Orthobunyaviridae | nsmgn=nsm+gn  | NC_001926.1   |
+-------------------+-------------------+---------------+---------------+

: RefID

I : Target proteins position

+-------------------+-------------------+-------------+--------------------------------+
| A                 | B                 | H           | Target /OG proteins position   |
+===================+===================+=============+================================+
| Orthobunyaviridae | Khurdun           |             | 23;292;526;944                 |
+-------------------+-------------------+-------------+--------------------------------+
| Orthobunyaviridae | Pacuvirus         |             | 1;1458;23;305;308;476;598;1386 |
+-------------------+-------------------+-------------+--------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | IDref       | 22;302;303;476;509;1372        |
+-------------------+-------------------+-------------+--------------------------------+

: RefID

J : Start of the proteins in the nucleic sequence

+-------------------+-------------------+------------+-----------------------------------------------+
| A                 | B                 | I          | Start of the proteins in the nucleic sequence |
+===================+===================+============+===============================================+
| Orthobunyaviridae | Khurdun           | pos        | 34                                            |
+-------------------+-------------------+------------+-----------------------------------------------+
| Orthobunyaviridae | Pacuvirus         | pos        | 53                                            |
+-------------------+-------------------+------------+-----------------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | pos        | 57                                            |
+-------------------+-------------------+------------+-----------------------------------------------+

: Start of the proteins in the nucleic sequence

K : DNAdist Threshold

+-------------------+-------------------+---------------+-------------------+
| A                 | B                 | J             | DNAdist Threshold |
+===================+===================+===============+===================+
| Orthobunyaviridae | Khurdun           | 34            | 0.01              |
+-------------------+-------------------+---------------+-------------------+
| Orthobunyaviridae | Pacuvirus         | 53            | 0.05              |
+-------------------+-------------------+---------------+-------------------+
| Orthobunyaviridae | Orthobunyaviridae | 57            | 0.01              |
+-------------------+-------------------+---------------+-------------------+

: DNAdist Threshold

L : Organism to keep ID

+-------------------+-------------------+-------------+---------------------------------+
| A                 | B                 | k           | Organism to keep ID             |
+===================+===================+=============+=================================+
| Orthobunyaviridae | Khurdun           | 0.01        |                                 |
+-------------------+-------------------+-------------+---------------------------------+
| Orthobunyaviridae | Pacuvirus         | 0.05        |                                 |
+-------------------+-------------------+-------------+---------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | 0.01        | MH484340.1;NC_005775;EF990927.1 |
+-------------------+-------------------+-------------+---------------------------------+

: ID Organism to keep

M : Organism with.out a particular protein ID

+-------------------+-------------------+------------+-------------------------------------------+
| A                 | B                 | L          | Organism with.out a particular protein ID |
+===================+===================+============+===========================================+
| Orthobunyaviridae | Khurdun           |            | OUT\*NSM=NC_078560.1                      |
+-------------------+-------------------+------------+-------------------------------------------+
| Orthobunyaviridae | Pacuvirus         |            | NONE\*NN=NN                               |
+-------------------+-------------------+------------+-------------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | IDs        | NONE\*NN=NN                               |
+-------------------+-------------------+------------+-------------------------------------------+

: Organism with.out a particular protein ID

N : Trimal

+-------------------+-------------------+---------------+---------------+
| A                 | B                 | M             | Trimal        |
+===================+===================+===============+===============+
| Orthobunyaviridae | Khurdun           | OUT\*NSM=ID   | Yes           |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Pacuvirus         | NONE          | Yes           |
+-------------------+-------------------+---------------+---------------+
| Orthobunyaviridae | Orthobunyaviridae | NONE          | Yes           |
+-------------------+-------------------+---------------+---------------+

: Organism with.out a particular protein ID

Utiliser trimal pour cette étude?

O : Different prot name

+-------------------+-------------------+-----------+-----------------------------------------------+
| A                 | B                 | N         | Different prot name                           |
+===================+===================+===========+===============================================+
| Orthobunyaviridae | Khurdun           | Y         | g1=gc;g2=gn;glycoproteinprecursor=polyprotein |
+-------------------+-------------------+-----------+-----------------------------------------------+
| Orthobunyaviridae | Pacuvirus         | Y         | g1=gc;g2=gn;glycoproteinprecursor=polyprotein |
+-------------------+-------------------+-----------+-----------------------------------------------+
| Orthobunyaviridae | Orthobunyaviridae | Y         | g1=gc;g2=gn;glycoproteinprecursor=polyprotein |
+-------------------+-------------------+-----------+-----------------------------------------------+

: Organism with.out a particular protein ID

### Metadata_files.txt

This file contains the following informations(if they are available in the NCBI genbank file format and the gb file format) on all the organims you choose.

``` txt
    Virus name  ID itol tree    Taxonomy    Serogroup   Protein accesion number Size
    Proteique sequence Start MLD   End MLD Sart GP38   End GP38    Start Gn    End Gn 
    Start NSm   End NSm Start Gc    End Gc  Nucleique accession number  Size Nucelique
    sequence Strand  Strain  Geographical region Host    Host source References
```

## Additional Instructions and warnings

### Proteins positions

Please pay attention to the ***Final_pos_prot.txt*** file. The blast alignement is only an estimation ce n'est pas une solution complétement fiable. Il est fortement déconseillé d'utiliser vertex si vous n'avez que très peu d'organisme avec des position dans les bases de données.

Il est également fortement recommandé de vérifier la cohérence entre les positions réelle de protéine et les positions prédite. Le fichier ***Pos_prot.txt***contient toutes les positions extraites des base de données.

Si les positions ne vous conviennent pas vous pouvez les modifier manuellement. Après quoi il vous suffit de relancer le pipeline à partir du ***script_mark_spe.py.***

### Region name in the database

A same proteins may have different region name d'un organismes à l'autre.

Si vous avez très peu de positions pour une protéines cela est peut-être dû à une "non reconnaissance du nom de la protéine.

You can use the file ***prot_all_Peribunyaviridae_Rpacu.gp*** to help you identified the different region name.

Que ce soit différent par rapport au nom que vous avez saisie dans les différents nom des proteines red rectangles.

![Gc and Gn different name](images/Diff_G1_G2.pdf)

ou que soit différent par rapport aux différents noms de la même protéine. blue rectangles

![Polyprotein different name](images/diff_nom_polyprot.pdf)

### Host/vectors caption
Images of the various orders of hosts are available on another public GitHub repository. Feel free to add new illustrations, this GitHub repository does not yet contain all the orders. The images were generated by AI.

## Pré-requis

python iqtree2 pipeline python phyml \## Installation \## Output files


