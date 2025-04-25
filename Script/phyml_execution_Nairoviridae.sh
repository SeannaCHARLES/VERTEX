#!/bin/bash
path="/home/seannacharles/Documents/"


    nb_line=$(wc -l < "${path}Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt")

    for (( i=1; i<=nb_line/2; i++ ))
                                       do
        
        ap=$(awk -F":" -v i="$i" 'NR==i {print $2}' "${path}Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt")
        if [[ -z "$ap" ]]; then
            protein=$(awk -F":" -v i="$i" 'NR==i {print $1}' "${path}Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt" | awk '{print $1}' )
        else
            model=$(awk -F":" -v i="$i" 'NR==i{print $2}' "${path}Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt"| awk '{print $1}')
        fi
        
        
        if [[ -n "$model" ]]; then
            echo Execution of iqtree2 for the protein "${protein}" of the study Nairoviridae with the model "${model}"
            echo s-"${path}"Nairoviridae/seg_M_phylo/"${protein}"_sdbl.phy d- aa m-"${model}" f- e a- e c- 4 b- 1000 
            iqtree2 -s "${path}"Nairoviridae/seg_M_phylo/"${protein}"_sdbl.phy -d aa -m "${model}" -B 1000
            model=""

        fi
    done

    for (( i=(nb_line/2)+1; i<=nb_line; i++ ))
                                        do
        
        ap=$(awk -F":" -v i="$i" 'NR==i {print $2}' "${path}Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt")
        if [[ -z "$ap" ]]; then
            protein=$(awk -F":" -v i="$i" 'NR==i {print $1}' "${path}Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt" | awk '{print $1}' )
        else
            model=$(awk -F":" -v i="$i" 'NR==i{print $2}' "${path}Nairoviridae/seg_M_phylo/Best_model_Nairoviridae.txt"| awk '{print $1}')
        fi
        
        
        if [[ -n "$model" ]]; then
            echo Execution of iqtree2 for the protein "${protein}" of the study Nairoviridae with the model "${model}" 
            echo s-"${path}"Nairoviridae/seg_M_phylo/"${protein}"_sdbl.phy d- nt m-"${model}" f- e a- e c- 4 b- 1000 
            #phyml -i "${path}"Nairoviridae/seg_M_phylo/"${protein}"_sdbl.phy -d nt -m"${model}" -f e -a e -c 4 -b 100 
            iqtree2 -s "${path}"CCHFV/seg_M_phylo/"${protein}"_sdbl.phy -d nt -m "${model}" -B 1000
            model=""
            
        fi
    done
echo "done" > /home/seannacharles/Documents/scriptNairoviridae_done.txt

