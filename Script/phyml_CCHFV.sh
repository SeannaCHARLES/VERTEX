#!/bin/bash





echo Protdist CCHFV


{ sed '1d' outfile | awk '$1 ~ /[a-zA-Z]/ {print $1}' | while read i; do
    echo -ne "\t${i}"
done
echo "" ;} > test_protdist.txt


nb_colonnes=$(head -n1 "test_protdist.txt" | awk -F"\t" '{print NF}')

sed '1d' outfile |awk -v max=$nb_colonnes '{
    for (i=1; i<=NF; i++) {
        printf "%s\t", $i;
        count++;
        if (count % max == 0) print "";  # Nouvelle ligne après 10 valeurs
    }
}
END { if (count % max != 0) print ""; }' >> test_protdist.txt