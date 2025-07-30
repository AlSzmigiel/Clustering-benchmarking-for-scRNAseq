#!/bin/bash

# mkdir -p /work/Master_Project/raw_data/Deng/data
# wget -O /work/Master_Project/raw_data/Deng/data.tar 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE45719&format=file'

# # Extract the data into the specified folder and unzip
# tar -C /work/Master_Project/raw_data/Deng/data -xvf /work/Master_Project/raw_data/Deng/data.tar
# gunzip /work/Master_Project/raw_data/Deng/data/*
## raw reads
paste $(ls data/GSM111*_zy* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
/work/Master_Project/raw_data/Deng/reads_zy.txt
for file in $(ls data/GSM111*_zy* | sort); do
    # Extract the number after "zy" in the filename
    number=$(echo "$file" | sed -n 's/.*_zy\([0-9]*\).*/\1/p')
    sample_name="zy_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_zy.txt
done
paste $(ls data/GSM111*_early2cell_* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_early2cell.txt
for file in $(ls data/GSM111*_early2cell_* | sort); do
    number=$(echo "$file" | sed -n 's/.*_early2cell_\([0-9]*\).*/\1/p')
    sample_name="early2cell_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_early2cell.txt
done
paste $(ls data/GSM111*_mid2cell_* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_mid2cell.txt
for file in $(ls data/GSM111*_mid2cell_* | sort); do
    number=$(echo "$file" | sed -n 's/.*_mid2cell_\([0-9]*\).*/\1/p')
    sample_name="mid2cell_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_mid2cell.txt
done
paste $(ls data/GSM111*_late2cell_* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_late2cell.txt
for file in $(ls data/GSM111*_late2cell_* | sort); do
    number=$(echo "$file" | sed -n 's/.*_late2cell_\([0-9]*\).*/\1/p')
    sample_name="late2cell_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_late2cell.txt
done
paste $(ls data/GSM111*_4cell_* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_4cell.txt
for file in $(ls data/GSM111*_4cell_* | sort); do
    number=$(echo "$file" | sed -n 's/.*_4cell_\([0-9]*\).*/\1/p')
    sample_name="4cell_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_4cell.txt
done
# 8cell files have a bit different notation
paste $(ls data/*_8cell_*-* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_8cell.txt
for file in $(ls data/*_8cell_*-* | sort); do
    number=$(echo "$file" | sed -n 's/.*_8cell_\([0-9]*\).*/\1/p')
    sample_name="8cell_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_8cell.txt
done
paste $(ls data/GSM111*_16cell_* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_16cell.txt
for file in $(ls data/GSM111*_16cell_* | sort); do
    number=$(echo "$file" | sed -n 's/.*_16cell_\([0-9]*\).*/\1/p')
    sample_name="16cell_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_16cell.txt
done

paste $(ls data/GSM111*_earlyblast_* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_earlyblast.txt

for file in $(ls data/GSM111*_earlyblast_* | sort); do
    number=$(echo "$file" | sed -n 's/.*_earlyblast_\([0-9]*\).*/\1/p')
    sample_name="earlyblast_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_earlyblast.txt
done

paste $(ls data/GSM111*_midblast_* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_midblast.txt

for file in $(ls data/GSM111*_midblast_* | sort); do
    number=$(echo "$file" | sed -n 's/.*_midblast_\([0-9]*\).*/\1/p')
    sample_name="midblast_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_midblast.txt
done

paste $(ls data/GSM111*_lateblast_* | sort) | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > /work/Master_Project/raw_data/Deng/reads_lateblast.txt

for file in $(ls data/GSM111*_lateblast_* | sort); do
    number=$(echo "$file" | sed -n 's/.*_lateblast_\([0-9]*\).*/\1/p')
    sample_name="lateblast_${number}"
    sed -i "0,/reads/s//${sample_name}/" /work/Master_Project/raw_data/Deng/reads_lateblast.txt
done

awk -F"\t" '{if ($1) print $1}' /work/Master_Project/raw_data/Deng/data/GSM1112767_zy2_expression.txt > /work/Master_Project/raw_data/Deng/gene-names.txt

paste  /work/Master_Project/raw_data/Deng/reads_*.txt > /work/Master_Project/raw_data/Deng/deng.txt
paste /work/Master_Project/raw_data/Deng/gene-names.txt /work/Master_Project/raw_data/Deng/deng.txt > /work/Master_Project/raw_data/Deng/deng-reads.txt
#sed -i '1s/^#//' /work/Master_Project/raw_data/Deng/deng-reads.txt
