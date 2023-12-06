#!/bin/bash
##################################################################
# get-protein-stability-data.sh meltomeFile
# queries the cross-species file from Meltome Atlas http://meltomeatlas.proteomics.wzw.tum.de:5003/
# (Jarzab et al., (2020), Nat. Methods.) for Arabidopsis
##################################################################

organismName="Arabidopsis"

outDir=AraTModel/protein-stability
cwd=$(pwd)

if [[ $# == 1 ]];
then
	meltomeFile=$1
else
	meltomeFile=AraTModel/protein-stability/meltome-atlas-cross-species.csv
fi

# processing cross-species table
awk -F"," '{if (NR==1 || $1 ~ name) print $2,$3,$4,$6,$7}' name="${organismName}" OFS="\t" $meltomeFile | sed 's/"//g' > $outDir/ath-protein-stability.csv


# downloading whole data set for A. thaliana
# cd $outDir

#wget -q ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/04/PXD011929/A.thaliana_R1.rar
#wget -q ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/04/PXD011929/A.thaliana_R2.rar
#wget -q ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/04/PXD011929/A.thaliana_R3.rar

# unpacking data sets
#unrar e A.thaliana_R1.rar results_TPP_TR_full_splitted_ids_links.xlsx
#mv results_TPP_TR_full_splitted_ids_links.xlsx ath-meltome-r1.xlsx
#rm A.thaliana_R1.rar

#unrar e A.thaliana_R2.rar results_TPP_TR_full_splitted_ids_links.xlsx
#mv results_TPP_TR_full_splitted_ids_links.xlsx ath-meltome-r2.xlsx
#rm A.thaliana_R2.rar

#unrar e A.thaliana_R3.rar results_TPP_TR_full_splitted_ids_links.xlsx
#mv results_TPP_TR_full_splitted_ids_links.xlsx ath-meltome-r3.xlsx
#rm A.thaliana_R3.rar

# melting curves:
#wget -q ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/04/PXD011929/TPP_A.thaliana_P023760.rar
#unrar e TPP_A.thaliana_P023760.rar TPP_A.thaliana_P023760.xlsx
#rm TPP_A.thaliana_P023760.rar

cd $cwd

