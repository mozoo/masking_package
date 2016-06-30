#!/bin/sh
cd ..
mkdir T-Coffee
#Uncomment the following line if you would like to eliminate dashes ("-").
#sed -i 's/-//g' *.fas
sed -i 's/*//g' *.fas
#Put X in the following line for aminoacids, N for nucleotides.
sed -i 's/?/X/g' *.fas
cp *.fas ./T-Coffee
cp ./masking_package/t_coffee_cmd.sh ./T-Coffee
cd T-Coffee
sh t_coffee_cmd.sh
sed -i 's/ /_/g' *.fasta_aln
cd ..
mkdir Aliscore
mkdir BMGE
mkdir GBlocks
mkdir Noisy
cp ./T-Coffee/*.fasta_aln ./Aliscore
cp ./T-Coffee/*.fasta_aln ./BMGE
cp ./T-Coffee/*.fasta_aln ./GBlocks
cp ./T-Coffee/*.fasta_aln ./Noisy
cp ./masking_package/aliscore_cmd.sh ./Aliscore
cp ./masking_package/BMGE_cmd.sh ./BMGE
cp ./masking_package/GBlocks_cmd.sh ./GBlocks
cp ./masking_package/noisy_cmd.sh ./Noisy
cd Aliscore
sh aliscore_cmd.sh
cd ..
cd BMGE
sh BMGE_cmd.sh
cd ..
cd GBlocks
sh GBlocks_cmd.sh
cd ..
cd Noisy
sh noisy_cmd.sh
cd ..
mkdir Masking
cd Masking
Rscript ../masking_package/masking.R
mv Rplots.pdf ../masking.pdf
Rscript ../masking_package/masking_agreement.R
mv Rplots.pdf ../masking_agreement.pdf
cd ..
