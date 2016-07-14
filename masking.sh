#!/bin/sh

cd ..

if [ ! -d T-Coffee ]
	then mkdir T-Coffee
	fi

#Comment the following line out if you don't want to eliminate dashes ("-").
sed -i 's/-//g' *.fas
sed -i 's/*//g' *.fas
#Put X in the following line for aminoacids, N for nucleotides.
sed -i 's/?/X/g' *.fas

cp *.fas ./T-Coffee
cp ./masking_package/t_coffee_cmd.sh ./T-Coffee
cd T-Coffee
sh t_coffee_cmd.sh
sed -i 's/ /_/g' *.fasta_aln
cd ..

if [ ! -d Aliscore ]
	then mkdir Aliscore
	fi

if [ ! -d BMGE ]
	then mkdir BMGE
	fi

if [ ! -d GBlocks ]
	then mkdir GBlocks
	fi

if [ ! -d Noisy ]
	then mkdir Noisy
	fi

if [ ! -d Zorro ]
	then mkdir Zorro
	fi

cp ./T-Coffee/*.fasta_aln ./Aliscore
cp ./T-Coffee/*.fasta_aln ./BMGE
cp ./T-Coffee/*.fasta_aln ./GBlocks
cp ./T-Coffee/*.fasta_aln ./Noisy
cp ./T-Coffee/*.fasta_aln ./Zorro

cp ./masking_package/aliscore_cmd.sh ./Aliscore
cp ./masking_package/BMGE_cmd.sh ./BMGE
cp ./masking_package/GBlocks_cmd.sh ./GBlocks
cp ./masking_package/noisy_cmd.sh ./Noisy
cp ./masking_package/zorro_cmd.sh ./Zorro

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

cd Zorro
sh zorro_cmd.sh
cd ..

if [ -d Masking ]
	then rm Masking/*
	else mkdir Masking
	fi

cd Masking
Rscript ../masking_package/masking.R
mv Rplots.pdf ../masking.pdf
Rscript ../masking_package/masking_agreement.R
mv Rplots.pdf ../masking_agreement.pdf
cd ..
