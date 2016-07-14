#!/bin/sh
for i in *_aa.fas
	do t_coffee -in=S"$i" -mode psicoffee -proxy -multi_core -n_core 16 -output=fasta_aln
	t_coffee -in=S"$i" -mode expresso -proxy -multi_core -n_core 16 -output=fasta_aln
	t_coffee -in=S"$i" -mode accurate -proxy -multi_core -n_core 16 -output=fasta_aln
	done
