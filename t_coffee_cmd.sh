#!/bin/sh
for i in atp*.fas
	do t_coffee -in=S"$i" -mode psicoffee -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	t_coffee -in=S"$i" -mode expresso -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	t_coffee -in=S"$i" -mode accurate -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	done
for i in cox*.fas
	do t_coffee -in=S"$i" -mode psicoffee -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	t_coffee -in=S"$i" -mode expresso -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	t_coffee -in=S"$i" -mode accurate -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	done
for i in cyt*.fas
	do t_coffee -in=S"$i" -mode psicoffee -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	t_coffee -in=S"$i" -mode expresso -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	t_coffee -in=S"$i" -mode accurate -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	done
for i in nad*.fas
	do t_coffee -in=S"$i" -mode psicoffee -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	t_coffee -in=S"$i" -mode expresso -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	t_coffee -in=S"$i" -mode accurate -email=federico.plazzi@unibo.it -proxy -multi_core -output=fasta_aln
	done
#for i in rrn*.fas
#	do t_coffee -in S"$i" -mode mrcoffee -email federico.plazzi@unibo.it -proxy -multi_core -output fasta_aln
#	done
#t_coffee -in StRNA.fas -mode mrcoffee -email federico.plazzi@unibo.it -proxy -multi_core -output fasta_aln
