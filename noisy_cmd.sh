#!/bin/sh
for i in *.fasta_aln
do noisy --cutoff 0.80 --missing ? --nogap --smooth 5 --seqtype P "$i"
done
