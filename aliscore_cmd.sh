#!/bin/sh
for i in *.fasta_aln
do perl -I ~/Aliscore_v.2.0 ~/Aliscore_v.2.0/Aliscore.02.2.pl -N -r 4753 -w 30 -i "$i"
done
