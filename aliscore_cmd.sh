#!/bin/sh
for i in *.fasta_aln
do perl -I /opt/Aliscore_v.2.0 /opt/Aliscore_v.2.0/Aliscore.02.2.pl -N -r 6555 -w 30 -i "$i"
done
