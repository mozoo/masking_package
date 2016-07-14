#!/bin/sh
for i in *.fasta_aln
do /opt/Gblocks_0.91b/Gblocks "$i" -t=p -b2=58 -b3=10 -b4=5 -b5=a -v=80
done
