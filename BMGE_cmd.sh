#!/bin/sh
for i in *.fasta_aln
do java -jar ~/BMGE-1.1/BMGE.jar -i "$i" -t AA -of "$i"_BMGE.fas -op "$i"_BMGE.phy -oh "$i"_BMGE.htm -m BLOSUM30 -h 0.75
done
