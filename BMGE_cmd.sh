#!/bin/sh
for i in *.fasta_aln
do java -jar /opt/BMGE-1.12/BMGE.jar -i "$i" -t AA -of "$i"_BMGE.fas -op "$i"_BMGE.phy -oh "$i"_BMGE.htm
done
