#!/bin/sh
for i in *.fasta_aln
	do zorro $i > $i.zorro
	done
