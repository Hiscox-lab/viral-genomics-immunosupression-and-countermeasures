#!/bin/bash/
#take alignment file to virus to determine the amino acid change derived from SNPs
#future script developement can be updated to have paths to data directory as opposed to running within the bam folder.


for i in $(ls *.bam | sed 's/.bam//' | uniq)
do
	mkdir -p ../DiversiTools
/home/xiaofeng/sources/perl_5.32.1/bin/perl /home/rebee/projects/hiscox-artic-pipeline/scripts/diversiutils_aa_details.pl -bam ${i}.bam -ref /home/rebee/tools/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta -orfs /home/rebee/projects/hiscox-artic-pipeline/scripts/CodingRegion.txt -stub ../DiversiTools/${i}
done

#move to DiversiTools output directory
echo moving to DiversiTools output directory
cd ../DiversiTools

echo concatenating AA outouts
cat *AA.txt > all.AA.txt

echo running R script to identify dominant mutants
R CMD BATCH /home/rebee/projects/hiscox-artic-pipeline/scripts/dt_dominant_mutations.R

echo script fin
