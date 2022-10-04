# Comparison of Molnupiravir with PF-07321332 treatment either individually or combination in immune supressed mice infected with SARS-CoV-2 as a model for severe COVID-19.

add author list and affiliations when finalised

## Raw data avaliable at Short Read Archive (SRA): Project Acession Number: PRJNA886870

http://www.ncbi.nlm.nih.gov/bioproject/886870

## Sequencing of throat swabs and respiratory tissue derived from mice using the Nimagen amplicon protocol

1. 116 samples were sequenced in total.
2. Amplicon library prep and QC carried out by Jack Pilgrim and Alex Rzeszutek 
3. Illumina sequencing conducted by Centre for Genomics Research (CGR)
4. Demultiplexing and initial data processing conducted by Richard Gregory and Sam Haldenby (CGR Informatics team)
5. Raw fastq files were processed with the [EasySeq v0.9 pipeline v0.9;](https://github.com/JordyCoolen/easyseq_covid19)
6. Consensus sequences were ran through nextclade to assess quality and to determine the breadth of coverage of the consensus sequences
7. Primer Trimmed bam files were inputted into [DiversiTools](http://josephhughes.github.io/DiversiTools/).
8. entropy.txt and AA.txt files were imported into RStudio to assess nucleotide and amino acid changes.


## scripts and contributions
  
1. nextclade.sh - RPR
2. plot-SH665Y.R RPR, IDB
3. plot-tvts.R RPR, IDB

Data visualisation reference: https://github.com/Hiscox-lab/AGILE-molnupiravir-viral-genomics

## Correspondence

Corresponde to: Rebekah Penrice-Randal (rebee@liverpool.ac.uk), Julian Hiscox (julianh@liverpool.ac.uk) or James Stewart (jpstewar@liverpool.ac.uk) 


