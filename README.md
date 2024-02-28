# The effect of molnupiravir and nirmatrelvir on SARS-CoV-2 genome diversity in infected and immune suppressed mice

Rebekah Penrice-Randal<sup>1*†</sup>, Eleanor G. Bentley<sup>1*</sup>, Parul Sharma<sup>1</sup>, Adam Kirby<sup>1</sup>, I’ah Donovan-Banfield<sup>1,2</sup>, Anja Kipar<sup>1,3</sup>, Daniele F. Mega<sup>1</sup>, Chloe Bramwell<sup>1,4</sup>, Joanne Sharp<sup>4</sup>, Andrew Owen<sup>4,5</sup>, Julian A. Hiscox<sup>1,2,6</sup> and James P. Stewart<sup>1,5</sup>.


<sup>1</sup>Department of Infection Biology and Microbiomes, University of Liverpool, Liverpool, UK.

<sup>2</sup>NIHR Health Protection Research Unit in Emerging and Zoonotic Infections, Liverpool, UK.

<sup>3</sup>Laboratory for Animal Model Pathology, Institute of Veterinary Pathology, Vetsuisse Faculty, University of Zurich, Switzerland.

<sup>4</sup>Department of Pharmacology and Therapeutics, University of Liverpool, UK.

<sup>5</sup>Centre of Excellence in Long-acting Therapeutics (CELT), University of Liverpool, UK.

<sup>6</sup>A*STAR Infectious Diseases Laboratories (A*STAR ID Labs), Agency for Science, Technology and Research (A*STAR), Singapore.



<sup>*</sup>These authors contributed equally.

<sup>†</sup> Corresponding author: rebee@liverpool.ac.uk 



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


