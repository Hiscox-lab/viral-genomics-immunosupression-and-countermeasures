#install nextclade with conda install nextclade
#assigns lineages but also provides qc of consensus sequences and amino acid changes and deletions
#future tasks may be to include this in the pipeline
#can use my paths for references and inputs
#only variable that needs changing is input fasta (will do getopts shortly)
#optional variables could be output-basename and output (where this is changed to a specific path) 
#usage for this script is path/to/script/nextclade.sh path/to/consensus/fasta path/to/desired/output/

fasta=$1
output=$2

mkdir -p $output

nextclade --input-fasta=${fasta} --input-root-seq=/home/rebee/tools/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta --genes=E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --input-gene-map=/home/rebee/tools/nextclade/data/sars-cov-2/genemap.gff --input-tree=/home/rebee/tools/nextclade/data/sars-cov-2/tree.json --input-qc-config=/home/rebee/tools/nextclade/data/sars-cov-2/qc.json --output-json=${output}/nextclade.json --output-csv=${output}/nextclade.csv --output-tsv=${output}/nextclade.tsv --output-tree=${output}/nextclade.auspice.json --output-dir=${output}/ --output-basename=nextclade
