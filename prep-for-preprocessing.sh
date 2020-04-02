# Run all the pre-processing
# Charlotte Capitanchik last edit 30.03.20

# Step 0. Make some directories we need
mkdir data
mkdir data/fastq
mkdir pre_processing/yeast-genome

# Step 1. Download the metadata from array express
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8895/E-MTAB-8895.sdrf.txt -O metadata.tsv
echo "Step 1 completed: Downloaded metadata from array express"

# Step 2. Download all the fastq into the fastq folder
# Note this is 93 files
files=$(cut -f34 metadata.tsv | tail -n +2)
for f in $files; do
    name=$(basename $f)
    echo "downloading ${name}"
    wget $f -O ${name}
done
rename fastq.gz fq.gz E*fastq.gz
mv *fq.gz data/fastq/.
n_files=$(ls data/fastq | wc -l)
echo "Step 2 completed: Downloaded ${n_files} fastqs"

# Step 3. Download yeast genome and annotation
wget ftp://ftp.ensembl.org/pub/release-99/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O pre_processing/yeast-genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz -O pre_processing/yeast-genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz
gunzip pre_processing/yeast-genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip pre_processing/yeast-genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz
echo "Step 3 completed: Downloaded yeast genome and annotation from Ensembl"

# Step 4. Set up a conda environment with the software you need
conda env create -f psiclip_environment.yml
echo "Step 4 completed: Created psiclip conda environment"

