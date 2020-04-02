# Run all the pre-processing
# Charlotte Capitanchik last edit 30.03.20

# Step 0. Make some directories we need
mkdir pre_processing/yeast-genome

# Step 3. Download yeast genome and annotation
wget ftp://ftp.ensembl.org/pub/release-99/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O pre_processing/yeast-genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz -O pre_processing/yeast-genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz
gunzip pre_processing/yeast-genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip pre_processing/yeast-genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz
echo "Downloaded yeast genome and annotation from Ensembl"

# Step 4. Set up a conda environment with the software you need
conda env create -f psiclip_environment.yml
echo "Created psiclip conda environment"
