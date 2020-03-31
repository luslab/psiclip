# psiCLIP
Scripts for the pre-processing and figures from "[PsiCLIP reveals dynamic RNA binding by DEAH-box helicases before and after exon ligation](https://www.biorxiv.org/content/10.1101/2020.03.15.992701v1.abstract)". The data for this paper is accessible from Array Express: E-MTAB-8895. 

## Dependencies
To install everything you need to run this code I recommend creating a conda environment. This is done for you in the pre-processing script, but does require that you have Anaconda or miniconda installed.
```
conda env create -f psiclip_environment.yml
```
Alternatively, ensure you have the packages in psiclip_environment.yml installed.

## Pre-processing
Mapping of demultiplexed fastq files and processing to crosslink sites is performed by a Snakemake workflow. This consists of a Snakefile and config.yaml. The chunks of code are contained in the Snakefile and the sample annotations are in the config.yaml. 
#### I want to reproduce the paper
*Note: This analysis has been tested and optimised for a SLURM computing cluster. The files are small enough that the analysis could feasibly run on a laptop. In which case the command to run Snakemake should be altered to remove the* `--cluster` *flag. For alternative cluster architectures, the cluster paramater for each step would be need to be changed to the correct format.*
1. Clone this repository and run `./RUN-ME-for-preprocessing.sh`. This script will download the raw fastq files from Array Express and the metadata information; it will download the yeast genome and annotation from Ensembl; it will create a conda environment called "psiclip" containing all the programs you need and then run the Snakemake pipeline.
#### I want to modify the code for new substrates or samples
In order to modify this code for new samples, you will need to annotate them as we have for our samples in the config.yaml. If you are using a different pre-mRNA substrate you will need a fasta file with your substrate and the snRNAs to map against and an annotation file with the intron position. You can copy the format we have used for our different substrates. The final step for adding a new substrate is to add a new rule in the Snakemake workflow to generate the index for your substrate, and add a new line to the mapping rule, such that it will know to check for your new substrate. The format is straightforward to copy from the substrates available already.
Note that post-"substratome" mapping is done to the yeast genome in this script, but this could be modified if your sample was in human extract for example.

## Scripts for downstream analysis
The main concept for the downstream analysis is that cDNAs mapping to the substratome are normalised to the total cDNAs that map to the yeast genome, because we found that this behaved a lot like a spike-in. We also gaussian smooth the signal and subtract the FLAG-tagged experimental condition from the non-tagged condition. 

#### I want to reproduce the paper
1. Run the steps for pre-processing.
2. Run relevant script from figure_scripts folder. Ensure you are in your psiclip conda environment/ or have the relevant R packages installed.

## Contact
Charlotte Capitanchik, charlotte.capitanchik@crick.ac.uk if you have questions.
