# psiCLIP
Scripts for the pre-processing and figures from "PsiCLIP reveals the arrangement of DEAH-box helicases before and after exon ligation".
In order to replicate this analysis you will also need to download the starting files from Array Express: E-MATB-XXXXXXXXX

## Dependencies

## Pre-processing
Mapping of demultiplexed fastq files and processing to crosslink sites is performed by a Snakemake workflow. This consists of a Snakefile and config.yaml. The chunks of code are contained in the Snakefile and the sample annotations are in the config.yaml. 
In order to modify this code for new samples, you will need to annotate them as we have for our samples in the config.yaml. If you are using a different pre-mRNA substrate you will need a fasta file with your substrate and the snRNAs to map against and an annotation file with the intron position. You can copy the format we have used for our different substrates. The final step for adding a new substrate is to add a new rule in the Snakemake workflow to generate the index for your substrate, and add a new line to the mapping rule, such that it will know to check for your new substrate. The format is straightforward to copy from the substrates available already.
Note that post-"substratome" mapping is done to the yeast genome in this script, but this could be modified if your sample was in human extract for example.

## Scripts for downstream analysis
The main concept for the downstream analysis is that cDNAs mapping to the substratome are normalised to the total cDNAs that map to the yeast genome, because we found that this behaved a lot like a spike-in. We also gaussian smooth the signal and subtract the FLAG-tagged experimental condition from the non-tagged condition. 

## Contact
Charlotte Capitanchik, charlotte.capitanchik@crick.ac.uk if you have questions.
