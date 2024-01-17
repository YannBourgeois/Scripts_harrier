# Scripts for landscape genetics analyses of the Réunion harrier (_papangue_)

This is a set of scripts that were used for a population genomic study of the Réunion harrier (_Circus maillardi_). The preprint can be found [here](https://www.authorea.com/doi/full/10.22541/au.166011136.60432214).

A [workshop](https://github.com/YannBourgeois/Scripts_harrier/tree/main/Workshop_Popgen_harrier) derived from this study and originally developed for students at the University of Portsmouth is also included.

In each folder there is a README file, which is simply a text file containing the annotated commands that were used for different parts of the analysis. 
These commands often call bash scripts (files ending with '.sh') that can be launched on a SLURM scheduler, but would probably require minor modifications to work on your own resource (e.g. queues names, path to where the data are located, number of cores requested, etc.). Make sure to examine and edit the scripts before launching. 

### Running Stacks on GBS data
You can find instructions in the folder [Scripts_stacks](https://github.com/YannBourgeois/Scripts_harrier/tree/main/Scripts_stacks)
A short report provided by the sequencing facility can be found in the Data_description.pdf file.
Quality statistics obtained with FASTQC can be found in the folder [Fastqc_reports](https://github.com/YannBourgeois/Scripts_harrier/tree/main/Fastqc_reports)

### Calling variants in RNA-seq data
Instructions to run GATK, call and annotate variants on RNA-seq data are in [this folder](https://github.com/YannBourgeois/Scripts_harrier/tree/main/Scripts_RNAseq)

### Running population genetics and landscape genetics analyses on GBS data

Scripts for population genetic analyses can be found in [this folder](https://github.com/YannBourgeois/Scripts_harrier/tree/main/Landscape_genetics_scripts)
A VCF file, as well as a Quarto (which can be opened in RStudio) and a HTML document are provided to reproduce the various maps and population genetics analyses in the paper.

We also ran a LD-based analysis using the reference genome of _Accipiter gentilis_.
The scripts to align reads to the reference genome, call variants with freebayes, and reconstruct recent changes in effective population sizes with GONE can be found in [this folder](https://github.com/YannBourgeois/Scripts_harrier/tree/main/Scripts_GONe_alignment_RAD_with_reference)

### Running analyses on RNA-seq data (mutation load)

You will find in [this folder](https://github.com/YannBourgeois/Scripts_harrier/tree/main/Scripts_SliM_mutation_load) a set of scripts to estimate PiN/PiS as well as SliM scripts to simulate the evolution of mutation load. 



Feel free to reuse and adapt to your own study system. For any question, you can post under the "Issues" tab.

