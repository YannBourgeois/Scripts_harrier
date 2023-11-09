# Scripts for landscape genetics analyses of the Réunion harrier (_papangue_)

This is a set of scripts that were used for a population genomic study of the Réunion harrier (_Circus maillardi_). The preprint can be found [here](https://www.authorea.com/doi/full/10.22541/au.166011136.60432214).

A [workshop](https://github.com/YannBourgeois/Scripts_harrier/tree/main/Workshop_Popgen_harrier) derived from this study and originally developed for students at the University of Portsmouth is also included.

In each folder there is a README file, which is simply a text file containing the annotated commands that were used for different parts of the analysis. 
These commands often call bash scripts (files ending with '.sh') that can be launched on a SLURM scheduler, but would probably require minor modifications to work on your own resource (e.g. queues names, path to where the data are located, number of cores requested, etc.). Make sure to examine and edit the scripts before launching. 

For landscape genetic analyses, a VCF file, as well as a Quarto (which can be opened in RStudio) and a HTML document are provided to reproduce the various maps and population genetics analyses in the paper.

Feel free to reuse and adapt to your own study system. For any question, you can post under the "Issues" tab.

