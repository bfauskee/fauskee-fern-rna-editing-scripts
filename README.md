# fauskee-fern-rna-editing-scripts
Data from Fauskee et al. 2025: Comparative phylogeneyic analyses of RNA editing in fern plastomes suggest possible adaptive innovations
Scripts: = scripts used for analysis, .Rev files are RevBayes scripts used for comparative analyses, RNA_editing_detection.sh maps RNA reads to plastid coding sequences, RNA_editing_detetion.Rmd actually detects the editing sites, interprets their amino acid change, and more
## Scripts
### Birth-death
* **BDP.Rev**: This scripts discribes our Birth-Death Process. We assign exponential priors to Speciation and Extinction rates. Then we assign moves to them so we can sample them in the MCMC. We set rho to the approximate sampling probability (your best guess). This script uses a value of 0.03 since we had 18 samples in the phylogeny out of an estimated 600. Since we are not actually "hard" dating this tree we set an arbitrary root age of 1. The present is 0 so we will essentially just adjust branch lengths within this 1->0 window, yielding informative branch lengths. dbBDP() puts it all together by creating a birth death treee distribution, We assign moves to adjust our branch lengths. We also use bd_tree.setValue(topology) to fix our topology based on some other phylogenetic analysis, in this case a ML estimation using ~80 chloroplast genes.
* **Morph_clock.Rev**: This script defines our morphological clock. The "morphological" traits here are RNA editing sites, which we code as binary. This can be easily adapted for other binary morphological traits as well. We draw the average rate of morphological evolution from an exponential prior. we also set an exponential prior on the standard deviation of the mean morph. change. We then add moves for these. We also create a deterministic variable, the variance of morphological change, which is calculated my squaring the standard deviation. We also assign each branch their specific rates (of morpho change) these are drawn from an uncorrelated log-normal distribution. We then add moves for the branch rates, contained in a for loop looping over each branch. This sets up the clock model, and is used by the master_mcmc script which actually runs the estimation.
* **Mol_clock.Rev**: This defines the molecular clock, here we use just the one chloroplast gene rbcL. It is done nearly the same way as the mohphological clock
* **Morph_model.Rev**: We define our morphological change model, often referred to as Q. The code here translates to a Jukes-Cantor model with gamma-distributed rate heterogeneity, descritized into 4 rate classes. we use fnJC(2) which calls for the Jukes-Cantor model (equal rates and state frequencies) and 2 character states (since this is binary data). We put an exponential prior on alpha. alpha is the shape parameter which determines the shape of the gamma distribution from which we will estimate our rates. Rates then become a determistic node from alpha. We assign moves and clamp our data.
* **Mol_model.Rev**: We define our molecular change model. We are going to use a GTR+I+G model, which means we use the general time-reversible model (unequal base frequencies and unequal exchangeability). Base frequencies are drawn from a Dirichlet distribution as are exchangeability rates. we then add moves to sample them from their respecive Dirichlet distributions. using the fnGTR() function with our site frequencies and ex. rate variables we have defined the GTR portion of our model. We draw a proportion of invariant sites (the +I part) from a beta distribution and add moves to sample. we then similarly draw the gamma shape parameter (alpha) from an exponential distribution. the rate classes are then sampled with the fnDiscretizeGamma() function where we define 4 rate classes. Like before, we assign moves, and clamp the model to rbcL (our molec. data)
* **master_mcmc.Rev**: This is what ultimately runs the MCMC. It can be done by exporting RevBayes into the user's path and running `rb master_mcmc.Rev` or that same command can be contained into a slurm script or something similar if running on a high-performance cluster. We need a few commands here such as reading in our data (molecular and morphological in this case). we also load in a phylogeny since we are fixing the topology. we set topology as a constant variable. We also have to tell Revbayes how many branches are in our tree (2 * the number of taxa -2). We then load in all other relevant scripts: the BDP tree prior, morphological and molecular clock priors, mohpological change model, and molecular substitution model. We then initialize and define out monitors (things that will be printed to the screen, log file, and tree posterior file). All that is left is to define how many independent runs we want to have and how many generations we want them to run for.

### RNA editing detection and characterization
* **RNA_editing_detection.sh** This script uses our chloroplast coding genes, maps RNA to them and outputs a .tsv noting the number of reads mapped to each site in each sample as well as the how many of each nucleotide map to the reference. This tsv output is used with the R script which will identify and characterize the RNA editing sites. This is a linux pipeline. The commands can easily be placed in a slurm (or other scheduler) script for use on a high-performance cluster. This is adapted from Edera and Sanchez-Puerta 2021 with minor modifications
* **RNA_editing_detection.Rmd**: This R markdown script is designed to be run locally (or on a server if needed) on the .tsv output from RNA_editing_detection.sh. It will output several .csv files. Importantly, it outputs one file with all detected RNA edits and includes the type of edit (C-to-U or U-to-C), the codon position of the edit, the codon sequence, the amino acid change produced by the editing event and also include the gene and postion of each edit.  

## edits-only
* Contains outputs form the RNA editing detection for each species. 


