# fauskee-fern-rna-editing-scripts
Data from Fauskee et al. 2025: Comparative phylogeneyic analyses of RNA editing in fern plastomes suggest possible adaptive innovations
scripts: = scripts used for analysis, .Rev files are RevBayes scripts used for comparative analyses, RNA_editing_detection.sh maps RNA reads to plastid coding sequences, RNA_editing_detetion.Rmd actually detects the editing sites, interprets their amino acid change, and more
##Scripts
-BDP.Rev: This scripts discribes our Birth-Death Process. We assign exponential priors to Speciation and Extinction rates. Then we assign moves to them so we can sample them in the MCMC. We set rho to the approximate sampling probability (your best guess). This script uses a value of 0.03 since we had 18 samples in the phylogeny out of an estimated 600. Since we are not actually "hard" dating this tree we set an arbitrary root age of 1. The present is 0 so we will essentially just adjust branch lengths within this 1->0 window, yielding informative branch lengths. dbBDP() puts it all together by creating a birth death treee distribution, We assign moves to adjust our branch lengths. We also use bd_tree.setValue(topology) to fix our topology based on some other phylogenetic analysis, in this case a ML estimation using ~80 chloroplast genes.


assemlbies.zip = annotated genbank format plastid assemblies used here

edits-only.zip = spreadsheets with all detected rna edits for all species analyzed here
