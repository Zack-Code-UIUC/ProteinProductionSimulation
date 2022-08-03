# TASEPnoise

This repository contains original Matlab functions used for our paper on the integrated bacterial gene expression modeling for protein production noise. The functions run simulations and analyze simulated data.

Each folder contains code used for different purposes.

<a href="https://github.com/JacobsWagnerLab/TASEPnoise/tree/master/code_for_gene_expression_distribution">code_for_gene_expression_distributions</a>: This folder contains code for measuring various statistical properties of gene expression, derived from many iterations. For example, use this code to obtain distributions of mRNA and protein numbers, distributions of RNAP headways upon loading and upon completion of transcription, distributions of RNAP headway change during transcription elongation, distributions of RNAP numbers on a DNA template, distributions of RNAP and ribosome loading intervals (transcription/translation initiation rate), and distributions of mRNA lifetimes, all from n = 1000 simulations. Also, use this code to show the temporal fluctuations in protein production from a DNA template. Output files from an examplary simulation are uploaded.

<a href="https://github.com/JacobsWagnerLab/TASEPnoise/tree/master/code_for_gene_expression_dist_no_elongation">code_for_gene_expression_distributions_no_elongation</a>: Same as code_for_gene_expression_distributions, but sans transcription and translation elongation, i.e., an elongation-free model. Output files from an examplary simulation are uploaded.

<a href="https://github.com/JacobsWagnerLab/TASEPnoise/tree/master/code_for_RNAP_trafficking">code_for_RNAP_trafficking</a>: This folder contains code for measuring properties of RNAP trafficking that do not require large number of iterations. For example, use this code to plot RNAP trajectories from an examplary simulation run (or a DNA template) and to examine RNAP headway distributions. Note that translation and mRNA degradation are not simulated in this code.

<a href="https://github.com/JacobsWagnerLab/TASEPnoise/tree/master/code_misc">code_misc</a>: This folder contains miscellaneous code used in the paper. Each code is stand alone and its use is explained in the header.

# general comments
If available, check masterscript.m in each folder before begin.

<a href="https://github.com/JacobsWagnerLab/TASEPnoise/tree/master/code_for_gene_expression_distribution">code_for_gene_expression_distributions</a> and <a href="https://github.com/JacobsWagnerLab/TASEPnoise/tree/master/code_for_gene_expression_dist_no_elongation">code_for_gene_expression_distributions_no_elongation</a> are based on parallel computing (parfor function). We used machines that allowed for 12 workers. If needed, change the number of workers in the line declaring "parallelobject = parpool(12);"
