# RACIPEwb_stemness
The code for RACIPEwb framework to generate simulated gene expressions a nine component stemness gene regulatory network and analyze the network. ( Please note, all code written or modified by me includes a header noting this. Any other code in this github was not written by me, Madeline Galbraith, but was used in the creation of this paper.)

If you use this code please cite: [Huang, B. et al. Decoding the mechanisms underlying cell-fate decision-making during stem cell differentiation by random circuit perturbation. J Roy Soc Interface 17, 20200500 (2020).](https://doi.org/10.1098/rsif.2020.0500)


This code been extended from the original RACIPE framework to include the binding and unbinding of OCT4, SOX2, and the OCT4-SOX2 complex. This code is specifically for the nine-component stemness network which includes - Nanog, Gcnf, Gata6, Pbx1, Klf4, Cdx2, OCT4, SOX2, and OCT4-SOX2.

The original RACIPE code (generalized for any network) can be found https://github.com/simonhb1990/RACIPE-1.0

## The code files
 + stem_originalRACIPE.c generates gene expression data for the nine-component stemness network using the original RACIPE framework
 + stem_racipe_withBinding.c generates gene expression data for the nine-component stemness network using the RACIPEwb framework
 + stem_threhold_for_racipeWB.c generates the thresholds for RACIPEwb

## Compiling RACIPEwb
mpic++ -std=c++11 stem_originalRACIPE.c -o racipeOG

mpic++ -std=c++11 stem_racipe_withBinding.c -o racipeWB

g++ stem_threhold_for_racipeWB.c -o thresholds

## Run RACIPEwb
mpirun -n 10 ./racipeOG og

mpirun -n 10 ./racipeWB wb

./thresholds 

## Folders
+ Results in the paper for the original framework use code and data in "originalFramework" written by co-authors
+ The "og" folder is the data generated by the original RACIPE framework and is used for analysis to compare with the RACIPEwb framework
+ The "sf*" folders have data generated by the RACIPEwb framework with different ranges of binding and unbinding parameters as outline in the SI
    * The "sf3" folder is the RACIPEwb generated data used in the main text and analysis
+ The "analyses" folder has the code used to analyze and compare the RACIPE and RACIPEwb frameworks

