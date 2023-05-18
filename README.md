# Crustacean hoods strive for power
---------------------------------------------------------------------------------------

Authors and their affilations:
- Alexandra Belyaeva [@BelyaevaAlex](https://github.com/BelyaevaAlex), MSU +BI
- Alexandr Zhuravlev [@Baragozin](https://github.com/Baragozin), PNPI + BI

Scientific supervisor and her affilation:
- Polina Drozdova, ISU · Institute of Biology

## Introduction
------------

## Data
------------
## Aim and objectives:
------------
The aim of this project was to repeat the phylogenetic research in certain articles and understand about reproducibility of results.

__Objectives:__
1. searching for the vault of data  (literature review, blastn, NCBI database search due to data errors, study of additional materials)
2. repeating them with modern tools
3. comparison of the results of using several tools at a specific step
4. justification of the non-reproducibility of some steps of the original pipeline
5. creation of “must have” list with essential for this kind of phylogenetic analysis 
6. creation of  a list of shortcomings, typos, inconsistencies in the articles reviewed
7. creation of “must have” list with essential for replication checkpoints

## Pipeline of the project:
------------
![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Pipeline.png) 

## Our analysis tools and their versions
------------
__Data preparation:__
- The GenBank database:  https://www.ncbi.nlm.nih.gov/genbank/ 
- Nucleotide BLAST : https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
__Alignment + trimming:__
- MEGA 11 (11.0.10) for macOS: https://www.megasoftware.net  
- Unipro UGENE вер. 46.0 (64-bit version): https://ugene.net/ru/ 
- Jalview (version: 2.11.2.6): https://www.jalview.org/download/macos/ 
- CLUSTAL 2.1 Multiple Sequence Alignments
- MUSCLE v3.8.1551 by Robert C. Edgar
- MAFFT v7.490 (2021/Oct/30)
- Kalign (3.3.1)
- T-COFFEE Version_13.41.0.28bdc39 (2019-11-30 10:21:53 - Revision 5d5a1c1 - Build 465)
- EMBOSS:6.6.0.0
- Clustal Omega - 1.2.4 (AndreaGiacomo)
- trimAl v1.2: http://trimal.cgenomics.org/downloads 
__Data concatenation and partitions:__
- Catfasta2phyml v.1.2.0 : https://github.com/nylander/catfasta2phyml 
__Best-fitting evolutionary models:__
- PartitionFinder 2.1.1 : https://github.com/brettc/partitionfinder 
- Modeltest-ng x.y.z : https://github.com/ddarriba/modeltest 
- IQ-TREE multicore version 2.0.3 for Linux 64-bit built Dec 20 2020: http://www.iqtree.org
- jModelTest v2.1.10 : https://github.com/ddarriba/jmodeltest2 
__Building tree + dating:__
- BEAST v2.7.4 : https://www.beast2.org 
- IQ-TREE multicore version 2.0.3 for Linux 64-bit built Dec 20 2020: http://www.iqtree.org
- MEGA 11 (11.0.10) for macOS: https://www.megasoftware.net  
- RAxML-NG v. 0.9.0: https://github.com/amkozlov/raxml-ng 

## Results
------------



## Summary
------------



## References
------------
- [Moskalenko, V. N.; Neretina, T. V.; Yampolsky, L. Y. To the origin of Lake Baikal endemic gammarid radiations, with description of two new Eulimnogammarus spp. Zootaxa 4766(3), 457-471 (2020). https://doi.org/10.11646/ZOOTAXA.4766.3.5](https://www.mapress.com/zt/article/view/zootaxa.4766.3.5)
- [Bukin, Y.S., Petunina, J.V. & Sherbakov, D.Y. The Mechanisms for Genetic Diversity of Baikal Endemic Amphipod Gmelinoides fasciatus: Relationships between the Population Processes and Paleoclimatic History of the Lake. Russ J Genet 54, 1059–1068 (2018). https://doi.org/10.1134/S1022795418090053](https://link.springer.com/article/10.1134/S1022795418090053#citeas)
- [Bystřický, P.K., Rutová, T., Brož, V., Gajdošová, M., Juračka, P.J., Copilaş-Ciocianu, D. and Petrusek, A. Distribution patterns at different spatial scales reveal reproductive isolation and frequent syntopy among divergent lineages of an amphipod species complex in Western Carpathian streams. Limnol Oceanogr 67, 2796-2808 (2022). https://doi.org/10.1002/lno.12239](https://aslopubs.onlinelibrary.wiley.com/doi/10.1002/lno.12239)
- [Wattier, R., Mamos, T., Copilaş-Ciocianu, D. et al. Continental-scale patterns of hyper-cryptic diversity within the freshwater model taxon Gammarus fossarum (Crustacea, Amphipoda). Sci Rep 10, 16536 (2020). https://doi.org/10.1038/s41598-020-73739-0](https://www.nature.com/articles/s41598-020-73739-0)



