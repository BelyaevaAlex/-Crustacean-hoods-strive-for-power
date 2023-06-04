# The reproducibility of the Crustacean phylogenetic analyses on open source data
---------------------------------------------------------------------------------------

**Authors and their affilations:**
- Alexandra Belyaeva [@BelyaevaAlex](https://github.com/BelyaevaAlex), MSU +BI
- Alexandr Zhuravlev [@Baragozin](https://github.com/Baragozin), PNPI + BI

**Scientific supervisor and her affilation:**
- Polina Drozdova, ISU · Institute of Biology

## Table of contents:
- [Introduction](#sec1) </br>
- [Repository structure](#sec2) </br> 
- [Data](#sec3) </br>
- [Aim and objectives](#sec4) </br>
- [Pipeline of the project](#sec5) </br>
- [Our analysis tools and their versions](#sec6) </br>
- [Results](#sec7) </br>
- [Summary](#sec8) </br> 
- [References](#sec9) </br>

<a name="sec1"></a>
## Introduction

In phylogenetic studies, the problem of reproducibility of data is a major issue. This problem arises due to several factors, including the complexity of the data, the use of different methods and software packages, and the lack of standardization in reporting results. The use of different methods and software packages can lead to different results, which makes it difficult to compare and reproduce results. Additionally, the lack of standardization in reporting results makes it difficult for other researchers to reproduce the same results. 

This project examined the issue of phylogenetic analysis reproducibility in a number of studies on freshwater crustaceans based on open data. 


<a name="sec2"></a>

## Repository structure
- [image](/image) contains figures for README.md you are currently reading
- [checklists](/checklists) in this folder you can find all our checklists in all types
- [Moskalenko](/Moskalenko) contains scripts for analyzing the reproduction of Moskalenko's article [^1]
- [Bystricky](/Bystricky) contains scripts for analyzing the reproduction of Bystricky's article [^2]

NB: There are no folders called Bukin and Watier, since software applications without code are selected for them as analysis tools.

<a name="sec3"></a>

## Data

__1. Moskalenko [^1]: Specimens were collected by the authors of the article in 2012 and 2013 field seasons on Lake Baikal  and preserved in 90% ethanol. A subset of species for sequencing was chosen to represent major lineages within Baikal amphipods, including lineages close to basal nodes. Non-Baikal freshwater palearctic _Gammarus species_ morphologically similar to Baikal amphipods were also considered. Sequences that branched out beyond _Gammarus species_, namely _Hyalella azteca_ and _Dikerogammarus villosus_ were used as outgroups.__

Reproducibly extracting data from NCBI with the help of Entrez Utilities (E-Utils) is possible:

```
from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'Your crustacean email@yandex.ru' 
accession_numbers = ["MN005083", "MN005110", "MN005067", "MN148361", "MT110188",
                    "MN005085", "MN005112", "MN005069", "MN148353", "MT110190",
                    "MN005075", "MN005102", "MN005059", "MN148355", "MT110180", "GERH00000000",
                    "MN005082", "MN005109", "MN005066", "AY926682", "MT110187",
                    "MN005072", "MN005099", "MN148354", "MT110177", "GEQP00000000",
                    "MN005087", "MN005114", "MN005071", "MN148352", "MT110192", "GEPN00000000",
                    "MN005079", "MN005106", "MN005063", "MN148358", "MT110184", "GEQG00000000",
                    "MN005076", "MN005103", "MN005060", "MN148356", "MT110181", "GEQM00000000",
                    "MN005086", "MN005113", "MN005070", "MT110191", "GEQL00000000",
                    "MN005080", "MN005107", "MN005064", "MN148359", "MT110185",
                    "MN005074", "MN005101", "MN005058", "MN148349", "MT110179",
                    "MN005084", "MN005111", "MN005068", "MN148351", "MT110189",
                    "MN005073", "MN005100", "MN005057", "MN148362", "MT110178", "GEPS00000000",
                    "MN005081", "MN005108", "MN005065", "MN148360", "MT110186",
                    "MN005077", "MN005104", "MN005061", "MN148357", "MT110182",
                    "KU056126", "JF966123", "KU056214", "JX899354", "KY378963", "JF966002",
                    "MN005078", "MN005105", "MN005062", "MN148350", "MT110183", "GERD00000000",
                    "AY926725", "AY926847", "AY529073", "MG320683", "KJ721817",
                    "KF521875", "KY618503", "KY618465", "KF521835", "JF966033",
                    "EF582873", "EF582919", "JF965739", "EF570327", "JF966065",
                    "KF824647", "EF582921", "JF965743", "EU146922", "JF966069",
                    "EF582864", "EF582912", "JF965747", "AB893342", "JF966071",
                    "EF582877", "AF202982", "JF965765", "JF965943", "JF966086",
                    "EF582879", "EF582925", "KP789769", "JF965984", "JF966111",
                    "EF582858", "EF582905", "EF582954", "MG734968",
                    "AJ440902", "EF582898", "KF478496", "MG986757",
                    "NC_039403", "AY743944", "DQ464742", "JX446367", "XM_018164572", "XM_018168833"]

query = " OR ".join(accession_numbers)
handle = Entrez.efetch(db="nucleotide", id=query, rettype="fasta", retmode="text")
sequences = handle.read()
```

However, transcriptome assemblies were sometimes mixed among the sequences, which had to be worked with using blastn.

__2. Bystřický [^2]: Individual lines of _G. fossarum_ in their contact zone in the Western Carpathians are being studied.__

There is no source data. But there is an archive already with ready alignments in supporting information: https://aslopubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Flno.12239&file=lno12239-sup-0007-DataS1.zip 

__3. Bukin [^3]: Specimens of the amphipod _G. fasciatus_ were collected in the littoral zone of Lake Baikal + previously published in pop-set data.__

URL-link for pop-set data:
https://www.ncbi.nlm.nih.gov/popset/376403838

The sequences obtained in the study were deposited in the database under the following accession numbers: MG214957, MG214958.

MG214957:
https://www.ncbi.nlm.nih.gov/nuccore/MG214957.1/

MG214958:
https://www.ncbi.nlm.nih.gov/nuccore/MG214958.1/

__4. Wattier [^4]: COI sequences for the DNA-barcode region were either newly generated or derived from the literature for 4926 _Gammarus fossarum_ individuals, from 498 sampling sites distributed throughout 19 European countries.__

All these sequences were derived into three groups.

1. All newly produced COI sequences are released in GenBank under accession numbers: MT978656–MT980726.
https://www.ncbi.nlm.nih.gov/nuccore/?term=MT978656:MT980726[pacc]

2. The 463 sequences from Lagrue et al. which were initially unreleased in Genbank are now available under accession numbers: MT411018-MT411480.
https://www.ncbi.nlm.nih.gov/nuccore/?term=MT411018:MT411480[pacc]

3. Sequences and metadata are available in BOLD dataset DS-GFOSCDEU.
http://www.boldsystems.org/index.php/Public_SearchTerms?query=DS-GFOSCDEU

<a name="sec4"></a>
## Aim and objectives

__Goal:__ The aim of this project was to repeat the phylogenetic research in certain articles and understand about reproducibility of results.

__Objectives:__
1. searching for the vault of data  (literature review, blastn, NCBI database search due to data errors, study of additional materials)
2. repeating them with modern tools
3. comparison of the results of using several tools at a specific step
4. justification of the non-reproducibility of some steps of the original pipeline
5. creation of “must have” list with essential for this kind of phylogenetic analysis 
6. creation of  a list of shortcomings, typos, inconsistencies in the articles reviewed
7. creation of “must have” list with essential for replication checkpoints

<a name="sec5"></a>
## Pipeline of the project

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Pipeline.png) 


<a name="sec6"></a>
## Our analysis tools and their versions

In view of the fact that, in addition to the reproducibility task, it was decided to try, test, and compare as many tools as possible, unfortunately, it is not possible to not create one virtual environment due to the inconsistency of the versions of the required programs. However, if you only want to replicate our research without comparison, then all the necessary programs can be found in the file `Crustacean_hoods.yml` in order to create a reproducible conda environment using this file. To visualize trees, it is better to use the R - package `ggtree`:
To install this package, start R (version "4.3") and enter:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
```

The rest of the programs (MEGA 11, BEAST v2.7.4 ) are installed as regular applications on your computer, the process varies depending on the operating system.

__Data preparation:__

- The GenBank database:  https://www.ncbi.nlm.nih.gov/genbank/ 
- Nucleotide BLAST : https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

__Alignment:__

- MEGA 11 (11.0.10) for macOS: https://www.megasoftware.net  
- Unipro UGENE вер. 46.0 (64-bit version): https://ugene.net/ru/ 
- Jalview (version: 2.11.2.6): https://www.jalview.org/download/macos/ 
- CLUSTAL 2.1 Multiple Sequence Alignments: http://www.clustal.org
- MUSCLE v3.8.1551 by Robert C. Edgar: https://apolo-docs.readthedocs.io/en/latest/software/applications/MUSCLE/index.html
- MAFFT v7.490 (2021/Oct/30): https://mafft.cbrc.jp/alignment/software/
- Kalign (3.3.1): https://github.com/TimoLassmann/kalign
- T-COFFEE Version_13.41.0.28bdc39 (2019-11-30 10:21:53 - Revision 5d5a1c1 - Build 465): https://tcoffee.org/Projects/tcoffee/index.html
- EMBOSS:6.6.0.0: http://en.bio-soft.net/format/emboss.html

__Trimming:__

- trimAl v1.2: http://trimal.cgenomics.org/downloads 

__Data concatenation and partitions:__

- Catfasta2phyml v.1.2.0 : https://github.com/nylander/catfasta2phyml 
- MEGA 11 (11.0.10) for macOS: https://www.megasoftware.net  

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

__Visualization of trees:__
- IQ-TREE multicore version 2.0.3 for Linux 64-bit built Dec 20 2020: http://www.iqtree.org
- ggtree 3.8.0: https://www.bioconductor.org/packages/release/bioc/html/ggtree.html

<a name="sec7"></a>
## Results

### Moskalenko:

**Preparing data for alignment:**

The article describes two new species of _Baikal gammarids_, _Eulimnogammarus trengove sp. nov._ and _Eulimnogammarus tchernykhi sp. nov._, with some morphology indicating a basal position within the radiation, but strangely, with further concatenation of the data, we have not to consider them due to the lack of data provided in the article. In addition, the description of datasets that was presented in the methods and materials section does not contain these two types of _Baikal gammarids_. In addition, a rather interesting fact is that by the type of phylogenetic trees in the article, it seems that the data was still divided into two datasets in a slightly different way (this can be judged due to the difference and the number of species under consideration).

Also to process the data from both the source and the tools obtained during the application, a few small scripts were written to generate the correct input data for the next tools in the queue of our pipeline tools. All of them are presented in the file [Datasets.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Moskalenko/Datasets.ipynb) file. For convenience, they are designed in the form of functions.

**Multiple alignment:**

According to one of the articles, a comparative analysis of the use of tools at different stages of our analysis was carried out. 

6 tools for multiple alignment were used, as a result, the tool (MUSCLE) used by the authors of the article did not show particularly great benefits in its use:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/muscle_compare.png)

Our analysis [Multiple_al.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Moskalenko/Multiple_al.ipynb) showed that kalign tool wins on all indicators based on our data.

**Data concatenation:**

Further, at the stage of data concatenation ([Data_concatenation.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Moskalenko/Data_concatenation.ipynb)), it turned out that a tool like catfasta2phyml is great for getting the output of a data splitting scheme in a file.cfg for PartitionFinder, however, concatenated our data incorrectly, which did not even allow us to trim sequences using the trimAl tool. The MEGA 11 program was chosen for concatenation.

__Dataset 1:__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/gr_1_UG.png)

__Dataset 2:__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/gr_2_UG.png)

__All sequences without gaps in the data:__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/gr_all_UG.png)

**Trimming:**

Although this was not stated in the article, but by the look of our alignments, it was decided not to skip the trimming stage. This step of our pipeline is presented in the [Trimming.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Moskalenko/Trimming.ipynb) file.

**The selection of the evolutionary model:**

In this article, this stage is skipped and it is simply written in the methods and materials section that an ML-tree is being built in MEGA7.

Our selection of the evolutionary model was performed using both ModelTest (ModelTest-NG) and ModelFinder in IQ-TREE ([Evolution_model.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Moskalenko/Evolution_model.ipynb)):

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/model_m.png)

In general, we believe that both tools have come to about the same thing. However, in the article, again, this point was omitted and it is not even known which of the tools was used to select the model, or at least which model of evolution they selected and by what criterion (BIC or AIC).

**Building ML-trees:**

Apart from the fact that the trees were built in MEGA 7, nothing else was said in the article.

Our ML-trees were built in RAxML-NG and in IQ-TREE [ML_trees.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Moskalenko/ML_trees.ipynb), but the optimal solution seemed to be to use IQ-TREE, since it worked faster.

Next, we will present for each of our datasets the trees we obtained using the evolutionary model predicted in ModelFinder in IQ-TREE, but visualized using two different tools: with the help of the library in R ggtree and the library in Python Phylo. In this case, the drawings from the original article are not given, since it is obvious according to the previous theses, the phylogenetic type of trees will be different, since they were built according to different data.

__Dataset 1:__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Group_1_ml.png)

__Dataset 2:__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Group_2_ml.png)

__All sequences without gaps in the data:__
![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Group_all_ml.png)

**Building BA-trees for Dataset 1:**

Phylogenies were also reconstructed by Bayesian analysis (BA) in BEAST2:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Group_1_1_ba.png)

But protein coding sequences (H3) entered as amino acid sequences: 

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Group_1_2_ba.png)

It is very strange that in the original article, the same authors obtained different results on phylogenetic analysis without adding new data, however, none of their articles indicate either evolutionary models for building trees, or parameters set in some tools.

__The main conclusions we have received:__

1) The article described two new species of _Baikal gammarids_, _Eulimnogammarus etingovae sp. nov._ and _Eulimnogammarus tchernykhi sp. nov._, and using the construction of phylogenetic trees, they tried to determine the phylogeny placement of these two new species. However, according to the table and description of datasets provided by the authors of this article, they throw them out of consideration, which is indicated by the type of trees we have built.
2) The hypothesis that _Hyalella azteca_ is an outgroup has been confirmed by ML methodology.
3) According to our analysis, such a unique result was not obtained, as the only species that is consistently attributed to the sister group in relation to the rest of the second Baikal radiation was _Pachyschesis branchialis_. 
4) The hypothesis was confirmed: The inferred relationship of the two Baikal cades with cosmopolitan palearctic _Gammarus species_ strongly depends on the species and genes included and phylogenetic methodology used.
5) Finally, ML analysis suggests that Baikal’s first clade isn't sister to the basal palearctic species which contradicts the results obtained in the article.
6) The lack of agreement between different phylogenetic methods with the low resolution allowed by limited sequence data may be interpreted as an evidence of fast ancient divergence. The hypothesis is confirmed.

### Bystricky:

**Preparing data for alignment:**

One of the main problems of this article is the lack of source data [Dataset.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Bystricky/Dataset.ipynb), we found only alignment files in the appendix to the article (supporting information).

Without having the source data, it is quite difficult to judge how the alignment files were obtained, besides, in the future we will face the difficulty that the files are not concatenated due to the fact that they contain different sequences and in different order.

**Multiple alignment and trimming:**

It is impossible to reproduce these results, the trimming stage was also not described in the article, so in order to avoid any errors, it was not performed.

**Data concatenation:**

This cannot be reproduced using the provided alignment files. In addition, according to the type of phylogenetic trees, they may have been analyzed separately.

**The selection of the evolutionary model:**

We also refuse to use the PartitionFinder tool, since it is quite difficult to install and in addition, it was not even indicated which model was chosen and by what criteria. Therefore, we chose an evolutionary model using two other tools (RAxML-NG and IQ-TREE) by BIC ([Ev_model.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Bystricky/Ev_model.ipynb)):

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/model_b.png)

Apart from the last file, the tools predicted similar models.

**Building ML-trees:**

Next, trees were built using RAxML-NG and in IQ-TREE [ML_trees_b.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Bystricky/ML_trees_b.ipynb), and using the python module from Bio import Phylo, they were visualized.

__16S_alignment.fas:__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/16S_alignment_ml.png)

__28S_alignment.fas:__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/28S_alignment_ml.png)

__EF1a_alignment_phased.fas:__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/EF1a_alignment_phased_ml.png)

In the best quality of the image can be found in the [ML_trees_b.ipynb](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/Bystricky/ML_trees_b.ipynb) file.

For 28S, which is homozygous and lineages exhibit reciprocal monophyly (except CWE D), we applied both a tree-based method. 

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/first_tools.png)

The EF1a marker is heterozygous and lineages are not reciprocally monophyletic. Therefore, tree- or distance-based species delimitation methods are not applicable. Consequently, we used the haploweb method since it is specifically designed to delimit species based on such markers . The haploweb distinguishes putative species clusters based on mutual allelic exclusivity, following Doyle’s view of species as fields for recombination. It was run on the HaplowebMaker web-server (https://eeg-ebe.github.io/HaplowebMaker/):

Species delimitation based on the haploweb method applied to the phased EF1a alleles. Each circle represents a unique haplotype. Haplotypes co-occurring within the same individual are connected by a colored curved line. Each tick mark represents one substitution:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Haplotypes.png)

__The main conclusions we have received:__

1) Maximum likelihood family trees for each analyzed marker. 
2) The same synthetic individuals belonging to different pedigrees were obtained, 
3) The outgroup is everywhere _Gammarus pulex_. 
4) Note that both alleles of the EF1a marker are present in the tree.
5) The haploweb method, applied on EF1a data, confirmed the lack of hybridization since each lineage harbored unique alleles.

### Bukin:

**Multiple alignment:**

The alignment was done using UGENE software Muscle tool. 

**Trimming:**

There was no information about trimming in the original paper. It was decided to make one with -automated1 setting.

Example of the command: __trimal -in RAW_DATA -out DATA_TRIM -automated1__

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/bukin_no_trim.png)

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/bukin_trim.png)

**The selection of the evolutionary model:**

_jmodeltest_ software and the BIC creteria was used for model choice. In the original article the best-fitting GTR model with gamma correction (GTR + G) was used for phylogenetic analysis. During the replication was got another model.

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/1_3.PNG)

**Building ML-trees:**

On the our picture we can see the original clades of Central, Northern and southwestern populations. Little difference may be the result of poor given data. 

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Bukin.PNG)

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/all_in_one_2_trim_lognorma_strict_trees_built.txt.png)

### Watier:

**Multiple alignment:**

The alignment was done using UGENE software Muscle tool. 

**Trimming:**

There was no information about trimming in the original paper. But there was problem with turning on pipeline due to the sample size. So we had to trim four times with -automated1 setting.

Example of the command: __trimal -in RAW_DATA -out DATA_TRIM -automated1__ x4

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/wat_no_trim.png)

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/wat_4_trim.png)

**The selection of the evolutionary model:**

_bmodeltest_ package from BEAST2 software was used for model choice. The best-fitting model of substitution was the transversion model (TVM) with gamma-distributed rate heterogeneity (G) and a given proportion of invariable sites (I): 

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/bmodeltest.PNG)

Received model matched with the original:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/modelIndexes.png)

**Building ML-trees:**

It isn’t possible to compare trees, due to the Jawa heap space error, which hasn’t been solved yet.

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Mistake.PNG)

The result of the original article:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/Watier.PNG)

<a name="sec8"></a>
## Summary

Two of the four articles should be considered complex and poorly reproducible for the following reasons:

1. Lack of code availability: One of the main challenges in reproducing phylogenetic analyses is the unavailability of the code used in the analysis. Without access to the original code, it becomes difficult for other researchers to understand and replicate the analysis. This lack of transparency can hinder the reproducibility of the results;

2. The absence of source data in our phylogenetic analysis studies was a significant barrier to reproducibility. Without access to the original source data, other researchers like us cannot replicate the analysis or verify the findings; 

3. Missing or incomplete data: Sometimes, articles do not provide the full dataset used in the analysis, making it impossible to reproduce the results accurately. Also, information about the source data and explanations to them are sometimes located throughout the article, which is why you have to go back through the previous steps of the pipeline;

4. Incomplete or insufficient method descriptions: Some articles may not provide detailed descriptions of the methods used in the phylogenetic analysis. This lack of information makes it challenging for other researchers to accurately reproduce the analysis. Insufficient method descriptions can include incomplete information about the software versions, parameter settings, and other critical details necessary for replication;

5. Ambiguities in data preprocessing: Preprocessing steps, such as sequence alignment or data filtering, are crucial in phylogenetic analysis. However, the specific details of these steps are often not well-documented or may be subject to interpretation. Without clear instructions or descriptions, reproducing the exact preprocessing steps becomes challenging and can lead to inconsistencies in results.The authors miss the moments of correcting output and, accordingly, input data during pipeline execution;

The first item on the list of shortcomings and the last two items also apply to the other two articles.

Summarizing all of the above, we have prepared the following checklists that can help in reproducing and writing articles on phylogenetic analysis:

**“Must have” list with essential for replication checkpoints:**

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/1.png)

**Our list of shortcomings, typos, inconsistencies in the articles reviewed:**

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/2.png)

**“Must have” list with essential for this kind of phylogenetic analysis:**

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/image/3.png)

<a name="sec9"></a>
## References

[^1]:[Moskalenko, V. N.; Neretina, T. V.; Yampolsky, L. Y. To the origin of Lake Baikal endemic gammarid radiations, with description of two new Eulimnogammarus spp. Zootaxa 4766(3), 457-471 (2020). https://doi.org/10.11646/ZOOTAXA.4766.3.5](https://www.mapress.com/zt/article/view/zootaxa.4766.3.5)
[^2]:[Bukin, Y.S., Petunina, J.V. & Sherbakov, D.Y. The Mechanisms for Genetic Diversity of Baikal Endemic Amphipod Gmelinoides fasciatus: Relationships between the Population Processes and Paleoclimatic History of the Lake. Russ J Genet 54, 1059–1068 (2018). https://doi.org/10.1134/S1022795418090053](https://link.springer.com/article/10.1134/S1022795418090053#citeas)
[^3]:[Bystřický, P.K., Rutová, T., Brož, V., Gajdošová, M., Juračka, P.J., Copilaş-Ciocianu, D. and Petrusek, A. Distribution patterns at different spatial scales reveal reproductive isolation and frequent syntopy among divergent lineages of an amphipod species complex in Western Carpathian streams. Limnol Oceanogr 67, 2796-2808 (2022). https://doi.org/10.1002/lno.12239](https://aslopubs.onlinelibrary.wiley.com/doi/10.1002/lno.12239)
[^4]:[Wattier, R., Mamos, T., Copilaş-Ciocianu, D. et al. Continental-scale patterns of hyper-cryptic diversity within the freshwater model taxon Gammarus fossarum (Crustacea, Amphipoda). Sci Rep 10, 16536 (2020). https://doi.org/10.1038/s41598-020-73739-0](https://www.nature.com/articles/s41598-020-73739-0)



