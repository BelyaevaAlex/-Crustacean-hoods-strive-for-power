# The reproducibility of the Crustacean phylogenetic analyses on open source data
---------------------------------------------------------------------------------------

Authors and their affilations:
- Alexandra Belyaeva [@BelyaevaAlex](https://github.com/BelyaevaAlex), MSU +BI
- Alexandr Zhuravlev [@Baragozin](https://github.com/Baragozin), PNPI + BI

Scientific supervisor and her affilation:
- Polina Drozdova, ISU · Institute of Biology

## Introduction
------------
In phylogenetic studies, the problem of reproducibility of data is a major issue. This problem arises due to several factors, including the complexity of the data, the use of different methods and software packages, and the lack of standardization in reporting results. The use of different methods and software packages can lead to different results, which makes it difficult to compare and reproduce results. Additionally, the lack of standardization in reporting results makes it difficult for other researchers to reproduce the same results. 

This project examined the issue of phylogenetic analysis reproducibility in a number of studies on freshwater crustaceans based on open data. 
## Data
------------
__1. Moskalenko: Specimens were collected by the authors of the article in 2012 and 2013 field seasons on Lake Baikal  and preserved in 90% ethanol. A subset of species for sequencing was chosen to represent major lineages within Baikal amphipods, including lineages close to basal nodes. Non-Baikal freshwater palearctic Gammarus species morphologically similar to Baikal amphipods were also considered. Sequences that branched out beyond Gammarus species, namely Hyalella azteca and Dikerogammarus villosus were used 
as outgroups.__

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

__2. Bystřický: Individual lines of G. fossarum in their contact zone in the Western Carpathians are being studied.__

There is no source data. But there is an archive already with ready alignments in supporting information: https://aslopubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Flno.12239&file=lno12239-sup-0007-DataS1.zip 

__3. Bukin: Specimens of the amphipod G. fasciatus were collected in the littoral zone of Lake Baikal + previously published in pop-set data.__

URL-link for pop-set data:
https://www.ncbi.nlm.nih.gov/popset/376403838

The sequences obtained in the study were deposited in the database under the following accession numbers: MG214957, MG214958.

MG214957:
https://www.ncbi.nlm.nih.gov/nuccore/MG214957.1/

MG214958:
https://www.ncbi.nlm.nih.gov/nuccore/MG214958.1/

__4. Wattier: COI sequences for the DNA-barcode region were either newly generated or derived from the literature for 4926 Gammarus fossarum individuals, from 498 sampling sites distributed throughout 19 European countries.__

All these sequences were derived into three groups.

1. All newly produced COI sequences are released in GenBank under accession numbers: MT978656–MT980726.
https://www.ncbi.nlm.nih.gov/nuccore/?term=MT978656:MT980726[pacc]

2.The 463 sequences from Lagrue et al. which were initially unreleased in Genbank are now available under accession numbers: MT411018-MT411480.
https://www.ncbi.nlm.nih.gov/nuccore/?term=MT411018:MT411480[pacc]

3.Sequences and metadata are available in BOLD dataset DS-GFOSCDEU.
http://www.boldsystems.org/index.php/Public_SearchTerms?query=DS-GFOSCDEU

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
“Must have” list with essential for replication checkpoints:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/check_list_3.png)

Our list of shortcomings, typos, inconsistencies in the articles reviewed:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/check_list_2.png)

“Must have” list with essential for this kind of phylogenetic analysis:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/check_list_1.png)

Also to process the data from both the source and the tools obtained during the application, 3 small scripts were written to generate the correct input data for the next tools in the queue. All of them are presented in the file Moskalenko_1.ipynb. For convenience, they are designed in the form of functions.

According to one of the articles, a comparative analysis of the use of tools at different stages of our analysis was carried out. 

6 tools for multiple alignment were used, as a result, the tool (MUSCLE) used by the authors of the article did not show particularly great benefits in its use:

![](https://github.com/BelyaevaAlex/-Crustacean-hoods-strive-for-power/blob/main/muscle_compare.png)

Our analysis showed that kalign tool wins on all indicators based on our data.

Further, at the stage of data concatenation, it turned out that a tool like catfasta2phyml is great for getting the output of a data splitting scheme in a file.cfg for PartitionFinder, however, concatenated our data incorrectly, which did not even allow us to trim sequences using the trimAl tool. The MEGA 11 program was chosen for concatenation: https://www.megasoftware.net/web_help_10/Concatenation_Utility.htm

Although this was not stated in the article, but by the look of our alignments, it was decided not to skip the trimming stage, as mentioned earlier.

## Summary
------------

Two of the four articles should be considered complex and poorly reproducible for the following reasons: 1) the absence of source data in some of them; 2) information about the source data and explanations to them are located throughout the article; 3) the tools are indicated without their versions and provided parameters; 4) the authors miss the moments of correcting output and, accordingly, input data during pipeline execution. The last two items of this list of shortcomings also apply to the other two articles.

## References
------------
- [Moskalenko, V. N.; Neretina, T. V.; Yampolsky, L. Y. To the origin of Lake Baikal endemic gammarid radiations, with description of two new Eulimnogammarus spp. Zootaxa 4766(3), 457-471 (2020). https://doi.org/10.11646/ZOOTAXA.4766.3.5](https://www.mapress.com/zt/article/view/zootaxa.4766.3.5)
- [Bukin, Y.S., Petunina, J.V. & Sherbakov, D.Y. The Mechanisms for Genetic Diversity of Baikal Endemic Amphipod Gmelinoides fasciatus: Relationships between the Population Processes and Paleoclimatic History of the Lake. Russ J Genet 54, 1059–1068 (2018). https://doi.org/10.1134/S1022795418090053](https://link.springer.com/article/10.1134/S1022795418090053#citeas)
- [Bystřický, P.K., Rutová, T., Brož, V., Gajdošová, M., Juračka, P.J., Copilaş-Ciocianu, D. and Petrusek, A. Distribution patterns at different spatial scales reveal reproductive isolation and frequent syntopy among divergent lineages of an amphipod species complex in Western Carpathian streams. Limnol Oceanogr 67, 2796-2808 (2022). https://doi.org/10.1002/lno.12239](https://aslopubs.onlinelibrary.wiley.com/doi/10.1002/lno.12239)
- [Wattier, R., Mamos, T., Copilaş-Ciocianu, D. et al. Continental-scale patterns of hyper-cryptic diversity within the freshwater model taxon Gammarus fossarum (Crustacea, Amphipoda). Sci Rep 10, 16536 (2020). https://doi.org/10.1038/s41598-020-73739-0](https://www.nature.com/articles/s41598-020-73739-0)



