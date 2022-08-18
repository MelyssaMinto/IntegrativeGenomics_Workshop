# Intro to Epigenomics (Day 1)



### Learning Outcomes 
- You will learn several types of epigenomic regulation and data
- You will understand how epigenomic sequencing data is processed for analysis 
- Practice generating reading in count data into R and performing differential peak analysis 

## The central dogma of biology
![Adapted from: Chromatin modification and the epigenetic landscape in a materials context. Namec et al, 2021](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41578-020-00238-z/MediaObjects/41578_2020_238_Fig1_HTML.png?as=webp)
> ℹ️ Adapted from Namec et al (2021) Nat Rev

The **central dogma of molecular biology** has long reigned as the transformation of DNA to RNA (transcription) and from RNA to Protien (translation). We now know that proteins can interact with DNA to modulate the production of RNAs via trans-regulation. Simimarly, modifcations to the DNA (i.e. DNA methylation, histone tail modifications) work via cis-regulation to modulate the production of RNAs. Further we will discuss different levels of epigentic control of the transcription. 

## What is Epigenomics?
**Epigenomics** is the study of molecular elements that regulate gene expression. The word Epigenome roughly translates to "on top of the genome" with "Epi" meaning top. Put simply, its the study of what is on our DNA. Historically, this has been studied by examining various modes of genomic data such as chromatin accessibility, Transcription Factor (TF) binding, mdifications to the DNA (i.e. DNA methylation, histone tail modifications), and the overla 3D conformation of the genome. Each of these layers have some functional role underlying the ability for RNAPolII to bind to promoters and initiate transcription. 



### Chromatin Accesibilty
![Adapted from: A continuum of accessibility states broadly reflects the distribution of chromatin dynamics across the genome, Klemm et al, 2019 ](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41576-018-0089-8/MediaObjects/41576_2018_89_Fig1_HTML.png?as=webp)
> ℹ️ Adapted from Klemm et al (2019) Nat Rev

The DNA is wrapped around **nucleosomes** which then can wrap around itself like a tangled telephone cord which leads to a compact chromatin state. Chromatin compaction is throught to be one of the largest driving forces in transcriptip. Transcriptional machinery and other epigenetic factors simply cannot bind the the DNA if the chromatin is "closed" or in a compact **heterochromatin** state. 


![Adapted from: Principal methods for measuring chromatin accessibility., Klemm et al, 2019 ](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41576-018-0089-8/MediaObjects/41576_2018_89_Fig2_HTML.png?as=webp)
> ℹ️ Adapted from Klemm et al (2019) Nat Rev

Recent and common methods to assay chromatin accessibility are ATAC-seq and DNase-seq.  Both methods are able to pull down regions of the genome that are more "open" and the sequencing readouts shows pileups or peaks of areas where the chromatin is accessible. 

![Adapted from: Population-scale measurements of chromatin accessibility reflect the average accessibility of a heterogeneous collection of single molecules., Klemm et al, 2019 ](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41576-018-0089-8/MediaObjects/41576_2018_89_Fig3_HTML.png?as=webp)
> ℹ️ Adapted from Klemm et al (2019) Nat Rev


| | |
|-|-|
|`NOTE` | Methods vary in the way they capture chromatin and careful consideration should be made in determining the method of choice and tools to best analyize the data. |
---

### 3D Functional Domains: TADs and LADs
Even though the supercoiling of a tangled telephone coord might seem random, the compation and 3D organization of the DNA is not. 

![Adapted from: Hierarchical organization of chromatin structure., Bonev & Cavalli, 2016 ](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnrg.2016.112/MediaObjects/41576_2016_Article_BFnrg2016112_Fig2_HTML.jpg?as=webp)
> ℹ️ Adapted from Bonev & Cavalli (2016) Nat Rev


### Transcription Factors & Chromatin Remodelors 
Proteins can interact with the DNA and have epigenomic control. **Transcription factors** (TFs) are a class of DNA binding proteins that bind ina sequence specific way. TFs are typically expressed in a precise spatiotemporal pattern when they are triggered by signalling cascades from biologucal cues like developmetal triggers or a particulat stimulus. TFs are known to bind to regulatory regions such as enhancers and promoters to activate and silence programs of gene expression.

**Chromatin remodelers** are proteins that interact with the DNA in an ATP-dependent manner to alter nucleosome positioning and thus the 3D archetecture of the genome. These remodelers are known to form complexes that mediate histone turnover by inducing nycleosome sliding. This shifting of nucleomes thus reveals sites for TF binding.

### DNA Methylation
Methyl groups can be added to bases of DNA by **Dnmt** proteins while **Tet** proteins are known to remvove DNA methylation. DNA methylation is associated with gene silencing. DNA methylation occurs widely at **CpG sites**, that is where Cytonsine and Guananine follow each other in a genomic sequence. However, **CpG Islands** are long stretches of DNA that is CG rich with less nucleosome desity that is not typically heavily methylated. More recent research has shown that non CpG methylation is more common in neuronal tissue (specifically at CpA sites) and plays an important role in gene regulation.

![](https://www.mdpi.com/genes/genes-08-00148/article_deploy/html/images/genes-08-00148-g002-550.jpg)
> ℹ️ Adapted from Jang et al. (2017) Genes 

### Histone Tail Modification
As DNA is wrapped aound nuclesomes comprised of histone protiens, those histones have flexible tails that can be postranslationally modified which reuslts is chnages in chromatin conformation and accesbility. For instance, acetylation of Hisone 3 Lysine at position 27 (H3K17ac) is associated with eurchromatin, active enhancers and transcription and tri-methylation of Hisotone 3 at Lysine 27 is associated with herterochromain and gene silencing. One can assay for different histone tail modifcations via ChIP-sequencing. Similar to TF CHip, one can align these reads to the genome to get pileups of reads and infer where specific modifactions occur in the geneome. Some histone modications tend to have a large spread whule others are kinown to have discrete pins of binding resulting in sharp peaks. This should be considered when determing the aligortihms and parameters used to determine loci from sequecing assays


<details>
<summary style="font-weight: bold;">  ✎ What is the difference between DNA Methylation and Histone Methylation? ✎ </summary>
<br>
<blockquote>
DNA methylation occurs directly on the nucleotide and is typically associated with chromatin repression. Histone methylation occurs on the tails of histone. Specific residues of the histone tail can me methylated up to three times (mono-, di-, tri- methylation). Histone methylation is associated with varying chromatin and gene regulation depeding on the residue that is methylated and the how many methyl groups are added (1-3). 
</blockquote>
</details>




---

## Raw Reads to count matrix
There are an increasing number of technologies to assay the different epigenomic marks. All of these technologies produce segments of sequences that was captured. Generally, the goal is to use alogirithms to estimate where an epigenitic marks on the genome by taking the raw sequences, align them to a reference genome, and estimate a quantitative "count" of where those sequences tend to pileup. 

![](https://www.encodeproject.org/images/c45f4d8c-0340-4fcb-abe3-e4ff0bb919be/@@download/attachment/EncodeDatatypes2013-7.png)
> ℹ️ Adapted from ENCODE, Image credits: Darryl Leja (NHGRI), Ian Dunham (EBI), Michael Pazin (NHGRI)

### Raw Sequencing data 
![](https://ars.els-cdn.com/content/image/1-s2.0-S1046202320300591-gr1.jpg)
>  ℹ️ Adapted from Nakato and Sakata (2021) Methods
Raw sequencing data is typically stored in a .fastq file and holds the raw sequenceing reads as well as important information about each read including the quality of each base. Using that data, one can use an alignner (e.g. Bowtie2, STAR) to map the reads to a reference genome. Then one can use a peak caller (e.g. MAC2) to find where many reads pile up in the genome. Once you know where the peaks of reads are in the genome there are various down stream analyses:
1. Differential peak analysis
2. Annotation of the peaks by genomic loci or chromatin state
3. Mapping peaks to genes to estimate which genes those marks may regulate.

In this course we will demonstrate a selection of methods on how to do differential peak analysis. The major steps include 1) calling peaks 2) generating the peak count matrix 3) analyzing count data 4) clustering peaks by count 5) Differential analysis of the peaks between two conditions

### Calling peaks
To call peaks, we will use [MACS2](https://pypi.org/project/MACS2/) a widely used tool to call peaks for ChIP-seq, ATAC-seq, and DNase-seq data. There are paramters like  `--broad` that can be used to tune the peak calling for epigeneitcs marks that may either sharp narrow peaks or broad-spread peaks. 

Here is an example usage for calling ATAC peaks as suggested by MACS2 documentation
``` bash
macs2 callpeak 
```

The peaks will likely be in some sort of .bed format with at minimum three columns representing the Chr Start and End postions of each peak. 

```bash
chr1 
```

Once peaks are called from each sample, the next step is to compare the peaks between samples or conditions. To do that, a consensus peak set must me made. This can be done in various ways from simply merging all the peaks into one file or using tools to statisitvally infer confident peaks from each set. `bedtools merge` is a common tool to merge peak sets.

![](https://bedtools.readthedocs.io/en/latest/_images/merge-glyph.png)

To merge peaks sets you need to first combine them into one file, sort it and then merge. One important step in this is to filter out regions of the genome that are higly mappable or that have a lot of repeats. This region is called the 'blacklist'. Common genomes have blacklists that can be found here: https://github.com/Boyle-Lab/Blacklist/tree/master/lists

```bash
cat peak_set1.bed peak_set2.bed peak_set3.bed > all_peaks.cat
bedtools sort all_peaks.cat > all_peaks.sort
bedtools merge all_peaks.sort > all_peaks.merge
bedtools filter all_peaks.merge blacklist.bed >  all_peaks.bed
```


| | |
|-|-|
|`NOTE` | the file extensions (eg .cat, .sort) are not specific to any program. Here I use them to help denote what step in merging process each file is. All the files that are made during this brocess are technically still .bed files which are tab delimited files text files with specific guidlines about what goes in each column. See more about the different .bed files here: [Bed File Format  - Definition and supported options](http://useast.ensembl.org/info/website/upload/bed.html)  |
---
<details>
<summary style="font-weight: bold;">  ✎ What does the <code>cat</code> function do? ✎ </summary>
<br>
<blockquote>
<strong>cat</strong> is short for con<strong>cat</strong>enate which simply means to combine together in a series. For example,

```bash
cat "string1" "string2"
```
<pre>
string1string2
</pre>
</blockquote>


</details>


### Generating count matrix 
Once there is a consensus list of regions to compare, the next step is to calculate the number of reads (or counts) from these regions in each sample. This is usually organzed in a table where each column is a sample and each region (or peak) is a row. 

| peak | Sample1 | Sample2 |  Sample3 | Sample4 |
| --- | --- |--- | --- |--- |
| peak1 | 30 | 24 | 15 |10 |
| peak2 | 45 | 60 | 22 |18 |
| peak3 | 0 | 1 | 0 |0 |
| peak4 | 127 | 88 | 90 |98 |
| ... | ... | ... | ... | ... |

This can easly be done using the R package `Rsubreads`
```R
# -------load libraries--------
library(readr)
library(dplry)
library(Rsubreads)

#----Read in consensus peak set-----
peak_set <- read_tsv("all_peaks.bed", col_names = F)
bam_paths = list.files(".", pattern = ".bam")

#-----Get Peak Counts ---------------
# the function needs two inputs: the peakset 
# and the paths to the samples. 
# The peakset needs to have five columns 
# GeneID | Chr | Start | End | Strand 

# formatting the consensus peak set
peak_set = peak_set %>%
mutate(GeneID = paste0("peak", 1:n()),
       Chr = X1,
       Start = X2,
       End = X3,
       Strand = ".") %>%
select(GeneID, Chr, Start, End, Strand)

# run function to get counts
peak_counts = featureCounts(bams = bam_paths,
                            feature.file = peak_counts,
                            nThreads = 4,
                            isPaired = FALSE)

# check results 
peak_counts@counts
```

<details>
<summary style="font-weight: bold;">  ✎ How can you learn more about a specific R function? ✎ </summary>
<br>
<blockquote>

1. You can google! 🔍 someFunction() R
2. Use the R help function within the console: > ?sumFunction() 
3. Search for the function in the Packages pane in Rstudio 

<img src="https://confluence.wustl.edu/download/attachments/189860067/rstudio-packages.png?version=1&modificationDate=1639065252927&api=v2" alt="R Studio Packages Pan" height="300" >


Good functions tend to have deailed descriptions of thier arguments and some nice examples to model your code off of. It is generally good practice to read over the different aurguments ( especially the defaults) that each function has to ensure that it is optimized to your data. 

</blockquote>
</details>


### Analyzing count data
Now that there is a count matricx we can ask different questions with the data. 
1. How do the samples cluster?
2. Where are these peaks in the genome? 
3. Which Peaks are different from each other?
4. Are there any patterns in the peaks that are change that we can see?

#### Clustering
A good quality assurance check when working with count data, especially between different conditions is to check how well conditions cluster together. This is commonly done using a PCA. PCA is unsupervised modeling technique that decreases the number of dimensions in the data. The first principle component (PC1) is a line/plane in the data that explains of the highest variation in that data. The second Principle Component (PC2) is perpendicular to PC1 that explains the 2nd highest variation in the data. In sequencing analysis PCA is used to plot the top variation in the data to show that conditions cluster together. 

We can do this in R
```R
# -------load libraries--------
library(PCA)
pca()
```


#### Differential Analysis 
Ultimatrly, most researchers want to know what the gene expression or binding differences are between two or more groups. There are many tools that can carry out differential analysis whhich vary in statistical approaches. The typical pupeline includes 1) getting the counts 2) filtering and normalization 3) statistical testing 4) select significant features. A large part of what varies between tools is 1) the underlying assumptions about how the count data is distributed and 2) how they model feature-wise and sample-wise variance/dispersion. This workshop will focus on the use of `DESeq2` but there are other tools like `diffbind` and `` that are worth looking into. 

With `DESeq2` it is reccomended that you input raw unnormalized counts. The tool will compute a log-scaled normalization and filtering bsaed on the distrubition of the average count value of each feature. In our case our features are peaks and our columns are samples. 

```R
#------Load libraries
#------seet up seseq object 
#-------make design matrix
#-------perform analysis
#-------visualizations
```



## Glossary
| | |
|-|-|
|**Nucleosome** | unit of chromatin that is comprrised of 147 bps of DNA wrapped around an octomer of histone proteins (H2A, H2B, H3, and H4)|
|**Heterochromatin** | Methods vary in the way they capture chromatin and careful consideration should be made in determining the method of choice and tools to best analyize the data. |
|**Euchromatin** | Methods vary in the way they capture chromatin and careful consideration should be made in determining the method of choice and tools to best analyize the data. |
---

## References

- [Nemec, S., Kilian, K.A. Materials control of the epigenetics underlying cell plasticity. Nat Rev Mater 6, 69–83 (2021).](https://doi.org/10.1038/s41578-020-00238-z)
- [Klemm, S.L., Shipony, Z. & Greenleaf, W.J. Chromatin accessibility and the regulatory epigenome. Nat Rev Genet 20, 207–220 (2019). ](https://doi.org/10.1038/s41576-018-0089-8)
- [Bonev, B., Cavalli, G. Organization and function of the 3D genome. Nat Rev Genet 17, 661–678 (2016).](https://doi.org/10.1038/nrg.2016.112)
- [Moore, L., Le, T. & Fan, G. DNA Methylation and Its Basic Function. Neuropsychopharmacol 38, 23–38 (2013).](https://doi.org/10.1038/npp.2012.112)
- [Jang HS, Shin WJ, Lee JE, Do JT. CpG and Non-CpG Methylation in Epigenetic Gene Regulation and Brain Function. Genes. 2017; 8(6):148.](https://doi.org/10.3390/genes8060148)
- [West AE, Greenberg M. Neuronal; Activity-Related Gene Transcription in Synapse Development and Cognitive Function](https://cshperspectives.cshlp.org/content/3/6/a005744)
- [Nakato R, Sakata T. Methods for ChIP-seq analysis: A practical workflow and advanced applications. Methods. 2021 Mar;187:44-53. Epub 2020 Mar 30.](https://www.sciencedirect.com/science/article/pii/S1046202320300591)

