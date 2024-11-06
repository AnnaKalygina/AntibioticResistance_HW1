# A bioinformatic approach to detecting mutations conferring ampicillin resistance in *Escherichia coli*

**Author**: Anna Kalygina  
**Affiliation**: Bioinformatics Institute, Saint Petersburg, Russia

## Abstract
Antibiotic resistance in bacteria, such as *Escherichia coli*, complicates the treatment of infections, often due to mutations in target proteins and active efflux mechanisms. In this study, we analyzed the genome of an ampicillin-resistant *E. coli* strain to identify key mutations involved in resistance. We found several mutations affecting cell wall synthesis and efflux pumps, which may contribute to the observed resistance, providing insights into potential targets for future treatments.

**Keywords**: ampicillin resistance in *E. coli*, bioinformatic analysis

## Introduction
Treating bacterial infections is increasingly complicated by the ability of bacteria to develop resistance to different antibiotics. Resistance mechanisms include mutations in the target of the antimicrobial agent, reducing its binding capacity, and active efflux systems that remove the agent from the cell [@erb2007prevalence]. These modifications are particularly effective against antibiotics like ampicillin, which rely on intracellular accumulation for their efficacy [@doi:10.2147/IDR.S221212].

Ampicillin is a $\beta$-lactam antibiotic commonly used to treat a range of bacterial infections by targeting cell wall synthesis [@rafailidis2007]. However, resistance to ampicillin in bacteria such as *Escherichia coli* is on the rise, often due to mutations affecting target proteins or efflux systems.

In this study, we analyzed the full genome sequence of an ampicillin-resistant *E. coli* strain to identify single nucleotide mutations contributing to resistance. We identified several mutations in genes related to efflux pumps and cell wall synthesis, which may play a crucial role in conferring resistance. Understanding these mutations provides insights into the mechanisms of antibiotic resistance and potential targets for future treatment strategies.

## Methods

### Data Preprocessing
The reference genome sequence for *Escherichia coli* K-12 MG1655 strain was obtained from NCBI with accession number NC 000913.3. The sequence data from ampicillin-resistant *E. coli*, produced by Illumina paired-end sequencing with a read length of 150 bp, was obtained from the open database [@figshare2019].

The reads were inspected using FastQC [@andrews2010] (version 0.12.0) to evaluate quality metrics, such as per-base quality scores, GC content, and adapter content. The reads were then filtered using Trimmomatic [@bolger2014] (version 0.40) with the following parameters:

| **Parameter**         | **Description**                                                            |
|-----------------------|----------------------------------------------------------------------------|
| `SLIDINGWINDOW:10:20` | Trim using a sliding window of size 10 and average quality threshold of 20 |
| `LEADING:20`          | Remove low-quality bases from the start if quality score is below 20       |
| `TRAILING:20`         | Remove low-quality bases from the end if quality score is below 20         |
| `MINLEN:20`           | Discard reads shorter than 20 bp after trimming                            |

**Table 1. Summary of Trimmomatic trimming parameters.**

### Read Alignment
The pipeline used the `BWA-MEM` algorithm [@li2010] (version 0.7.17) to map each read to the K-12 MG1655 reference genome. BWA-MEM was selected because it is well-suited for handling reads of 70 bp or longer. SAM files were converted to BAM format and sorted using `samtools` (version 1.13). Sorted BAM files were indexed for efficient querying of reads by genomic position. For resulting read statistics, refer to the Supplementary Information (Table S1).

### Variant Calling
Variant calling was performed on the aligned sequences. Reads were first piled up using `samtools mpileup` [@10.1093/gigascience/giab008].

Variants were called using VarScan [@koboldt2012] (version 2.4.0). A minimum variant frequency of 30% (`min-var-freq = 0.3`) was used. This threshold was chosen to prioritize identifying clonal or near-clonal variants.

Single nucleotide polymorphisms (SNPs) and small insertions and deletions (indels) were automatically annotated using SnpEff [@cingolani2012program] (version 5.2c) with the K-12 MG1655 reference genome database.

### Data and Software Availability
All data processing steps were performed in a conda environment to ensure reproducibility. The environment configuration and the analysis scripts are available in the project's [GitHub repository](https://github.com/AnnaKalygina/AntibioticResistance_HW1).

## Results

### Several Polymorphisms with Moderate and High Impact Were Identified
Following the alignment of the sequencing reads and variant calling, a total of six significant variants were identified in the *E. coli* K-12 MG1655 strain (Table 2). These variants were characterized by their chromosomal location, impact on gene function, and predicted protein changes.

The majority of the variants (four out of six) were classified as having a moderate impact on gene function, indicating potential changes in the protein product that could affect biological activity. Three of these variants were missense mutations in the *ftsI*, *acrB*, and *mntP* genes, resulting in amino acid substitutions at p.Ala544Gly, p.Gln569Leu, and p.Gly25Asp, respectively. These mutations may alter the functional domains of the encoded proteins, potentially affecting their activity.

One variant, classified as a modifier, was identified upstream of the *glnH* gene. This upstream gene variant might influence the regulation of the gene, although its effect on protein production is not direct and likely to be of lower impact.

Additionally, a synonymous variant was found in the *rsgA* gene at position c.756C>A (p.Ala252Ala), which was categorized as having a low impact. Synonymous variants do not result in an amino acid change and are therefore not expected to significantly impact protein function.

The depth of coverage for each identified variant ranged between 13 and 17, indicating the reliability of results. Further functional analysis may be necessary to determine the precise phenotypic consequences of the identified missense and modifying variants.

| **Gene Name** | **Position** | **Ref. all** | **Alt. all.** | **Impact** | **Variant Type**  | **Amino Acid Change** |
|---------------|--------------|--------------|---------------|------------|-------------------|-----------------------|
| *ftsI*        | 93043        | C            | G             | MODERATE   | missense          | p.Ala544Gly           |
| *acrB*        | 482698       | T            | A             | MODERATE   | missense          | p.Gln569Leu           |
| *glnH*        | 852762       | A            | G             | MODIFIER   | upstream gene     | -                     |
| *mntP*        | 1905761      | G            | A             | MODERATE   | missense          | p.Gly25Asp            |
| *envZ*        | 3535147      | A            | C             | MODERATE   | missense          | p.Val241Gly           |
| *rsgA*        | 4390754      | G            | T             | LOW        | synonymous        | p.Ala252Ala           |

**Table 2. Variants identified in the *E. coli* K-12 MG1655 strain, detailing the reference and alternative alleles, gene impact, and protein change.**

## Discussion
While some mutations, such as those in genes like *envZ*, may be neutral or have limited influence on antibiotic resistance, several of the identified variants could be strongly associated with antibiotic-resistant phenotypes in *E. coli*. We suggest that missense mutations in *ftsI*, *acrB*, and *mntP* are strong candidates contributing to ampicillin resistance.

The *ftsI* gene encodes penicillin-binding protein 3 (PBP3), which plays a crucial role in synthesizing the peptidoglycan septum during cell division. PBP3 has a high affinity for $\beta$-lactam antibiotics, including penicillins, and mutations in *ftsI* could reduce binding affinity, thereby diminishing the efficacy of ampicillin [@GAbotta]. Specifically, the identified p.Ala544Gly substitution may alter the structure of the active site, affecting PBP3's interaction with ampicillin.

Both *mntP* and *acrB* encode efflux pumps. *MntP* functions as a manganese exporter, while *AcrB* is part of the AcrAB-TolC efflux system, known for its broad substrate specificity, including various antibiotics. The p.Gln569Leu mutation in *acrB* and p.Gly25Asp mutation in *mntP* could alter the conformation of these pumps, potentially increasing the efficiency of ampicillin efflux. Mutations in efflux systems like AcrB have been previously associated with multidrug resistance in *E. coli* [@sennhauser2007, @marciano2022].

The presence of these missense mutations suggests a possible multitarget mechanism of resistance, where modifications in cell wall synthesis and efflux contribute synergistically to ampicillin resistance. The synonymous variant identified in *rsgA* is likely to have minimal effect on resistance due to the lack of amino acid change, indicating that it may not play a significant role in the resistant phenotype.

Our findings indicate that treatment in similar resistance profiles should focus on alternative antimicrobials that target distinct pathways, such as those not affected by efflux mechanisms or $\beta$-lactam target modifications. Experimental validation through biochemical assays and phenotypic susceptibility tests is required to elucidate the specific roles of these mutations in conferring antibiotic resistance.

## Supplementary Information

### Read Processing Results
The data was preprocessed and aligned with 97.34% reads surviving.

| **Step**             | **Forward reads count** | **Reverse reads count** | **% surviving of original** |
|----------------------|-------------------------|-------------------------|-----------------------------|
| Before processing    | 455,876                 | 455,876                 | 100                         |
| Trimming             | 445,689                 | 445,689                 | 97.77                       |
| Alignment            | 443,765                 | 443,765                 | 97.34                       |

**Table S1. Read count across analysis steps.**

---

### Bibliography
The bibliography can be included using a reference section with BibTeX or manual entries, as used in LaTeX.
