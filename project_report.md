# What causes antibiotic resistance in *E.coli*?
Project report

Authours: Anna Kalygina and Alexandr Lavrov

## Abstract
 The elucidation of the molecular details of antibiotic resistance will lead to improvements inextending the efficacy of current antimicrobials. For this purpose, the combined bioinformatics approach was used to locate mutations associated with antimicrobial resistance in the chromosomal sequences of *Escherichia coli* K-12 strain. We identified single nucleotide variants associated with the ampicillin resistance using a combination of VarScan and snpEff.

**Keywords:** ampicillin resistance in *E.coli*, bioinformatic analysis.

## Introduction
Treating bacterial infections is increasingly complicated by the ability of bacteria to develop resistance to different antibiotics. Resistance to antibiotics can be caused by a variety of mechanisms: (i) the presence of an enzyme that inactivates the antimicrobial agent; (ii) a mutation in the target of the antimicrobial agent that reduces its binding capacity; (iii) post-transcriptional and post-translational modification of the target of the antimicrobial agent, which reduces its binding capacity; (iv) reduced uptake of the antimicrobial agent; and (v) active efflux of the antimicrobial agent.

In our study, we analyse full genome sequence of an ampicillin resistant *E.coli* strain to identify single nucleotide mutations responsible for acquiry of antibiotic resistance. Description of mutations improves our understanding of universal pathways thatform barriers for antimicrobial agents and helps to identify the possible targets for bacterial treatment.

## Methods
#### Alignment
 The reference genome sequence for *Escherichia coli* K-12 MG1655 strain was obtained from NCBI with accession number NC 000913.3. The sequence data from ampiccilin resistant *E.coli*, produced by Illumina paired-end squencing, obtained from the open database.

 The reads were processed using Trimmomatic (Bolger, Lohse, & Usadel ,2014).
 The pipeline uses BWA-MEM algorithm (Li and Durbin, 2010) to map each read end separately to the  K-12 MG165 reference genome.

#### Variant calling


## Results
resistant to ampicillin (Ampr10) carry mutations in a locus which we have designated ampA

## Discussion

## References



