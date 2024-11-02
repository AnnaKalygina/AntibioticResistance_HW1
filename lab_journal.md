# Lab journal for the project on antibiotic resistance

### 1. Load necessary files and insepct raw sequences
The [reference genome](/GCF_000005845.2_ASM584v2_genomic.fna.gz) of K-12 *E.coli* strain and [annotation file](/GCF_000005845.2_ASM584v2_genomic.gff.gz) were downloaded from NCBI FTP with the following accession number: GCF_000005845.2_ASM584v2.

The reads from pair-end Illumina shortgun sequencing ([R1](/amp_res_1.fastq.gz) and [R2](/amp_res_2.fastq.gz))of an ampicillin-resistant strain were downloaded from the [following database](https://figshare.com/articles/dataset/amp_res_2_fastq_zip/10006541/3).

To inspect the raw sequences without unpacking the zipped file we use command:
``` bash
gzcat amp_res_1.fastq.gz | head
```

Each read has 4 lines of information. The first line starts with the @ symbol, and contains identifiers and information about this read. The second line contains the actual sequence, then on line three there is a ‘+’, which may sometimes have the identifier and info repeated. Line 4 contains the quality string, where ASCII characters encode the quality score for each base. The quality score ranges from 0 to about 40; the higher the number, the greater the accuracy of the base call. 

To check how many reads we have we run: 
``` bash
gzcat amp_res_1.fastq.gz| grep -v ">" | wc -l
```
There are 1467964 reads for the R1 and 1497339 reads for the R2.

file                format  type  num_seqs     sum_len  min_len  avg_len  max_len
amp_res_1.fastq.gz  FASTQ   DNA    455,876  46,043,476      101      101      101

file                format  type  num_seqs     sum_len  min_len  avg_len  max_len
amp_res_2.fastq.gz  FASTQ   DNA    455,876  46,043,476      101      101      101

Next, we inspect files using FastQC.
According to the FASTQC the per base sequence quality was abnormal (red). Per tile sequence quality, per base equence content and per base GC content were unusual (orange).

[Per base sequence quality](/img/per_base_sq.png)
[Per base sequence content](/img/per_base_sequence_content.png)


### 2. Filtering reads
For reads filtering the Trommomatic package was used. We run Trimmomatic in paired end mode, with following parameters:
- Cut bases off the start of a read if quality below 20
- Cut bases off the end of a read if quality below 20
- Trim reads using a sliding window approach, with window size 10 and average quality  within the window 20. 
- Drop the read if it is below length 20.


``` bash
java -jar trimmomatic-0.35.jar PE -phred33 -trimlog trim.log amp_res_1.fastq.gz amp_res_2.fastq.gz paired_output_amp_res_1.fastq.gz unpaired_output_amp_res_1.fastqc.gz paired_output_amp_res_2.fastq.gz unpaired_output_amp_res_2.fastqc.gz SLIDINGWINDOW:10:20 LEADING:20 TRAILING:20 MINLEN:20
```
We got the following output: 
```
Input Read Pairs: 455876 Both Surviving: 445689 (97.77%) Forward Only Surviving: 9758 (2.14%) Reverse Only Surviving: 284 (0.06%) Dropped: 145 (0.03%)
TrimmomaticPE: Completed successfully
```

Now we repeat the fastqc analysis on trimmed data.


### 3. Aligning sequences to reference
#### BWA MEM alignment 
Here is the alignment statistics: 
891635 + 0 in total (QC-passed reads + QC-failed reads)
891378 + 0 primary
0 + 0 secondary
257 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
890569 + 0 mapped (99.88% : N/A)
890312 + 0 primary mapped (99.88% : N/A)
891378 + 0 paired in sequencing
445689 + 0 read1
445689 + 0 read2
887530 + 0 properly paired (99.57% : N/A)
889384 + 0 with itself and mate mapped
928 + 0 singletons (0.10% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)


#### Compressing SAM file to BAM
The alignment files were then compressed, sorted and indexed: 
``` bash 
samtools view -S -b alignment.sam > alignment.bam
samtools sort alignment.bam alignment_sorted.bam
samtools index alignment_sorted.bam
```
#### Visualizing BAM file
We use Desktop IGV for alignment visualisation. Additionally, the index file for the reference fasta was generated: 
``` bash 
samtools faidx GCF_000005845.2_ASM584v2_genomic.fna
```
### Variant calling 
First, we pileup the bases in the reads that do not match reference genome: 

``` bash
samtools mpileup -f GCF_000005845.2_ASM584v2_genomic.fna alignment_sorted.bam > my.mpileup
```
Second, we use VarScan to scan for variants in the reads. 
The threshold for the minimal variation frequency is 0.30. 6 SNPs were identified (The same result was for --min-var-freq set at 0.50). It allows for the detection of both high-frequency mutations and significant subpopulations, while reducing the likelihood of errors being called as true variants.

``` bash
varscan mpileup2snp my.mpileup --min-var-freq 0.30 --variants --output-vcf 1 > VarScan_results_var_freq_0.3.vcf
```

Only SNPs will be reported
Warning: No p-value threshold provided, so p-values will not be calculated
Min coverage:   8
Min reads2:     2
Min var freq:   0.3
Min avg qual:   15
P-value thresh: 0.01
Reading input from my.mpileup
4641476 bases in pileup file
9 variant positions (6 SNP, 3 indel)
0 were failed by the strand-filter
6 variant positions reported (6 SNP, 0 indel)

The following positions were recognised as stable mutations:
- 93,043 C -> G (in the gene-b0085, murE)
- 482,698 T -> A (gene-b0462, acrB)
- 852,762 A -> G (gene-b4416, small RNA RybA)
- 1,905,761 G -> A (gene-b1821, mntP)
- 3,535,147 A -> C (gene-b3404, envZ)
- 4,390,754 G -> T (gene-b4161, rsgA)