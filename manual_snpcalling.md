
# NGS for genetic variation

---

## Overview

* Goal
  * Understand the workflow
  * Pitfalls
  * Quality Control
  * Data formats

---

* Steps:
  * Create raw sequencing data
  * Quality Control
  * Mapping to a reference genome
  * Identify variants
  * Annotate variants

---

![Many sequence fragments](attachments/overview1.png)

---

![Map to a reference genome](attachments/overview2.png)

---

![Identify differences to the reference genome](attachments/overview3.png)

---

## Quality Control Everywhere !


* Garbage in - Garbage out
* Low tech
  * Files exist
  * Files are of realistic size
  * Files are not truncated
  * Files 'look' ok
  * Use checksums if possible
* High tech - depends on file format

* Important in automation!

---

## Preparation

* Get a sample from a patient, extract DNA
* Create a sequencing library (could be exome)
* Send to a sequencing provider
* **Beware**: paired end sequencing

![Library prep|400](attachments/libraryprep.png){=40%}

---

## Illumina sequencing


![](attachments/illumina_throughput.png)

---

# Raw FASTQ data

## Sequencing data in FASTQ format

~~~
@1121358 read ID
TGAATCTGGGAGGCGGAGGTTGCAGTGAGAGTGAGGCGAGATC sequence
+                           
<D<?><DIBB??BF;@;AB@CHG<A?F<@>@;F>>6?8??@@B sequence quality, every
@1121358
AGCCTCAAACTCCTGGGCTCAAGGGATCCTCACTTCTTGACCT
+
A:EFFACACFE?FCG5CG9>C@@B@<@=@FACFEE<F9C>E6F
@FCB06B3ABXX:8:1105:13917:3018
TCCCTTGAGCCCAGGAGTTTGAGGCTGCATTGAGCTATGATCA
+
DDHGHFJBCIIIEFECFCEFJCDGHHKHE?FECFHHBAJ?DCF
~~~

---

## Sequencing data in FASTQ format

Every four lines are a record

~~~
@1121358      <-- Read ID
TGAATCTGGGAG  <-- DNA sequence
+             <-- Separator (could contain the ID again)
<D<?><DIBB??  <-- Sequencing quality
~~~

---
## Quality score in PHRED format

~~~text
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
 |                         |    |        |
33                        59   64        73
 |                                       |
 0.2......................26...31........41

Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
~~~

see also: [Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format)

---

## Quality score in PHRED format

![Phred encoding](attachments/phred.png)

---

## Quality control of raw FASTQ

* Garbage in - Garbage out
* Low tech - check if fastq files:
  * Exist
  * Realistic size & number of lines/records
  * Look like fastq
  * Are not truncated
  * Checksums
* High tech - [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

---

## FASTQC

![Part of a good FastQC report](attachments/fastqc.png)

---

![Part of a bad FastQC report - we will see more later](attachments/fastqc_bad.png)

---

## What to do with bad FASTQC?


* Impact downstream

  * Fewer reads map
  * Reads map with lower quality
  * Less high quality SNPs called

* Be careful with batch effects
* Distinguish between:
  * An uncalled SNP
  * A locus that does not differ from the reference

---

* Fix quality:
  * Trim reads - remove bad quality
  * FastX-toolkit, Cutadapt, Flexbar, etc

* Practically - not required
  * Mappers deal with this
  * "Soft clipping"

---

## Mapping reads

* Using a reference genome
  * The reference has to be available

* Other methods
  * Assembly based
  * Alignment free methods

---

## Mapping reads - to a reference genome

![DOI:10.1093/bfgp/elu042](attachments/refalignment.png){width=80%}

---

## Human Reference Genome

![Welcome Collection London|300](attachments/welcome.png)

---

## Human Reference Genome

* Genome is complex
  * Large - 3,000,000,000 bases
  * Many repetitive regions
  * Gaps - the genome is not finished
  * **Multiple versions!**
    * Human - GRCh37/hg19
    * Human - GRCh38

---
## Mapping to a Human Reference Genome

* Lots of data
  * NovaSeq 6000 system outputs 6 Tb of data <2 days (~48 whole genomes)
  * Quality is good, however with 1 in 1000 mistakes - still many errors

---

## How does an aligner work?

* Build an reference index of k-mers
  * link a k-mer to a genome position


* Split the read in k-mers
  * Look each k-mer up in the index
  * Find genome positions where multiple k-mers point to
  * Fully build the alignment (slow)

For more information, read [Trapnell et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2836519/)

---

## How does an aligner work?

![Aligner|400](attachments/hashes2.png)

---


## Software

* Many tools available ([See here](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software#Short-Read_Sequence_Alignment))

* Procedure is often the same
  * Build an index
  * Map reads

We will use [BWA](http://bio-bwa.sourceforge.net/)

---
## Aligerns output SAM/BAM format

* Column based format
* Each line describes *one* read alignment to the reference genome

![](attachments/samformat1.png)

![](attachments/samformat2.png)

---

## SAM columns

![Sam table|400](attachments/sam_table.png)

[http://bio-bwa.sourceforge.net/bwa.shtml](http://bio-bwa.sourceforge.net/bwa.shtml)

---

## SAM flags

SAM Column 2 

e.g `16` -  `0b 0000 0001 000` -  `0x0010` => `strand of the query (1 for reverse)`


![flags|400](attachments/flags.png)

[http://bio-bwa.sourceforge.net/bwa.shtml](http://bio-bwa.sourceforge.net/bwa.shtml)

---

## SAM -> BAM

### SAM - Sequence Alignment/Map format
  * Uncompressed - Large
  * Slow lookup

### BAM - Binary Alignment/Map format
  * Compressed - Smaller
  * Can be index - fast lookup

Manipulate using [`samtools`](https://www.htslib.org/doc/samtools.html)

- Inspect & convert: `samtools view`
- Sort: `samtools sort`
- Create an index: `samtools index`

---

## Quality Control

* Low tech - size, basic formatting
* High tech
  * `samtools faidx` - no reads mapped to each chromosome
  * `samtools flagstat` - counts the number of alignments for each FLAG type

---

# Variant Calling

## What does a variant look like?

![SNPs in an alignment|800](attachments/alignment.png)

---

## Why again do we want to annotate variants?

* Diagnosis
* Prognosis
* GWAS, population genetics, etc

---

## Types of variants

* Point mutation
  * SNP, Single Nucleotide Polymorphism
    (Population Frequence > 1%)
* Insertion/Deletion
  * Small: indel
  * Large: genomic rearrangements
  * Repeat alterations
* Large rearrangments
  * Copy number variants
  * Translocations

---
## What does a variant look like?

![Not all mismatches are variants](attachments/alignment.png)

---

## Possible sources of problems

* Library preparion
* Sequencing error
* Mapping error & problems with the reference genome
  * Gaps and errors in the reference
  * Low complexity & repetitive regions

---

## How to call a variant?

Which tool to choose?

* Quality (e.g. false detection rates)
* Speed & memory usage
* Handle diploid vs polyploid vs complex mixtures?
* Sequencing platform - different platforms have different error models
* Different methods - heuristic - Bayesian - Deep learning

---

## Variant calling tools

Many different algorithms!

* [DeepVariant](https://github.com/google/deepvariant), uses deep learning
* [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)
* [Freebayes](https://github.com/freebayes/freebayes)
* SNVer, Varscan, Platypus...

We will use: [Bcftools](https://samtools.github.io/bcftools/bcftools.html)

---

## Variant Call Format (VCF)

Generic format for variants

* Can be indexed & compressed (using `tabix`)
* Can contain much additional information:
  * Indels, SNPs
  * Multi-sample
  * Multi-allele
  * genotype calls & statistics
  * Sequencing depth

* Format description: [https://samtools.github.io/hts-specs/VCFv4.2.pdf](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
---

## Variant call format (VCF)

1. **#CHROM**: Chromosome of the variant call
2. **POS**: Nucleotide position
3. **ID**: SNP ID (dbSNP) - if not (yet) known: a dot.
4. **REF**: Reference allele - sequence in the reference genome.
5. **ALT**: Alternative allele(s) - Variant sequence.
6. **QUAL**: Quality score - phred scaled (*e.g.* 20 -> 1% chance of error)
7. **FILTER**: Is this SNP filtered out? (needs manual posthoc filtering - `.` means no filter)
8. **INFO**: Extra information on the SNP call in question
9. **FORMAT ...**: Information on the SNP calls.

---

## Variant call format (VCF)

![VCF example](attachments/vcf.png)

---

## VCF - INFO fields

~~~
##INFO=<ID=INDEL,
        Number=0,
        Type=Flag,
        Description="Indicates that the variant is an INDEL.">
##INFO=<ID=DP,
        Number=1,
        Type=Integer,
        Description="Raw read depth">

---

DP=100;VDB=0.466639;SGB=-1.38624;MQSBZ=0;FS=0;MQ0F=0;AC=4;AN=4;DP4=0,0,3,97;MQ=60

~~~

---

## VCF - FORMAT fields

~~~
##FORMAT=<ID=PL,
          Number=G,
          Type=Integer,
          Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,
          Number=1,
          Type=String,
          Description="Genotype">

---

FORMAT	  TLE66_N.bam	    TLE66_T.bam
GT:PL	  1/1:255,178,0	    1/1:253,123,0

~~~


---

## Variant call quality control

**Low tech**

  * File based
  * Checksum?
  * Grep on the quality column..
  * Find a few interesting SNPs and check IGV

**High tech**

  * `bcftools stats` & `plot-vcfstats` (read the manual)


## Variant cleaning


![Not all variants look spectacular](attachments/igv_bad.png)

---

## Variant cleaning

Further cleaning is required - we use the [`vt`](https://genome.sph.umich.edu/wiki/Vt) toolkit

- **filter** remove all calls with a `QUAL<20`
- **decompose** takes multiallelic variants and splits them into
  multiple monoallelic variants. This makes it easier to assess the
  importance of either variant.
- **normalize** ensures that nucleotides around small inserts and
  deletions are properly trimmed & aligned
- **uniq** remove duplicate variants

---

## Variant normalization

![Variant normalization](attachments/snpnorm.png)

---

# Variant Annotation

## Variant annotaion

We have many variants - which ones are interesting?

* Mutation effect
* Population frequency
* GWAS significant
* Known clinical association
* Gene expression
* Protein interaction

## Disrupt protein coding genes?

![Effect on protein coding genes](attachments/effects.png)

---

## Sources of annotation

![Sources of variant information](attachments/asources.png)

---

## Variant Annotation Tools

* [Annovar](http://www.openbioinformatics.org/annovar)
* [Ensembl Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html)
* [**SnpEff**](https://pcingola.github.io/SnpEff/)

---

## SpnEff

* **Per SNP & splice variant**
* **Gene affected**
* **Severity**

    * HIGH
    * MODERATE
    * LOW
    * MODIFIER

---
## SpnEff

* **Transcript Biotype**
* **cDNA & protein position**
* **Effect Type**

    * missense_variant
    * nonsense_variant
    * intron_variant
    * frameshift_variant
    * [more](https://pcingola.github.io/SnpEff/se_inputoutput/)


---

# Using IGV

## Using IGV

Go to the [IGV web app](https://igv.org/app/)

![IGV web app](attachments/igv1.png)

* **Ensure the database is correct (hg38!)**

---

## Using IGV

* Download the `bam`, `bam.bai` and `vcf` files to your computer
* Select `Tracks` / `Local Files`

![IGV select files to upload|200](attachments/igv2.png)

* And select all `bam`, `bai` (and `vcf` if you have them) files in one go.
* Note - you need to zoom in to see something - find a location to
  look at in the `bam` or `vcf`
* In our case - check out the `NOTCH1` gene

---

# Now with real data

## T-ALL case

![T-Cell Acute Lymphoblastic Leukemia](attachments/case.png)

---

## T-ALL case

* Two samples - Tumor cells & Normal cells
* Only data from a small region from chromosome 9 (to keep it manageable)
* Using Grch38 genome build

**To Start:**

* Start a jupyterhub job
* Get from toledo: `manual_snp_calling_workflow.ipynb`
* Upload to you jupyterhub - follow the notebook
* Goal
  * Understand the workflow
  * Study the in- and output of every step
  * Think about *controls* - how do you know everything is still ok?

* Future: we will automate this exact same workflow using [SnakeMake](https://snakemake.readthedocs.io/en/stable/)
