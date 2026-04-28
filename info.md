* **Quality string**: GIves the base-bu-base confidence of the sequencer's calls
* A → ASCII 65 → Q = 32
* F → ASCII 70 → Q = 37
* K → ASCII 75 → Q = 42

ASCII-encoded **Phred quality scores**: Higher scores = higher confidence

$$Q=\frac{-10}{\log_{10}(p)}$$

* $Q20 \rightarrow 1$% error
* $Q30 \rightarrow 0.1$% error
* $Q40 \rightarrow 0.01$% error

Used for:
- QC assessment
- Timming low-quality ends
- filtering poor reads
- helping aligners weigh confidence

---

**Adapter contamination**: when a sequencing read contains part of the synthetic adapter sequence that was added during klibrary preparation (instead of only biological sequence from the sample). Occurs when DNA/RNA fragment is shorter than the read length so the sequencer reads through the real insert and into the adapter on the end. It can also happen it the adapter dimerizes or if poorly size-selected fragments get sequenced.

- Adapter contamination can reduce alignment rates and increase mismatches
- Can distort downstream analyses
- Often appear as rising adapter content signal toward 3-prime end of reads

Illumina adapter carryover e.g., `GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAAATCTCGTATGC`

homopolymer-rich sequences e.g.,
`CCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT`

---

### Alignment-based quantification 
**Tools**:
- `STAR`: Commonly used: fast and performrs well for splice aware alignmenmt
- `HISAT2`

e.g.,

```
trimmed FASTQ
    ↓
STAR alignment to reference genome
    ↓
BAM file
    ↓
alignment QC
    ↓
gene-level counting
```

```
STAR --genomeDir STAR_index \
     --readFilesIn sample_trimmed.fastq.gz \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix sample_
```

### Lightweight transcript quantification
**Tools**:
- `Salmon`
- `Kallisto`

e.g., 
```
trimmed FASTQ
    ↓
Salmon quantification against transcriptome
    ↓
transcript-level abundance estimates
    ↓
gene-level summarization with tximport
    ↓
DESeq2 / edgeR / limma-voom
```

Forward: `--libType SF`

Reverse: `--libType SR`

---

Single-end RNA-seq = sequencer reads only one end of each cDNA fragment:

```
RNA fragment / cDNA fragment:
5' ----------------------------- 3'

Single-end sequencing:
READ →
5' -----> 
```
  * Lower mapping confidence
  * Lower isoform/splice-junction resolution
  * Weaker fusion detection
  * Gene-level differentiation expression is fine (same for paired-end sequencing)

Paired-end sequencing = both ends of the same fragment are sequenced:

```
Paired-end sequencing:
READ 1 →                 ← READ 2
5' ----->             <----- 3'
```

  * Higher mappoing confidence
  * Higher isoform/splice-junction resolution
  * Better fusion detection

Stranded-RNA-seq lets aligner interpret whether a read supports gene on the same strand or the opposite strand

---

**Lane**:Physical sequencing lane on Illumina flow cell. A single biological sample library is often split across multiple lanes to get more depth and balance the run:

e.g., 
```
genelab/extracted_rnaseq/C-FLT-4/C-FLT-4_S62_L005_R1_001.fastq.gz
genelab/extracted_rnaseq/C-FLT-4/C-FLT-4_S62_L006_R1_001.fastq.gz
genelab/extracted_rnaseq/C-FLT-4/C-FLT-4_S62_L007_R1_001.fastq.gz 
```

---

1. Build HIST2 genome index
2. Merge trimmed lane FASTQs into FASTQ per biological sample before alignment
3. Alignment:

e.g Output., 

* `C-Ba-1.hisat2.summary.txt`: $35,633,737$ input reads
* $29,943,499$ reads ($84.03$%) aligned *once*
* $5,164,864$ reads ($14.49$%), aligned once (multimapping)
* Only $525,374$ reads ($1.47$) failed to align
* Overall alignment rate $=98.53$%

- In short-read single-end RNA-seq, especially at 50 bp, some reads will end up ambiguously in repetitive regions, gene families, pseudogenes, low-complexity sequence (common with 50bp)
  - Mouse transcriptomes contain repetitive and homologous regions
  - `525374 (1.47%) aligned 0 times` Little leftover contamination, adapter artifact, major reference mismatch

---

`-U` Single-end reads

`--rna-strandness R` Dataset is reverse-stranded

---

4. Feature count: Count single-end reads against Ensembl exon features and summariz them to genes

* Assigned counts: $25.8M$ to $30.7M$ (assigned fractions: $56.3$% and $58.1$%)
* Default option: not count multiparameters
* Wrong strands settings would drop instead of clustering around same value in every BAM.
* Short reads are more likely to match multiple loci, paralogs, repeats, pseudogenes, or homologous transcript regions
* Samples such as `C-FLT-4` have more assignedm reads because they're deeper libraries

5. Normalization

* Normalization in DESeq2: Accounts for sequencing depth and library composition $\rightarrow$ ensure gene expression levels are comparable across samples.
  * Geometric average is calculated with logs
  * Filter out `infinity` $\rightarrow$ Helps focus scaling factors on house keeping genes - scaling factors are based only on 'stable' genes that aren't jumping between 0 and high expression. Once scaling factors are determined they are then applied to **all** genes including tissue-specific ones
  * Substract average log value from `log(counts)` $\rightarrow$ Ratio reads in each sample to avg. across all samples
  
6. Calculate median of ratios for each sample (avoids large/outliers taking over)

* DESeq2 and EdgeR are both based on negative bionomial modeling
  * EdgeR better at analyzing genes with low expression counts (dispersion estimation captures variability in t he *sparse* count data). 

---

Multicellular program = Coordinated gene-expression pattern shared across multiple interacting cell types in a tissue, organoid, tumor, developmental system.

---

### Proteomics

* Thermo Fisher instruments produce `.raw`
  * contents:
    * spectra
    * scan metadata
    * chromatographic information
    * precursor/fragment ion information
    * instrument settings
* `.raw` $\rightarrow$ `.mzML`

`.mzML` $\rightarrow$ Fragpipe