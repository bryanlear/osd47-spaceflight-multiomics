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