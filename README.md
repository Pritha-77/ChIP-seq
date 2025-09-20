
# PPARG ChIP-Seq Analysis Pipeline (Human hg38)

This repository contains a complete ChIP-Seq data analysis pipeline for **PPARG** using R and human genome reference **hg38**. It includes read filtering, alignment, peak calling, quality control, peak annotation, and visualization.


````

## 🔧 Requirements

### R Packages

```R
BiocManager::install(c(
  "biomaRt", "readxl", "writexl", "Rsubread", "Rsamtools", "ChIPQC",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "GenomeInfoDb",
  "ChIPseeker", "rtracklayer", "GenomicRanges", "ggplot2", "openxlsx", "ggupset", "viridis"
))
````

### Python (External)

* [MACS2](https://github.com/macs3-project/MACS) for peak calling
  Install via pip:

  ```bash
  pip install macs2
  ```

## 🔬 Pipeline Overview

### 1. Read Filtering

Using `FastqStreamer()` to filter reads by quality and presence of "N".

```R
# Filters reads with alphabetScore > 300 and < 10 'N's
writeFastq(filt2, "filtered.fastq", mode = "a")
```

### 2. Genome Indexing

```R
library(BSgenome.Hsapiens.UCSC.hg38)
buildindex("hg38_mainchrs", "hg38.mainChrs.fa", memory = 8000)
```

### 3. Read Alignment

```R
align("hg38_mainchrs", "filtered.fastq", output_format = "BAM")
```

### 4. Peak Calling (MACS2)

```bash
macs2 callpeak -t aligned.bam -c control.bam -n SampleName --outdir PeakDirectory
```

### 5. Quality Control (ChIPQC)

```R
ChIPQCsample("aligned.bam", annotation = "hg38", blacklist = "blacklist.bed")
```

### 6. Peak Annotation

```R
annotatePeak(peaks, tssRegion = c(-10000, 20000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
```

### 7. Visualization

* Annotation Pie Chart
* Annotation Bar Plot
* Distance to TSS Plot
* UpSet Plot

```R
plotAnnoPie(peakAnno)
ggplot(...) + geom_bar()
```

## 📊 Example Output

* `peaks_with_gene_name1.xlsx`: Annotated peak list
* `Annotation_Pie.png`, `BarPlot.png`: Distribution plots
