# ==============================================================================
# 0. Setup and Initial Settings
# ==============================================================================
cat("--- STEP 0: Setup and Library Loading ---\n")

# setwd("D:/peaks/analysis")

# --- Load Core Libraries ---
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(ggplot2)
library(ggrepel)
library(limma)     # Differential analysis (limma-voom)
library(edgeR)     # Differential analysis (DGEList and TMM normalization)
library(DESeq2)    # Differential analysis (Negative Binomial model)
library(Rsamtools) # Read BAM files
library(GenomicAlignments) # Count reads
library(ChIPseeker) # Peak Annotation
library(clusterProfiler) # Gene Set Enrichment
library(enrichplot)      # Gene Set Enrichment Visualization
library(msigdbr)     # MSigDB collections (Hallmark, GO)
library(TFBSTools)   # Motif retrieval
library(JASPAR2020)  # Motif database
library(motifmatchr) # Motif scanning
library(ggfortify)   # PCA plotting

# --- Load Genome-Specific Libraries (Mouse mm10) ---
# Ensure these are installed: BiocManager::install(c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db", "BSgenome.Mmusculus.UCSC.mm10"))
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)


# ==============================================================================
# 1. Peak Merging and Read Counting (Unified Consensus Peak Set)
# ==============================================================================
cat("\n--- STEP 1: Peak Merging and Read Counting ---\n")

# A. List your MACS2 narrowPeak files
peak_files <- c(
  "Mel_Day0_1_peaks.narrowPeak",
  "Mel_Day0_2_peaks.narrowPeak",
  "Mel_Day2_1_peaks.narrowPeak",
  "Mel_Day2_2_peaks.narrowPeak"
)

# B. Read each as a GRanges object and combine
peak_list <- lapply(peak_files, rtracklayer::readPeakFile)
all_peaks <- do.call(c, peak_list)

# C. Reduce: Merge all overlapping/nearby peaks into one consensus region
# Use min.gapwidth=100 to merge peaks within 100 bp of each other
merged_peaks <- GenomicRanges::reduce(all_peaks, min.gapwidth = 100)

# D. Export consensus peaks to BED format
rtracklayer::export(merged_peaks, "Merged_peaks.bed", format = "BED")
cat("✅ Merged_peaks.bed generated successfully! Number of peaks:", length(merged_peaks), "\n")

# E. Define BAM files and create BamFileList
bam_files <- c("Day0_1.bam", "Day0_2.bam", "Day2_1.bam", "Day2_2.bam")
bam_list <- Rsamtools::BamFileList(bam_files, yieldSize = 2000000)

# F. Count reads overlapping each merged peak
peak_counts_SO <- GenomicAlignments::summarizeOverlaps(
  features = merged_peaks,
  reads = bam_list,
  mode = "Union",
  singleEnd = TRUE,      # Change to FALSE if paired-end
  ignore.strand = TRUE
)

# G. Extract count matrix and add identifiers
count_matrix <- assay(peak_counts_SO)
colnames(count_matrix) <- c("Day0_1", "Day0_2", "Day2_1", "Day2_2")
peak_names <- paste0(seqnames(merged_peaks), "_", start(merged_peaks), "_", end(merged_peaks))
rownames(count_matrix) <- peak_names

# H. Save the raw count table
write.csv(count_matrix, "MACS2_Merged_Peak_ReadCounts.csv", quote = FALSE)
peak_counts <- count_matrix # Assign for consistency

# ==============================================================================
# 2. Data Preparation, Filtering, and Annotation (Shared steps)
# ==============================================================================
cat("\n--- STEP 2: Data Preparation, Filtering, and Annotation ---\n")

# A. Create peak annotation table
peak_annot <- data.frame(
  PeakID = rownames(peak_counts),
  Chr = seqnames(merged_peaks),
  Start = start(merged_peaks),
  End = end(merged_peaks)
)
rownames(peak_annot) <- peak_annot$PeakID

# B. Create sample metadata
meta_peak <- data.frame(
  SampleID = sample_names,
  Condition = conditions,
  stringsAsFactors = FALSE
)
rownames(meta_peak) <- meta_peak$SampleID

# C. Filter low-count peaks (keep peaks with counts > 0 in at least 50% of samples)
perc_keep <- 0.5
peak_keep <- rowSums(peak_counts > 0) >= ceiling(perc_keep * ncol(peak_counts))
peak_counts_filtered <- peak_counts[peak_keep, ]
peak_annot_filtered <- peak_annot[rownames(peak_counts_filtered), ]

cat("Peaks after filtering:", nrow(peak_counts_filtered), "\n")

# D. Peak Annotation (Run once for all analyses)
peak_gr_filtered <- GRanges(
  seqnames = peak_annot_filtered$Chr,
  ranges = IRanges(start = peak_annot_filtered$Start, end = peak_annot_filtered$End),
  names = peak_annot_filtered$PeakID
)
peak_anno <- ChIPseeker::annotatePeak(
  peak_gr_filtered,
  tssRegion = c(-10000, 20000),
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  annoDb = "org.Mm.eg.db"
)
anno_df <- as.data.frame(peak_anno)
anno_df$PeakKey <- paste0(anno_df$seqnames, "_", anno_df$start, "_", anno_df$end)


# ==============================================================================
# 3A. Differential Binding Analysis (limma-voom/edgeR)
# ==============================================================================
cat("\n--- STEP 3A: Differential Binding Analysis (limma-voom/edgeR) ---\n")
padj_cutoff <- 0.05
lfc_cutoff <- 1

# 1. CREATE DGEList, NORMALIZE (TMM), and VOOM
dge_peak <- edgeR::DGEList(counts = peak_counts_filtered, samples = meta_peak, genes = peak_annot_filtered)
dge_peak <- edgeR::calcNormFactors(dge_peak, method = "TMM")
pdf("limma_voom_mean_variance_plot.pdf")
dge_peak_voom <- limma::voom(dge_peak, plot = TRUE)
dev.off()

# 2. Fit linear model
design <- model.matrix(~ 0 + Condition, data = dge_peak_voom$targets)
colnames(design) <- gsub("Condition", "", colnames(design))
contrast_matrix <- limma::makeContrasts(comparison = Day2 - Day0, levels = design)
fit <- limma::lmFit(dge_peak_voom, design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

# 3. Get results and annotate
depeak_tbl_limma <- limma::topTable(fit2, coef = "comparison", n = Inf, sort.by = "P")
depeak_tbl_limma$PeakKey_merge <- rownames(depeak_tbl_limma)
depeak_tbl_limma <- merge(
  depeak_tbl_limma,
  anno_df[, c("PeakKey", "SYMBOL", "annotation")],
  by.x = "PeakKey_merge", by.y = "PeakKey", all.x = TRUE
)
colnames(depeak_tbl_limma)[colnames(depeak_tbl_limma) == "SYMBOL"] <- "Gene"

# 4. Save results
fwrite(depeak_tbl_limma, "MACS2_DE_Peaks_Day2_vs_Day0_limma_Annotated.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 5. PCA Plot
count_matrix_t <- as.data.frame(t(dge_peak_voom$E))
pca_result <- prcomp(count_matrix_t, scale. = TRUE)
var_explained <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
pca_plot_limma <- ggfortify::autoplot(
  pca_result, data = meta_peak, colour = "Condition", shape = "Condition", label = TRUE
) + labs(title = "PCA of limma-voom Normalized Counts")
ggsave("PCA_MACS2_Samples_limma.pdf", pca_plot_limma, width = 20, height = 16, units = "cm")

# 6. Volcano Plot
depeak_tbl_limma$Regulation <- "Not Significant"
depeak_tbl_limma$Regulation[depeak_tbl_limma$logFC > lfc_cutoff & depeak_tbl_limma$adj.P.Val < padj_cutoff] <- "Up (Day2 > Day0)"
depeak_tbl_limma$Regulation[depeak_tbl_limma$logFC < -lfc_cutoff & depeak_tbl_limma$adj.P.Val < padj_cutoff] <- "Down (Day2 < Day0)"
top_genes <- rbind(
  head(depeak_tbl_limma[depeak_tbl_limma$Regulation == "Up (Day2 > Day0)" & !is.na(depeak_tbl_limma$Gene), ], 10),
  head(depeak_tbl_limma[depeak_tbl_limma$Regulation == "Down (Day2 < Day0)" & !is.na(depeak_tbl_limma$Gene), ], 10)
)
depeak_tbl_limma$Label <- ifelse(depeak_tbl_limma$PeakKey_merge %in% top_genes$PeakKey_merge, depeak_tbl_limma$Gene, NA)

volcano_plot_limma <- ggplot(depeak_tbl_limma, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey70", "Up (Day2 > Day0)" = "firebrick2", "Down (Day2 < Day0)" = "dodgerblue3")) +
  ggrepel::geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  labs(title = "Differential Binding: Day2 vs Day0 (limma-voom)") +
  theme_bw(base_size = 14)
ggsave("MACS2_VolcanoPlot_limma_GeneAnnotated.pdf", volcano_plot_limma, width = 10, height = 8)


# ==============================================================================
# 3B. Differential Binding Analysis (DESeq2)
# ==============================================================================
cat("\n--- STEP 3B: Differential Binding Analysis (DESeq2) ---\n")

# 1. Prepare DESeq2 dataset and run analysis
dds <- DESeqDataSetFromMatrix(countData = round(peak_counts_filtered), colData = meta_peak, design = ~ Condition)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds$Condition <- relevel(dds$Condition, ref = "Day0")
dds <- DESeq(dds)

# 2. Extract results and apply LFC shrinkage
res <- results(dds, contrast = c("Condition", "Day2", "Day0"))
res <- lfcShrink(dds, coef = "Condition_Day2_vs_Day0", res = res, type = "ashr")
res_df <- as.data.frame(res)
res_df$PeakKey_merge <- rownames(res_df)

# 3. Annotate and save results
depeak_tbl_deseq2 <- merge(
  res_df,
  anno_df[, c("PeakKey", "SYMBOL", "annotation")],
  by.x = "PeakKey_merge", by.y = "PeakKey", all.x = TRUE
)
colnames(depeak_tbl_deseq2)[colnames(depeak_tbl_deseq2) == "SYMBOL"] <- "Gene"
fwrite(depeak_tbl_deseq2, "MACS2_DESeq2_Day2_vs_Day0_Annotated.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 4. PCA Plot
rld <- rlog(dds, blind = TRUE)
pcaData <- plotPCA(rld, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot_deseq2 <- ggplot(pcaData, aes(PC1, PC2, color = Condition, shape = Condition)) +
  geom_point(size = 3) +
  labs(title = "PCA of DESeq2 rlog Normalized Counts",
       x = paste0("PC1 (", percentVar[1], "%)"), y = paste0("PC2 (", percentVar[2], "%)")) +
  theme_bw(base_size = 14)
ggsave("PCA_MACS2_DESeq2.pdf", pca_plot_deseq2, width = 8, height = 6)

# 5. Volcano Plot
depeak_tbl_deseq2$Regulation <- "Not Significant"
depeak_tbl_deseq2$Regulation[depeak_tbl_deseq2$log2FoldChange > lfc_cutoff & depeak_tbl_deseq2$padj < padj_cutoff] <- "Up (Day2 > Day0)"
depeak_tbl_deseq2$Regulation[depeak_tbl_deseq2$log2FoldChange < -lfc_cutoff & depeak_tbl_deseq2$padj < padj_cutoff] <- "Down (Day2 < Day0)"
top_genes_deseq2 <- rbind(
  head(depeak_tbl_deseq2[depeak_tbl_deseq2$Regulation == "Up (Day2 > Day0)" & !is.na(depeak_tbl_deseq2$Gene), ], 10),
  head(depeak_tbl_deseq2[depeak_tbl_deseq2$Regulation == "Down (Day2 < Day0)" & !is.na(depeak_tbl_deseq2$Gene), ], 10)
)
depeak_tbl_deseq2$Label <- ifelse(depeak_tbl_deseq2$PeakKey_merge %in% top_genes_deseq2$PeakKey_merge, depeak_tbl_deseq2$Gene, NA)

volcano_plot_deseq2 <- ggplot(depeak_tbl_deseq2, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey70", "Up (Day2 > Day0)" = "firebrick2", "Down (Day2 < Day0)" = "dodgerblue3")) +
  ggrepel::geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  labs(title = "Differential Binding: Day2 vs Day0 (DESeq2)") +
  theme_bw(base_size = 14)
ggsave("MACS2_VolcanoPlot_DESeq2_GeneAnnotated.pdf", volcano_plot_deseq2, width = 10, height = 8)


# ==============================================================================
# 4. Motif Enrichment Setup (Preparation for HOMER)
# ==============================================================================
cat("\n--- STEP 4: Motif Enrichment Setup (for HOMER) ---\n")
dir.create("Motifs", showWarnings = FALSE)

# --- A. Setup for limma-voom DE peaks ---
up_peaks_limma <- subset(depeak_tbl_limma, logFC > lfc_cutoff & adj.P.Val < padj_cutoff)
down_peaks_limma <- subset(depeak_tbl_limma, logFC < -lfc_cutoff & adj.P.Val < padj_cutoff)
up_gr_limma <- GRanges(seqnames = up_peaks_limma$Chr, ranges = IRanges(start = up_peaks_limma$Start, end = up_peaks_limma$End))
down_gr_limma <- GRanges(seqnames = down_peaks_limma$Chr, ranges = IRanges(start = down_peaks_limma$Start, end = down_peaks_limma$End))
rtracklayer::export(up_gr_limma, "Motifs/diffbind_up_peaks_limma.bed", format = "BED")
rtracklayer::export(down_gr_limma, "Motifs/diffbind_down_peaks_limma.bed", format = "BED")

# --- B. Setup for DESeq2 DE peaks ---
up_peaks_deseq2 <- subset(depeak_tbl_deseq2, log2FoldChange > lfc_cutoff & padj < padj_cutoff)
down_peaks_deseq2 <- subset(depeak_tbl_deseq2, log2FoldChange < -lfc_cutoff & padj < padj_cutoff)
up_gr_deseq2 <- GRanges(seqnames = up_peaks_deseq2$Chr, ranges = IRanges(start = up_peaks_deseq2$Start, end = up_peaks_deseq2$End))
down_gr_deseq2 <- GRanges(seqnames = down_peaks_deseq2$Chr, ranges = IRanges(start = down_peaks_deseq2$Start, end = down_peaks_deseq2$End))
rtracklayer::export(up_gr_deseq2, "Motifs/diffbind_up_peaks_DESeq2.bed", format = "BED")
rtracklayer::export(down_gr_deseq2, "Motifs/diffbind_down_peaks_DESeq2.bed", format = "BED")

cat("\n💡 HOMER commands for motif discovery (run in terminal/shell with mm10 genome):\n")
cat('findMotifsGenome.pl Motifs/diffbind_up_peaks_limma.bed mm10 Motifs/Up_Peaks_limma/ -size 200 -len 8,10,12 -mask -p 4\n')
cat('findMotifsGenome.pl Motifs/diffbind_down_peaks_limma.bed mm10 Motifs/Down_Peaks_limma/ -size 200 -len 8,10,12 -mask -p 4\n')


# ==============================================================================
# 5. Gene Set Enrichment Analysis (using limma-voom results)
# ==============================================================================
cat("\n--- STEP 5: Gene Set Enrichment Analysis ---\n")

# A. Extract Entrez IDs for differentially bound genes (Up-regulated binding)
up_genes <- subset(depeak_tbl_limma, Regulation == "Up (Day2 > Day0)")
up_genes_entrez <- na.omit(unique(up_genes$ENTREZID))
# Ensure ENTREZID column is properly included in annotation step (via ChIPseeker)

if (length(up_genes_entrez) > 0) {
  # B. Perform GO Enrichment
  ego_GO <- clusterProfiler::enrichGO(
    gene = up_genes_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP", # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  if (!is.null(ego_GO)) {
    pdf("GO_BP_Enrichment_Up_Peaks.pdf", width = 10, height = 8)
    print(clusterProfiler::dotplot(ego_GO, showCategory = 15, title = "GO Biological Process Enrichment (Up Peaks)"))
    dev.off()
  }
  
  # C. Perform Hallmark Gene Set Enrichment (using msigdbr)
  m_t2g <- msigdbr::msigdbr(species = "Mus musculus", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)
  hallmark_enrich <- clusterProfiler::enricher(
    gene = up_genes_entrez,
    TERM2GENE = m_t2g
  )
  if (!is.null(hallmark_enrich)) {
    pdf("Hallmark_Enrichment_Up_Peaks.pdf", width = 10, height = 8)
    print(clusterProfiler::dotplot(hallmark_enrich, showCategory = 15, title = "Hallmark Gene Set Enrichment (Up Peaks)"))
    dev.off()
  }
  cat("✅ GO and Hallmark enrichment performed on Up-regulated peaks.\n")
} else {
  cat("⚠️ No significant up-regulated peaks to perform Gene Set Enrichment.\n")
}


# ==============================================================================
# 6. Example: Single Motif Search (PPARG)
# ==============================================================================
cat("\n--- STEP 6: Example: Single Motif Search (PPARG) ---\n")

# Use a single sample peak file's summits as example for motif search
example_peak_file <- peak_files[1]
example_peaks_df <- read.delim(example_peak_file, comment.char="#")

# Get peak summits and extend them (e.g., to 100bp)
example_summits_gr <- GRanges(
  seqnames = example_peaks_df$chr,
  IRanges(example_peaks_df$abs_summit, example_peaks_df$abs_summit)
)
example_summits_gr <- GenomicRanges::resize(example_summits_gr, 100, fix = "center")

# Get PPARG PFM from JASPAR2020 and match motif
opts <- list(name = "PPARG")
pfm_list <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
pfm <- pfm_list[[1]]

PPARG_Motifs <- motifmatchr::matchMotifs(
  pfm, example_summits_gr, BSgenome.Mmusculus.UCSC.mm10, out = "matches", genome = "mm10"
)

PPARG_positive_peaks <- example_summits_gr[which(motifmatchr::motifMatches(PPARG_Motifs))]
cat("PPARG motif found in", length(PPARG_positive_peaks), "peaks from", basename(example_peak_file), "\n")


# ==============================================================================
# 7. Peak Overlap and Consensus Analysis
# ==============================================================================
cat("\n--- STEP 7: Peak Overlap and Consensus Analysis ---\n")

# A. Load all individual peak files into a GRangesList
macsPeaks_GRL <- GRangesList(lapply(peak_files, function(f) {
  peakDFtemp <- read.delim(f, comment.char = "#")
  GRanges(seqnames = peakDFtemp[, "chr"], IRanges(peakDFtemp[, "start"], peakDFtemp[, "end"]))
}))
names(macsPeaks_GRL) <- gsub("_peaks.narrowPeak", "", peak_files)

# B. Get the non-redundant set of all peaks (same as merged_peaks if min.gapwidth=0 was used)
allPeaksSet_nR <- GenomicRanges::reduce(unlist(macsPeaks_GRL))

# C. Create overlap matrix and store in the non-redundant peak set
overlap <- lapply(macsPeaks_GRL, function(x) allPeaksSet_nR %over% x)
overlapMatrix <- do.call(cbind, overlap)
colnames(overlapMatrix) <- names(macsPeaks_GRL)
mcols(allPeaksSet_nR) <- overlapMatrix

# D. Visualize overlaps with a Venn Diagram
pdf("VennDiagram_PeakOverlaps.pdf")
limma::vennDiagram(mcols(allPeaksSet_nR), main = "Peak Overlaps")
dev.off()

# E. Define High-Confidence (HC) Peaks (found in both replicates of a condition)
Day0_HC_Peaks <- allPeaksSet_nR[
  rowSums(as.data.frame(mcols(allPeaksSet_nR)[, c("Mel_Day0_1", "Mel_Day0_2")])) >= 2
]
Day2_HC_Peaks <- allPeaksSet_nR[
  rowSums(as.data.frame(mcols(allPeaksSet_nR)[, c("Mel_Day2_1", "Mel_Day2_2")])) >= 2
]

cat("High-Confidence Day0 peaks:", length(Day0_HC_Peaks), "\n")
cat("High-Confidence Day2 peaks:", length(Day2_HC_Peaks), "\n")

rtracklayer::export(Day0_HC_Peaks, "Day0_HC_Peaks.bed", format = "BED")
rtracklayer::export(Day2_HC_Peaks, "Day2_HC_Peaks.bed", format = "BED")

cat("\n--- SCRIPT EXECUTION COMPLETE ---\n")