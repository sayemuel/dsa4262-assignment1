# Task 4: transcript discovery and quantification with Bambu on the 4 Task 3 bam files
suppressMessages({
  library(bambu)
  library(GenomicRanges)
})

wd      <- "~/task4"
fa.file <- "~/task3/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
gtf.file <- "~/task3/reference/Homo_sapiens.GRCh38.91.gtf"
bams    <- list.files("~/task3/bam", pattern = "\\.bam$", full.names = TRUE)
dir.create(file.path(wd, "bambu_out"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(wd, "rc_cache"), showWarnings = FALSE)
print(bams)

# ---- run Bambu (default settings, with reference annotations) ----
t0 <- Sys.time()
annotations <- prepareAnnotations(gtf.file)
se <- bambu(reads = bams, annotations = annotations, genome = fa.file,
            ncore = as.integer(Sys.getenv("NCORE", "4")))
cat("bambu runtime:", format(Sys.time() - t0), "\n")
saveRDS(se, file.path(wd, "se.rds"))
writeBambuOutput(se, path = file.path(wd, "bambu_out"))

# ---- 4.1 ----
cat("\n==== 4.1 show(se) ====\n"); show(se)
cat("assays:", paste(assayNames(se), collapse = ", "), "\n")

# gene names from the GTF (gene_id -> gene_name), for readable answers
g <- rtracklayer::import(gtf.file, feature.type = "gene")
gene_names <- setNames(g$gene_name, g$gene_id)

rd     <- rowData(se)
counts <- assays(se)$counts
total  <- rowSums(counts)

# ---- 4.2a novel transcripts ----
novel <- rd$novelTranscript
cat("\n==== 4.2a ====\n")
cat("novel transcripts:", sum(novel), "\n")
cat("  of which in novel genes:", sum(novel & rd$novelGene), "\n")
cat("  of which in known genes:", sum(novel & !rd$novelGene), "\n")
cat("novel genes:", length(unique(rd$GENEID[rd$novelGene])), "\n")

# ---- 4.2b top 5 genes (sum of counts over all samples) ----
se_gene <- transcriptToGeneExpression(se)
saveRDS(se_gene, file.path(wd, "se_gene.rds"))
gtot <- sort(rowSums(assays(se_gene)$counts), decreasing = TRUE)
cat("\n==== 4.2b top 5 genes (total counts, all samples) ====\n")
top <- head(gtot, 5)
print(data.frame(gene_id = names(top), gene_name = gene_names[names(top)], total_counts = round(top, 1)))
cat("per-sample counts for these genes:\n")
print(round(assays(se_gene)$counts[names(top), ], 1))

# ---- 4.2c / 4.2d transcripts per gene ----
tx_per_gene <- function(keep) {
  n <- table(rd$GENEID[keep])
  c(min = min(n), max = max(n), mean = round(mean(n), 3), genes = length(n), transcripts = sum(n),
    max_gene = names(n)[which.max(n)])
}
cat("\n==== 4.2c transcripts per gene (all annotations after discovery) ====\n")
print(tx_per_gene(rep(TRUE, nrow(se))))
cat("\n==== 4.2d transcripts per gene (expressed: total counts >= 10) ====\n")
print(tx_per_gene(total >= 10))
cat("\n(alternative: >= 10 counts in at least one sample)\n")
print(tx_per_gene(apply(counts, 1, max) >= 10))

# ---- 4.2e single-exon novel transcripts ----
n_exons <- elementNROWS(rowRanges(se))
cat("\n==== 4.2e single-exon novel transcripts:", sum(novel & n_exons == 1), "====\n")

# ---- 4.2f highest expressed transcript with only unique read counts ----
cat("\n==== 4.2f ====\n")
if ("uniqueCounts" %in% assayNames(se)) {
  uniq <- rowSums(assays(se)$uniqueCounts)
  only_unique <- total > 0 & abs(total - uniq) < 1e-6
  cat("transcripts whose counts are all unique reads:", sum(only_unique), "\n")
  cand <- sort(total[only_unique], decreasing = TRUE)
  top_u <- head(cand, 5)
  print(data.frame(tx = names(top_u), gene_id = rd[names(top_u), "GENEID"],
                   gene_name = gene_names[rd[names(top_u), "GENEID"]],
                   total_counts = round(top_u, 1), unique_counts = round(uniq[names(top_u)], 1)))
} else cat("uniqueCounts assay not found\n")

# ---- material for 4.4: novel transcripts with evidence ----
nov_idx <- which(novel)
nov <- data.frame(
  tx = rownames(se)[nov_idx],
  gene_id = rd$GENEID[nov_idx],
  gene_name = gene_names[rd$GENEID[nov_idx]],
  novelGene = rd$novelGene[nov_idx],
  class = if ("txClassDescription" %in% colnames(rd)) rd$txClassDescription[nov_idx] else NA,
  n_exons = n_exons[nov_idx],
  coords = as.character(unlist(range(rowRanges(se)[nov_idx]))),
  total_counts = round(total[nov_idx], 1),
  full_length = round(rowSums(assays(se)$fullLengthCounts)[nov_idx], 1),
  round(counts[nov_idx, , drop = FALSE], 1),
  check.names = FALSE)
if ("txNDR" %in% colnames(rd)) nov$txNDR <- rd$txNDR[nov_idx]
nov <- nov[order(-nov$total_counts), ]
write.csv(nov, file.path(wd, "novel_transcripts.csv"), row.names = FALSE)
cat("\n==== top novel transcripts ====\n"); print(head(nov, 15))
cat("\nDONE\n")
