options(stringsAsFactors = FALSE)

dir.create("results/hpa_validation", recursive = TRUE, showWarnings = FALSE)

hpa <- read.delim("data/manual_validation/STAB1_HPA_validation.tsv")

write.table(
  hpa,
  "results/hpa_validation/STAB1_HPA_validation_clean.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf("results/hpa_validation/Figure_HPA_validation_summary.pdf", width = 10, height = 6)
par(mar = c(4, 2, 4, 2))
plot.new()
title("Human Protein Atlas validation summary for STAB1")

y <- seq(0.9, 0.15, length.out = nrow(hpa))
text(0.02, 0.98, "Evidence", adj = 0, font = 2)
text(0.30, 0.98, "Result", adj = 0, font = 2)
text(0.58, 0.98, "Interpretation", adj = 0, font = 2)

text(0.02, y, hpa$evidence_type, adj = 0, font = 2, cex = 0.85)
text(0.30, y, hpa$result, adj = 0, cex = 0.85)
text(0.58, y, hpa$interpretation, adj = 0, cex = 0.80)
dev.off()

cat("HPA validation complete.\n")
