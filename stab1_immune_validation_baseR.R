options(stringsAsFactors = FALSE)

dir.create("results/immune_validation", recursive = TRUE, showWarnings = FALSE)

immune <- read.delim("data/manual_validation/STAB1_TIMER3_macrophage_key_cancers.tsv")

write.table(
  immune,
  "results/immune_validation/STAB1_TIMER3_macrophage_key_cancers_clean.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf("results/immune_validation/Figure_TIMER3_macrophage_correlations.pdf", width = 8, height = 6)
par(mar = c(5, 5, 4, 2))
ylim <- range(c(immune$Mac_TIMER_rho, immune$MacM2_QUANT_rho, immune$Mac_EPIC_rho), na.rm = TRUE)

plot(
  seq_len(nrow(immune)),
  immune$Mac_TIMER_rho,
  type = "b",
  pch = 19,
  ylim = ylim,
  xaxt = "n",
  xlab = "Cancer type",
  ylab = "Spearman correlation with STAB1",
  main = "STAB1 correlation with macrophage infiltration"
)
lines(seq_len(nrow(immune)), immune$MacM2_QUANT_rho, type = "b", pch = 17)
lines(seq_len(nrow(immune)), immune$Mac_EPIC_rho, type = "b", pch = 15)
axis(1, at = seq_len(nrow(immune)), labels = immune$cancer, las = 2)
abline(h = 0, lty = 2)
legend(
  "topleft",
  legend = c("TIMER macrophage", "QUANTISEQ M2 macrophage", "EPIC macrophage"),
  pch = c(19, 17, 15),
  lty = 1,
  bty = "n"
)
dev.off()

mean_rho <- data.frame(
  cancer = immune$cancer,
  mean_macrophage_rho = rowMeans(
    immune[, c("Mac_TIMER_rho", "MacM2_QUANT_rho", "Mac_EPIC_rho")],
    na.rm = TRUE
  )
)

mean_rho <- mean_rho[order(mean_rho$mean_macrophage_rho), ]

pdf("results/immune_validation/Figure_TIMER3_mean_macrophage_ranking.pdf", width = 8, height = 6)
par(mar = c(5, 6, 4, 2), las = 1)
barplot(
  mean_rho$mean_macrophage_rho,
  names.arg = mean_rho$cancer,
  horiz = TRUE,
  xlab = "Mean macrophage correlation across TIMER/QUANTISEQ/EPIC",
  main = "Ranking of STAB1-macrophage association"
)
abline(v = 0, lty = 2)
dev.off()

cat("Immune validation complete.\n")
