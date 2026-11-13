options(stringsAsFactors = FALSE)

dir.create("results/geo_validation", recursive = TRUE, showWarnings = FALSE)

geo <- read.delim("data/manual_validation/STAB1_GEO_validation.tsv")
geo$HR <- as.numeric(geo$HR)
geo$lower95 <- as.numeric(geo$lower95)
geo$upper95 <- as.numeric(geo$upper95)
geo$p_value <- as.numeric(geo$p_value)

write.table(
  geo,
  "results/geo_validation/STAB1_GEO_validation_clean.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

geo2 <- geo[!is.na(geo$HR), ]
geo2$label <- paste(geo2$cancer, geo2$dataset, geo2$endpoint, sep = " | ")
geo2 <- geo2[order(geo2$HR), ]

pdf("results/geo_validation/Figure_GEO_validation_HR_plot.pdf", width = 9, height = 6)
par(mar = c(5, 10, 4, 2), las = 1)

y <- seq_len(nrow(geo2))
xlim <- range(c(geo2$HR, geo2$lower95, geo2$upper95, 1), na.rm = TRUE)
xlim[1] <- max(0.1, xlim[1] * 0.7)
xlim[2] <- xlim[2] * 1.2

plot(
  geo2$HR,
  y,
  log = "x",
  xlim = xlim,
  yaxt = "n",
  pch = 19,
  xlab = "Hazard ratio",
  ylab = "",
  main = "Independent GEO validation of STAB1 survival association"
)
segments(geo2$lower95, y, geo2$upper95, y)
abline(v = 1, lty = 2)
axis(2, at = y, labels = geo2$label, cex.axis = 0.7)
text(xlim[2], y, paste0("p=", signif(geo2$p_value, 2)), pos = 2, cex = 0.7)
dev.off()

cat("GEO validation complete.\n")
