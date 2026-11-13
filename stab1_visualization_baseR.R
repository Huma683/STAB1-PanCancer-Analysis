options(stringsAsFactors = FALSE)

message("Starting STAB1 visualization script...")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The R package 'survival' is missing.")
}

library(survival)

# -----------------------------
# Helper functions
# -----------------------------

safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

theme_pdf <- function() {
  par(
    mar = c(5, 10, 4, 2),
    las = 1,
    cex.axis = 0.85,
    cex.lab = 0.9,
    cex.main = 1.0
  )
}

forest_plot <- function(tab, outfile, title, label_col = "label") {
  if (nrow(tab) == 0) stop("No data for forest plot: ", outfile)
  
  tab$HR_high_vs_low <- safe_num(tab$HR_high_vs_low)
  tab$lower95 <- safe_num(tab$lower95)
  tab$upper95 <- safe_num(tab$upper95)
  tab$cox_p_value <- safe_num(tab$cox_p_value)
  
  tab <- tab[!is.na(tab$HR_high_vs_low) & !is.na(tab$lower95) & !is.na(tab$upper95), ]
  tab <- tab[order(tab$HR_high_vs_low), ]
  
  y <- seq_len(nrow(tab))
  xlim <- range(c(tab$lower95, tab$upper95, 1), na.rm = TRUE)
  xlim[1] <- max(0.05, xlim[1] * 0.8)
  xlim[2] <- xlim[2] * 1.2
  
  pdf(outfile, width = 9, height = max(5, 0.35 * nrow(tab) + 2))
  theme_pdf()
  plot(
    tab$HR_high_vs_low,
    y,
    xlim = xlim,
    ylim = c(0.5, nrow(tab) + 0.5),
    log = "x",
    yaxt = "n",
    pch = 19,
    xlab = "Hazard ratio, high STAB1 vs low STAB1",
    ylab = "",
    main = title
  )
  segments(tab$lower95, y, tab$upper95, y, lwd = 1.6)
  abline(v = 1, lty = 2)
  axis(2, at = y, labels = tab[[label_col]], cex.axis = 0.75)
  
  ptxt <- ifelse(
    is.na(tab$cox_p_value),
    "",
    paste0("p=", format(tab$cox_p_value, digits = 2, scientific = TRUE))
  )
  text(xlim[2], y, ptxt, pos = 2, cex = 0.65)
  dev.off()
}

# -----------------------------
# Read existing result tables
# -----------------------------

pan_file <- "results/full/STAB1_pan_cancer_all_endpoints_cox.tsv"
cancer_file <- "results/full/STAB1_cancer_wise_all_endpoints_cox.tsv"
key_file <- "results/full/STAB1_key_cancers_all_endpoints_summary.tsv"
expr_file <- "results/full/STAB1_TCGA_tumor_vs_normal_wilcox.tsv"
merged_file <- "results/full/STAB1_merged_expression_clinical_all_endpoints.tsv"

if (!file.exists(pan_file)) stop("Missing: ", pan_file)
if (!file.exists(cancer_file)) stop("Missing: ", cancer_file)
if (!file.exists(key_file)) stop("Missing: ", key_file)
if (!file.exists(expr_file)) stop("Missing: ", expr_file)
if (!file.exists(merged_file)) stop("Missing: ", merged_file)

pan <- read.delim(pan_file, check.names = FALSE)
cancer <- read.delim(cancer_file, check.names = FALSE)
key <- read.delim(key_file, check.names = FALSE)
expr <- read.delim(expr_file, check.names = FALSE)
merged <- read.delim(merged_file, check.names = FALSE)

# Convert numeric columns
num_cols <- c("n", "events", "median_cutoff_STAB1", "HR_high_vs_low",
              "lower95", "upper95", "cox_p_value", "logrank_p_value", "FDR_p_value")

for (v in intersect(num_cols, names(pan))) pan[[v]] <- safe_num(pan[[v]])
for (v in intersect(num_cols, names(cancer))) cancer[[v]] <- safe_num(cancer[[v]])
for (v in intersect(num_cols, names(key))) key[[v]] <- safe_num(key[[v]])

# -----------------------------
# Figure 1: Pan-cancer forest plot
# -----------------------------

pan$label <- paste(pan$endpoint, pan$model, sep = " | ")

forest_plot(
  pan,
  "results/figures/Figure1_pan_cancer_all_endpoints_forest.pdf",
  "STAB1 pan-cancer survival across endpoints",
  label_col = "label"
)

# -----------------------------
# Figure 2: Key-cancer OS forest plot
# -----------------------------

key_os <- key[key$endpoint == "OS", ]
key_os$label <- paste(key_os$cancer_type, key_os$model, sep = " | ")

forest_plot(
  key_os,
  "results/figures/Figure2_key_cancers_OS_forest.pdf",
  "STAB1 overall survival in key cancers",
  label_col = "label"
)

# -----------------------------
# Figure 3: Cancer-wise HR heatmap
# -----------------------------

message("Creating HR heatmap...")

heat <- cancer[cancer$model == "Univariate Cox", ]
heat <- heat[!is.na(heat$HR_high_vs_low), ]

endpoints <- c("OS", "DSS", "DFI", "PFI")
cancers <- sort(unique(heat$cancer_type))

mat <- matrix(NA_real_, nrow = length(cancers), ncol = length(endpoints))
rownames(mat) <- cancers
colnames(mat) <- endpoints

for (i in seq_len(nrow(heat))) {
  ct <- heat$cancer_type[i]
  ep <- heat$endpoint[i]
  if (ct %in% rownames(mat) && ep %in% colnames(mat)) {
    mat[ct, ep] <- log2(heat$HR_high_vs_low[i])
  }
}

# Keep cancers with at least one endpoint
mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

# Order cancers by OS log2 HR if available
if ("OS" %in% colnames(mat)) {
  ord <- order(mat[, "OS"], na.last = TRUE)
  mat <- mat[ord, , drop = FALSE]
}

pdf("results/figures/Figure3_cancer_wise_HR_heatmap.pdf", width = 7, height = max(6, 0.25 * nrow(mat) + 2))
par(mar = c(5, 8, 4, 5), las = 1)
zlim <- c(-2, 2)
cols <- colorRampPalette(c("blue", "white", "red"))(101)

image(
  x = seq_len(ncol(mat)),
  y = seq_len(nrow(mat)),
  z = t(mat),
  col = cols,
  zlim = zlim,
  xaxt = "n",
  yaxt = "n",
  xlab = "Survival endpoint",
  ylab = "",
  main = "Cancer-wise STAB1 hazard ratio heatmap\nlog2(HR), high vs low STAB1"
)
axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
axis(2, at = seq_len(nrow(mat)), labels = rownames(mat), cex.axis = 0.65)

# Add values
for (i in seq_len(nrow(mat))) {
  for (j in seq_len(ncol(mat))) {
    if (!is.na(mat[i, j])) {
      text(j, i, sprintf("%.2f", mat[i, j]), cex = 0.55)
    }
  }
}

legend_x <- ncol(mat) + 0.8
legend_y <- seq(1, nrow(mat), length.out = 5)
legend_labels <- c("HR 0.25", "HR 0.5", "HR 1", "HR 2", "HR 4")
legend_values <- c(-2, -1, 0, 1, 2)
dev.off()

# -----------------------------
# Figure 4: Tumour-normal STAB1 difference
# -----------------------------

message("Creating tumour-normal expression plot...")

expr$median_tumor <- safe_num(expr$median_tumor)
expr$median_normal <- safe_num(expr$median_normal)
expr$FDR_p_value <- safe_num(expr$FDR_p_value)
expr$delta_median <- expr$median_tumor - expr$median_normal

expr <- expr[order(expr$delta_median), ]

pdf("results/figures/Figure4_TCGA_tumor_normal_STAB1_difference.pdf", width = 8, height = max(5, 0.35 * nrow(expr) + 2))
par(mar = c(5, 8, 4, 2), las = 1)
barplot(
  expr$delta_median,
  names.arg = expr$cancer_type,
  horiz = TRUE,
  xlab = "Median STAB1 difference: tumour minus normal",
  main = "STAB1 tumour-normal expression difference in TCGA",
  cex.names = 0.75
)
abline(v = 0, lty = 2)
dev.off()

# -----------------------------
# Figure 5: Key-cancer KM plots for OS
# -----------------------------

message("Creating key-cancer Kaplan-Meier plots...")

key_cancers <- c("KIRC", "UVM", "OV", "GBM", "LUSC", "ACC", "PRAD", "SKCM")
merged$OS_time <- safe_num(merged$OS_time)
merged$OS_event <- safe_num(merged$OS_event)
merged$STAB1 <- safe_num(merged$STAB1)

pdf("results/figures/Figure5_key_cancers_OS_KM_plots.pdf", width = 7, height = 6)

for (ct in key_cancers) {
  d <- merged[merged$cancer_type == ct &
                !is.na(merged$STAB1) &
                !is.na(merged$OS_time) &
                !is.na(merged$OS_event) &
                merged$OS_time > 0, ]
  
  if (nrow(d) < 40 || sum(d$OS_event, na.rm = TRUE) < 10) next
  
  med <- median(d$STAB1, na.rm = TRUE)
  d$STAB1_group <- ifelse(d$STAB1 >= med, "High", "Low")
  d$STAB1_group <- factor(d$STAB1_group, levels = c("Low", "High"))
  
  fit <- survfit(Surv(OS_time, OS_event) ~ STAB1_group, data = d)
  lr <- survdiff(Surv(OS_time, OS_event) ~ STAB1_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  
  plot(
    fit,
    col = c("blue", "red"),
    lwd = 2,
    xlab = "Days",
    ylab = "Overall survival probability",
    main = paste0(ct, ": STAB1 high vs low OS\nlog-rank p = ", signif(lr_p, 3))
  )
  legend(
    "topright",
    legend = c("Low STAB1", "High STAB1"),
    col = c("blue", "red"),
    lwd = 2,
    bty = "n"
  )
}

dev.off()

# -----------------------------
# Figure 6: Pan-cancer KM plot for OS
# -----------------------------

message("Creating pan-cancer Kaplan-Meier plot...")

d <- merged[
  !is.na(merged$STAB1) &
    !is.na(merged$OS_time) &
    !is.na(merged$OS_event) &
    merged$OS_time > 0,
]

med <- median(d$STAB1, na.rm = TRUE)
d$STAB1_group <- ifelse(d$STAB1 >= med, "High", "Low")
d$STAB1_group <- factor(d$STAB1_group, levels = c("Low", "High"))

fit <- survfit(Surv(OS_time, OS_event) ~ STAB1_group, data = d)
lr <- survdiff(Surv(OS_time, OS_event) ~ STAB1_group, data = d)
lr_p <- 1 - pchisq(lr$chisq, length(lr$n) - 1)

pdf("results/figures/Figure6_pan_cancer_OS_KM.pdf", width = 7, height = 6)
plot(
  fit,
  col = c("blue", "red"),
  lwd = 2,
  xlab = "Days",
  ylab = "Overall survival probability",
  main = paste0("Pan-cancer TCGA: STAB1 high vs low OS\nlog-rank p = ", signif(lr_p, 3))
)
legend(
  "topright",
  legend = c("Low STAB1", "High STAB1"),
  col = c("blue", "red"),
  lwd = 2,
  bty = "n"
)
dev.off()

# -----------------------------
# Save figure index
# -----------------------------

fig_index <- data.frame(
  figure = paste0("Figure ", 1:6),
  file = c(
    "Figure1_pan_cancer_all_endpoints_forest.pdf",
    "Figure2_key_cancers_OS_forest.pdf",
    "Figure3_cancer_wise_HR_heatmap.pdf",
    "Figure4_TCGA_tumor_normal_STAB1_difference.pdf",
    "Figure5_key_cancers_OS_KM_plots.pdf",
    "Figure6_pan_cancer_OS_KM.pdf"
  ),
  purpose = c(
    "Pan-cancer HRs for OS, DSS, DFI, and PFI",
    "Key-cancer OS HRs",
    "Cancer-wise endpoint heatmap of log2 HR",
    "Tumour-normal STAB1 median expression difference",
    "Kaplan-Meier OS curves for key cancers",
    "Pan-cancer Kaplan-Meier OS curve"
  )
)

write.table(
  fig_index,
  "results/figures/STAB1_figure_index.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Visualization complete.")
message("Figures created in: results/figures")
