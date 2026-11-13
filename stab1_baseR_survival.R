options(stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

message("Starting STAB1 base-R survival analysis...")

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The R package 'survival' is missing. This package is usually installed with R.")
}

library(survival)

expr_file <- "data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
surv_file <- "data/Survival_SupplementalTable_S1_20171025_xena_sp"

if (!file.exists(expr_file)) stop("Missing expression file: ", expr_file)
if (!file.exists(surv_file)) stop("Missing survival file: ", surv_file)

if (file.info(expr_file)$size < 1000000) {
  stop("Expression file is too small. It may be a failed download.")
}

if (file.info(surv_file)$size < 10000) {
  stop("Survival file is too small. It may be a failed download.")
}

message("Extracting STAB1 row from expression matrix...")

stab1_file <- "data/STAB1_expression_row.tsv"

if (!file.exists(stab1_file) || file.info(stab1_file)$size == 0) {
  cmd <- paste(
    "zcat", shQuote(expr_file),
    "| awk 'NR==1 || $1==\"STAB1\"' >",
    shQuote(stab1_file)
  )
  status <- system(cmd)
  if (status != 0) stop("Failed to extract STAB1 row.")
}

stab1 <- read.delim(stab1_file, check.names = FALSE)

if (nrow(stab1) < 1) {
  stop("STAB1 row was not found. Check whether gene names are plain symbols.")
}

sample_ids <- colnames(stab1)[-1]
stab1_values <- as.numeric(stab1[1, -1])

expr_df <- data.frame(
  sample = sample_ids,
  STAB1 = stab1_values,
  stringsAsFactors = FALSE
)

expr_df$patient <- substr(expr_df$sample, 1, 12)
expr_df$sample_type_code <- substr(expr_df$sample, 14, 15)

expr_df$sample_type <- ifelse(
  expr_df$sample_type_code == "01", "Tumor",
  ifelse(expr_df$sample_type_code == "11", "Normal", "Other")
)

message("Reading survival file...")

surv <- read.delim(surv_file, check.names = FALSE)

message("Survival columns are:")
print(colnames(surv))

sample_col <- "sample"
os_event_col <- "OS"
os_time_col <- "OS.time"
cancer_col <- "cancer type abbreviation"

if (!sample_col %in% colnames(surv)) stop("Column 'sample' not found.")
if (!os_event_col %in% colnames(surv)) stop("Column 'OS' not found.")
if (!os_time_col %in% colnames(surv)) stop("Column 'OS.time' not found.")
if (!cancer_col %in% colnames(surv)) stop("Column 'cancer type abbreviation' not found.")

surv_df <- data.frame(
  patient = substr(surv[[sample_col]], 1, 12),
  cancer_type = surv[[cancer_col]],
  OS_event = as.numeric(surv[[os_event_col]]),
  OS_time = as.numeric(surv[[os_time_col]]),
  stringsAsFactors = FALSE
)

message("Merging STAB1 expression with survival data...")

tumor_expr <- expr_df[expr_df$sample_type == "Tumor", ]

main <- merge(
  tumor_expr,
  surv_df,
  by = "patient",
  all.x = FALSE,
  all.y = FALSE
)

main <- main[
  !is.na(main$STAB1) &
    !is.na(main$OS_time) &
    !is.na(main$OS_event) &
    main$OS_time > 0,
]

message("Merged tumor samples with survival data: ", nrow(main))

if (nrow(main) < 100) {
  stop("Too few samples after merging. Check sample IDs.")
}

write.table(
  main,
  "results/STAB1_TCGA_tumor_survival_dataset.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Running pan-cancer median-split survival analysis...")

median_stab1 <- median(main$STAB1, na.rm = TRUE)
main$STAB1_group <- ifelse(main$STAB1 >= median_stab1, "High", "Low")
main$STAB1_group <- factor(main$STAB1_group, levels = c("Low", "High"))

fit <- survfit(Surv(OS_time, OS_event) ~ STAB1_group, data = main)
cox <- coxph(Surv(OS_time, OS_event) ~ STAB1_group, data = main)
cox_summary <- summary(cox)

logrank <- survdiff(Surv(OS_time, OS_event) ~ STAB1_group, data = main)
logrank_p <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)

pan_result <- data.frame(
  analysis = "Pan-cancer OS",
  n = nrow(main),
  events = sum(main$OS_event, na.rm = TRUE),
  median_cutoff_STAB1 = median_stab1,
  HR_high_vs_low = cox_summary$coefficients[1, "exp(coef)"],
  lower95 = cox_summary$conf.int[1, "lower .95"],
  upper95 = cox_summary$conf.int[1, "upper .95"],
  cox_p_value = cox_summary$coefficients[1, "Pr(>|z|)"],
  logrank_p_value = logrank_p
)

write.table(
  pan_result,
  "results/STAB1_pan_cancer_OS_cox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf("results/STAB1_pan_cancer_OS_KM.pdf", width = 7, height = 6)
plot(
  fit,
  col = c("blue", "red"),
  lwd = 2,
  xlab = "Days",
  ylab = "Overall survival probability",
  main = paste0("Pan-cancer TCGA: STAB1 high vs low\nLog-rank p = ", signif(logrank_p, 3))
)
legend(
  "topright",
  legend = c("Low STAB1", "High STAB1"),
  col = c("blue", "red"),
  lwd = 2,
  bty = "n"
)
dev.off()

message("Running cancer-wise OS Cox analysis...")

cancers <- sort(unique(main$cancer_type))
cancer_results <- list()

for (ct in cancers) {
  dat <- main[main$cancer_type == ct, ]
  
  if (nrow(dat) < 40) next
  if (sum(dat$OS_event, na.rm = TRUE) < 10) next
  
  med <- median(dat$STAB1, na.rm = TRUE)
  dat$STAB1_group <- ifelse(dat$STAB1 >= med, "High", "Low")
  dat$STAB1_group <- factor(dat$STAB1_group, levels = c("Low", "High"))
  
  cx <- try(coxph(Surv(OS_time, OS_event) ~ STAB1_group, data = dat), silent = TRUE)
  if (inherits(cx, "try-error")) next
  
  sf <- survfit(Surv(OS_time, OS_event) ~ STAB1_group, data = dat)
  lr <- survdiff(Surv(OS_time, OS_event) ~ STAB1_group, data = dat)
  lr_p <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  
  s <- summary(cx)
  
  cancer_results[[ct]] <- data.frame(
    cancer_type = ct,
    n = nrow(dat),
    events = sum(dat$OS_event, na.rm = TRUE),
    median_cutoff_STAB1 = med,
    HR_high_vs_low = s$coefficients[1, "exp(coef)"],
    lower95 = s$conf.int[1, "lower .95"],
    upper95 = s$conf.int[1, "upper .95"],
    cox_p_value = s$coefficients[1, "Pr(>|z|)"],
    logrank_p_value = lr_p,
    stringsAsFactors = FALSE
  )
  
  safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  pdf(file.path("results", paste0("STAB1_", safe_ct, "_OS_KM.pdf")), width = 7, height = 6)
  plot(
    sf,
    col = c("blue", "red"),
    lwd = 2,
    xlab = "Days",
    ylab = "Overall survival probability",
    main = paste0(ct, ": STAB1 high vs low\nLog-rank p = ", signif(lr_p, 3))
  )
  legend(
    "topright",
    legend = c("Low STAB1", "High STAB1"),
    col = c("blue", "red"),
    lwd = 2,
    bty = "n"
  )
  dev.off()
}

cancer_table <- do.call(rbind, cancer_results)

cancer_table$FDR_p_value <- p.adjust(cancer_table$cox_p_value, method = "BH")

write.table(
  cancer_table,
  "results/STAB1_cancer_wise_OS_cox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Creating simple tumor/normal STAB1 summary...")

expr_df$cancer_code <- substr(expr_df$sample, 6, 7)

expr_summary <- aggregate(
  STAB1 ~ cancer_code + sample_type,
  data = expr_df[expr_df$sample_type %in% c("Tumor", "Normal"), ],
  FUN = function(x) c(n = length(x), mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE))
)

expr_summary_fixed <- data.frame(
  cancer_code = expr_summary$cancer_code,
  sample_type = expr_summary$sample_type,
  n = expr_summary$STAB1[, "n"],
  mean_STAB1 = expr_summary$STAB1[, "mean"],
  median_STAB1 = expr_summary$STAB1[, "median"]
)

write.table(
  expr_summary_fixed,
  "results/STAB1_simple_tumor_normal_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Analysis complete.")
message("Created files:")
message("results/STAB1_TCGA_tumor_survival_dataset.tsv")
message("results/STAB1_pan_cancer_OS_cox.tsv")
message("results/STAB1_cancer_wise_OS_cox.tsv")
message("results/STAB1_simple_tumor_normal_summary.tsv")
message("results/STAB1_pan_cancer_OS_KM.pdf")
