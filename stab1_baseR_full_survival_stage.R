options(stringsAsFactors = FALSE)

dir.create("results/full", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

message("Starting full STAB1 base-R survival and stage-adjusted analysis...")

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The R package 'survival' is missing.")
}

library(survival)

expr_file <- "data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
surv_file <- "data/Survival_SupplementalTable_S1_20171025_xena_sp"
stab1_file <- "data/STAB1_expression_row.tsv"

if (!file.exists(expr_file)) stop("Missing expression file: ", expr_file)
if (!file.exists(surv_file)) stop("Missing survival file: ", surv_file)

if (!file.exists(stab1_file) || file.info(stab1_file)$size == 0) {
  message("Extracting STAB1 row...")
  cmd <- paste(
    "zcat", shQuote(expr_file),
    "| awk 'NR==1 || $1==\"STAB1\"' >",
    shQuote(stab1_file)
  )
  status <- system(cmd)
  if (status != 0) stop("Failed to extract STAB1 row.")
}

message("Reading STAB1 expression...")
stab1 <- read.delim(stab1_file, check.names = FALSE)

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

message("Reading survival and clinical file...")
surv <- read.delim(surv_file, check.names = FALSE)

clean_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("STAGE IV", x)] <- "Stage IV"
  out[grepl("STAGE III", x)] <- "Stage III"
  out[grepl("STAGE II", x)] <- "Stage II"
  out[grepl("STAGE I", x) & is.na(out)] <- "Stage I"
  factor(out, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
}

surv_df <- data.frame(
  patient = substr(surv[["sample"]], 1, 12),
  sample = surv[["sample"]],
  cancer_type = surv[["cancer type abbreviation"]],
  stage_raw = surv[["ajcc_pathologic_tumor_stage"]],
  stage_simple = clean_stage(surv[["ajcc_pathologic_tumor_stage"]]),
  OS_event = as.numeric(surv[["OS"]]),
  OS_time = as.numeric(surv[["OS.time"]]),
  DSS_event = as.numeric(surv[["DSS"]]),
  DSS_time = as.numeric(surv[["DSS.time"]]),
  DFI_event = as.numeric(surv[["DFI"]]),
  DFI_time = as.numeric(surv[["DFI.time"]]),
  PFI_event = as.numeric(surv[["PFI"]]),
  PFI_time = as.numeric(surv[["PFI.time"]]),
  stringsAsFactors = FALSE
)

message("Merging expression and clinical data...")

tumor_expr <- expr_df[expr_df$sample_type == "Tumor", ]

main <- merge(
  tumor_expr,
  surv_df,
  by = "patient",
  all.x = FALSE,
  all.y = FALSE
)

write.table(
  main,
  "results/full/STAB1_merged_expression_clinical_all_endpoints.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Merged samples: ", nrow(main))

run_model <- function(dat, endpoint, cancer_label, stage_adjust = FALSE) {
  time_col <- paste0(endpoint, "_time")
  event_col <- paste0(endpoint, "_event")
  
  if (!time_col %in% names(dat) || !event_col %in% names(dat)) return(NULL)
  
  d <- dat[
    !is.na(dat$STAB1) &
      !is.na(dat[[time_col]]) &
      !is.na(dat[[event_col]]) &
      dat[[time_col]] > 0,
  ]
  
  if (nrow(d) < 40) return(NULL)
  if (sum(d[[event_col]], na.rm = TRUE) < 10) return(NULL)
  
  med <- median(d$STAB1, na.rm = TRUE)
  d$STAB1_group <- ifelse(d$STAB1 >= med, "High", "Low")
  d$STAB1_group <- factor(d$STAB1_group, levels = c("Low", "High"))
  
  if (length(unique(d$STAB1_group)) < 2) return(NULL)
  
  if (stage_adjust) {
    d <- d[!is.na(d$stage_simple), ]
    d$stage_simple <- droplevels(d$stage_simple)
    
    if (nrow(d) < 40) return(NULL)
    if (sum(d[[event_col]], na.rm = TRUE) < 10) return(NULL)
    if (length(unique(d$stage_simple)) < 2) return(NULL)
    
    form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ STAB1_group + stage_simple"))
    model_type <- "Stage-adjusted Cox"
  } else {
    form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ STAB1_group"))
    model_type <- "Univariate Cox"
  }
  
  cx <- try(coxph(form, data = d), silent = TRUE)
  if (inherits(cx, "try-error")) return(NULL)
  
  s <- summary(cx)
  
  coef_name <- "STAB1_groupHigh"
  if (!coef_name %in% rownames(s$coefficients)) return(NULL)
  
  logrank_p <- NA_real_
  if (!stage_adjust) {
    lr <- try(survdiff(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ STAB1_group")), data = d), silent = TRUE)
    if (!inherits(lr, "try-error")) {
      logrank_p <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
    }
  }
  
  data.frame(
    cancer_type = cancer_label,
    endpoint = endpoint,
    model = model_type,
    n = nrow(d),
    events = sum(d[[event_col]], na.rm = TRUE),
    median_cutoff_STAB1 = med,
    HR_high_vs_low = s$coefficients[coef_name, "exp(coef)"],
    lower95 = s$conf.int[coef_name, "lower .95"],
    upper95 = s$conf.int[coef_name, "upper .95"],
    cox_p_value = s$coefficients[coef_name, "Pr(>|z|)"],
    logrank_p_value = logrank_p,
    stringsAsFactors = FALSE
  )
}

add_fdr <- function(tab) {
  if (is.null(tab) || nrow(tab) == 0) return(tab)
  tab$FDR_p_value <- NA_real_
  groups <- unique(paste(tab$endpoint, tab$model, sep = "___"))
  for (g in groups) {
    idx <- paste(tab$endpoint, tab$model, sep = "___") == g
    tab$FDR_p_value[idx] <- p.adjust(tab$cox_p_value[idx], method = "BH")
  }
  tab
}

endpoints <- c("OS", "DSS", "DFI", "PFI")

message("Running pan-cancer models...")
pan_results <- list()

for (ep in endpoints) {
  pan_results[[paste(ep, "univ", sep = "_")]] <- run_model(main, ep, "Pan-cancer", stage_adjust = FALSE)
  pan_results[[paste(ep, "stage", sep = "_")]] <- run_model(main, ep, "Pan-cancer", stage_adjust = TRUE)
}

pan_table <- do.call(rbind, pan_results)
pan_table <- add_fdr(pan_table)

write.table(
  pan_table,
  "results/full/STAB1_pan_cancer_all_endpoints_cox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Running cancer-wise models...")
cancer_results <- list()
cancers <- sort(unique(main$cancer_type))

for (ct in cancers) {
  dat_ct <- main[main$cancer_type == ct, ]
  for (ep in endpoints) {
    cancer_results[[paste(ct, ep, "univ", sep = "_")]] <- run_model(dat_ct, ep, ct, stage_adjust = FALSE)
    cancer_results[[paste(ct, ep, "stage", sep = "_")]] <- run_model(dat_ct, ep, ct, stage_adjust = TRUE)
  }
}

cancer_table <- do.call(rbind, cancer_results)
cancer_table <- add_fdr(cancer_table)

write.table(
  cancer_table,
  "results/full/STAB1_cancer_wise_all_endpoints_cox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Creating key cancer summary...")
key_cancers <- c("KIRC", "GBM", "PRAD", "ACC", "THYM", "UVM", "LUSC", "SKCM", "OV")

key_table <- cancer_table[cancer_table$cancer_type %in% key_cancers, ]
key_table <- key_table[order(key_table$cancer_type, key_table$endpoint, key_table$model), ]

write.table(
  key_table,
  "results/full/STAB1_key_cancers_all_endpoints_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Running TCGA tumor-normal expression comparison where available...")

patient_cancer <- unique(surv_df[, c("patient", "cancer_type")])

expr_annot <- merge(
  expr_df,
  patient_cancer,
  by = "patient",
  all.x = FALSE,
  all.y = FALSE
)

expr_annot <- expr_annot[expr_annot$sample_type %in% c("Tumor", "Normal"), ]

expr_results <- list()

for (ct in sort(unique(expr_annot$cancer_type))) {
  d <- expr_annot[expr_annot$cancer_type == ct, ]
  if (!all(c("Tumor", "Normal") %in% unique(d$sample_type))) next
  
  n_tumor <- sum(d$sample_type == "Tumor")
  n_normal <- sum(d$sample_type == "Normal")
  
  if (n_tumor < 5 || n_normal < 5) next
  
  wt <- try(wilcox.test(STAB1 ~ sample_type, data = d), silent = TRUE)
  if (inherits(wt, "try-error")) next
  
  med_tumor <- median(d$STAB1[d$sample_type == "Tumor"], na.rm = TRUE)
  med_normal <- median(d$STAB1[d$sample_type == "Normal"], na.rm = TRUE)
  
  direction <- ifelse(med_tumor > med_normal, "Higher in tumor", "Lower in tumor")
  
  expr_results[[ct]] <- data.frame(
    cancer_type = ct,
    n_tumor = n_tumor,
    n_normal = n_normal,
    median_tumor = med_tumor,
    median_normal = med_normal,
    direction = direction,
    wilcox_p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
}

expr_table <- do.call(rbind, expr_results)

if (!is.null(expr_table) && nrow(expr_table) > 0) {
  expr_table$FDR_p_value <- p.adjust(expr_table$wilcox_p_value, method = "BH")
}

write.table(
  expr_table,
  "results/full/STAB1_TCGA_tumor_vs_normal_wilcox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Full analysis complete.")
message("Created:")
message("results/full/STAB1_pan_cancer_all_endpoints_cox.tsv")
message("results/full/STAB1_cancer_wise_all_endpoints_cox.tsv")
message("results/full/STAB1_key_cancers_all_endpoints_summary.tsv")
message("results/full/STAB1_TCGA_tumor_vs_normal_wilcox.tsv")
