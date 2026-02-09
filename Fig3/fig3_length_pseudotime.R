#!/usr/bin/env Rscript
# Fig.3: Length (Diameter) ~ Pseudotime smoothing plot + Pearson correlation
#
# Usage:
#   Rscript fig3_length_pseudotime_plot.R --config Fig3/configs/fig3_monocle3.yaml
#
# Expects Seurat object meta.data contains:
#   - pseudotime
#   - Length_px
#
suppressPackageStartupMessages({
  library(Seurat)
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    k <- args[[i]]
    if (!startsWith(k, "--")) stop("Unknown arg: ", k)
    key <- sub("^--", "", k)
    if (key %in% c("help")) {
      out[[key]] <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) stop("Missing value for ", k)
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

help_msg <- function() {
  cat(
"fig3_length_pseudotime_plot.R

Required:
  --config <yaml>     YAML config path (reuse Fig3 config)

Example:
  Rscript fig3_length_pseudotime_plot.R --config Fig3/configs/fig3_monocle3.yaml
", sep = "")
}

need_yaml <- function() {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Please run: install.packages('yaml')", call. = FALSE)
  }
}

write_text <- function(path, txt) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(txt, con = con)
}

args <- parse_args()
if (!is.null(args$help) && isTRUE(args$help)) { help_msg(); quit(status = 0) }
if (is.null(args$config)) { help_msg(); stop("Missing --config", call. = FALSE) }

need_yaml()
cfg <- yaml::read_yaml(args$config)

# Inputs / outputs
in_rds <- if (!is.null(cfg$out_obj_rds)) cfg$out_obj_rds else if (!is.null(cfg$input_rds)) cfg$input_rds else "obj_oo.rds"
out_dir <- if (!is.null(cfg$out_dir)) cfg$out_dir else "results/Fig3"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- if (!is.null(cfg$plot_out_pdf)) cfg$plot_out_pdf else file.path(out_dir, "length_vs_pseudotime.pdf")
out_png <- if (!is.null(cfg$plot_out_png)) cfg$plot_out_png else file.path(out_dir, "length_vs_pseudotime.png")
out_cor_txt <- if (!is.null(cfg$plot_out_cor_txt)) cfg$plot_out_cor_txt else file.path(out_dir, "length_vs_pseudotime_correlation.txt")

span <- if (!is.null(cfg$plot_loess_span)) as.numeric(cfg$plot_loess_span) else 0.5
line_width <- if (!is.null(cfg$plot_linewidth)) as.numeric(cfg$plot_linewidth) else 1
alpha_fill <- if (!is.null(cfg$plot_alpha)) as.numeric(cfg$plot_alpha) else 0.15
color_hex <- if (!is.null(cfg$plot_color)) as.character(cfg$plot_color) else "#e64b35"

x_label <- if (!is.null(cfg$plot_xlabel)) as.character(cfg$plot_xlabel) else "Pseudotime"
y_label <- if (!is.null(cfg$plot_ylabel)) as.character(cfg$plot_ylabel) else "Diameter(µm)"

cor_method <- if (!is.null(cfg$plot_cor_method)) as.character(cfg$plot_cor_method) else "pearson"

message("[fig3-plot] Loading: ", in_rds)
obj <- readRDS(in_rds)

if (is.null(obj$pseudotime)) stop("Missing meta.data column: pseudotime", call. = FALSE)
if (is.null(obj$Length_px)) stop("Missing meta.data column: Length_px", call. = FALSE)

df <- obj@meta.data[, c("pseudotime", "Length_px")]
df <- df[!is.na(df$pseudotime) & !is.na(df$Length_px), , drop = FALSE]
if (nrow(df) == 0) stop("No valid rows after removing NA pseudotime/Length_px", call. = FALSE)

# Correlation (requested)
cor_value <- suppressWarnings(cor(df$pseudotime, df$Length_px, method = cor_method))
cor_value_round <- round(cor_value, 3)
message("[fig3-plot] ", cor_method, " correlation (pseudotime vs Length_px) = ", cor_value_round)

write_text(out_cor_txt, c(
  paste0("method: ", cor_method),
  paste0("correlation: ", cor_value_round),
  paste0("n: ", nrow(df))
))
message("[fig3-plot] Saved correlation: ", out_cor_txt)

df.long <- reshape2::melt(df, id.vars = "pseudotime", variable.name = "Variable", value.name = "Value")

df.long <- df.long %>%
  group_by(Variable) %>%
  mutate(ScaledValue = scale(Value)[, 1]) %>%
  ungroup()

p <- ggplot(df.long, aes(x = pseudotime, y = Value, colour = Variable, fill = Variable)) +
  stat_smooth(method = "loess", span = span, se = TRUE, linewidth = line_width, alpha = alpha_fill) +
  scale_colour_manual(values = color_hex) +
  scale_fill_manual(values = color_hex) +
  theme_classic() +
  labs(x = x_label, y = y_label) +
  guides(colour = "none", fill = "none")

ggsave(out_pdf, plot = p, width = 5, height = 4, units = "in")
ggsave(out_png, plot = p, width = 5, height = 4, units = "in", dpi = 300)

message("[fig3-plot] Saved: ", out_pdf)
message("[fig3-plot] Saved: ", out_png)
message("[fig3-plot] Done.")
