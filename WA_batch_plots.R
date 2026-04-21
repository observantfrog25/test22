library(here)
library(dplyr)
library(ggplot2)
library(tidyr)

chunk_dir <- "results/chunks/WA_para_testing/"

chunk_files <- list.files(chunk_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(chunk_files) == 0) stop("No chunk files found in ", chunk_dir)

all_metrics <- lapply(chunk_files, function(f) {
  chunk <- readRDS(f)
  chunk$metrics
})

combined <- dplyr::bind_rows(all_metrics)

long_df <- combined %>%
  pivot_longer(
    cols      = c(mean_aadt, median_aadt, min_aadt),
    names_to  = "summary_method",
    values_to = "aadt_exposure"
  ) %>%
  mutate(
    summary_method = recode(summary_method,
      mean_aadt   = "Mean",
      median_aadt = "Median",
      min_aadt    = "Min"
    )
  )

p_hist <- ggplot(long_df, aes(x = aadt_exposure, fill = summary_method)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~summary_method, scales = "free_y", ncol = 1) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Distribution of AADT exposure across WA destination cells",
    x     = "AADT exposure (sum of per-origin path summaries, scaled 1e-6)",
    y     = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

p_density <- ggplot(long_df, aes(x = aadt_exposure, color = summary_method)) +
  geom_density(linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Density comparison: Mean vs Median vs Min AADT exposure",
    x     = "AADT exposure",
    y     = "Density",
    color = "Summary"
  ) +
  theme_minimal()

p_box <- ggplot(long_df, aes(x = summary_method, y = aadt_exposure,
                             fill = summary_method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Boxplot comparison: Mean vs Median vs Min AADT exposure",
    x     = "",
    y     = "AADT exposure"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

if (!dir.exists("Plots/AADT examples/")) {
  dir.create("Plots/AADT examples/", recursive = TRUE)
}

ggsave("Plots/AADT examples/WA_aadt_histograms.png", p_hist,
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave("Plots/AADT examples/WA_aadt_density.png", p_density,
       width = 10, height = 6, dpi = 300, bg = "white")
ggsave("Plots/AADT examples/WA_aadt_boxplot.png", p_box,
       width = 8, height = 6, dpi = 300, bg = "white")

message("Plots saved to Plots/AADT examples/")

summary_stats <- long_df %>%
  group_by(summary_method) %>%
  summarise(
    n        = n(),
    mean     = mean(aadt_exposure, na.rm = TRUE),
    median   = median(aadt_exposure, na.rm = TRUE),
    sd       = sd(aadt_exposure, na.rm = TRUE),
    min      = min(aadt_exposure, na.rm = TRUE),
    max      = max(aadt_exposure, na.rm = TRUE),
    .groups  = "drop"
  )

print(summary_stats)
