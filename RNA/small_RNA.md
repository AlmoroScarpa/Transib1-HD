small_RNA
================
Matthew Beaumont
2025-10-13

# piRNA distribution

``` r
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
 
pirna_dist <- read.table("/Volumes/Temp/almoRNA/sRNA/results/pirna_dist/all/all_pirna_dist_fixed.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(pirna_dist) <- c("SampleID", "TE", "strand", "pos", "count")

# "p-element_dna/p" "hobo_dna/hat" "transib_riccardo"

pirna_dist <- pirna_dist %>%
  mutate(
    TE = as.factor(TE)
  )

pirna_dist_fil <- pirna_dist %>%
  filter(TE == "hobo_dna/hat")

pirna_dist_fil <- pirna_dist_fil %>%
  mutate(count = ifelse(strand == "-", -count, count))

max_count <- max(abs(pirna_dist_fil$count))

p <- ggplot(pirna_dist_fil) +
  geom_segment(aes(x = pos, xend = pos, y = 0, yend = count), color = "dimgrey") +
  geom_hline(yintercept = 0, color = "black") +
  scale_x_continuous(breaks = c(1, 1000, 2000, 3000),) +
  labs(
    title = NULL,
    x = "position",
    y = "counts (per million piRNAs)"
  ) +
  facet_wrap(~ SampleID, scales = "fixed", ncol = 4) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black"),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "none"
  )

p
```

![](small_RNA_files/figure-gfm/pirna_dist-1.png)<!-- -->

``` r
file <- "hobo_pirna_dist"

ggsave(paste0("figs/", file, ".png"), plot = p, width = 8, height = 4, dpi = 400)
ggsave(paste0("figs/", file, ".svg"), plot = p, width = 8, height = 4)
ggsave(paste0("figs/", file, ".eps"), plot = p, width = 8, height = 4)
```

# Pingpong

``` r
library(ggplot2)
library(tidyr)
library(dplyr)

pingpong <- "/Volumes/Temp/almoRNA/sRNA/results/pingpong/all/all_pingpong_fixed.txt"
pingpong <- read.table(pingpong, header = TRUE, sep = "\t")
colnames(pingpong) <- c("SampleID", "TE", "Strand", "Position", "Count", "Normalised", "Zscore")

# "p-element_dna/p" "hobo_dna/hat" "transib_riccardo"

pingpong_fil <- pingpong %>%
  filter(TE == "hobo_dna/hat")

pingpong_fil <- pingpong_fil %>%
  filter(Strand == "+")

pp <- ggplot(pingpong_fil, aes(x = Position, y = Normalised,
                           fill = Position == 10)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "red3")) +
  scale_y_continuous(breaks = c(0, 0.5), limits = c(0, 0.7)) +
  scale_x_continuous(breaks = c(1, 10, 19), limits = c(1, 19)) +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(~ SampleID, scales = "fixed", ncol = 4) +
  labs(title = NULL,
       x = "position",
       y = "normalised count") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black"),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "none"
  )

pp
```

![](small_RNA_files/figure-gfm/pingpong-1.png)<!-- -->

``` r
file <- "hobo_pingpong"

ggsave(paste0("figs/", file, ".png"), plot = pp, width = 8, height = 4, dpi = 400)
ggsave(paste0("figs/", file, ".svg"), plot = pp, width = 8, height = 4)
ggsave(paste0("figs/", file, ".eps"), plot = pp, width = 8, height = 4)
```

# Pingpong comp

``` r
library(tidyverse)
library(ggrepel)

pp <- read.table("/Volumes/Temp/almoRNA/sRNA/results/pingpong/all/all_pingpong_fixed.txt", header = FALSE)
colnames(pp) <- c("SampleID", "TE", "Strand", "Position", "Count", "Normalised", "Zscore")

pp10 <- pp %>%
  filter(Position == 10) %>%
  filter(Strand == "+") %>%
  mutate(TE = sub("_.*", "", TE)) %>%
  mutate(TE = case_when(
    TE == "transib" ~ "transib1 (consensus)",
    TE == "transib1" ~ "transib1 (degraded)",
    TRUE ~ TE
  ))

pp_mean <- pp10 %>%
  filter(SampleID %in% c("TEN015", "TOM008")) %>%
  group_by(SampleID, TE) %>%
  summarise(Zscore = mean(Zscore), .groups = "drop") %>%
  pivot_wider(
    id_cols = TE,
    names_from = SampleID,
    values_from = Zscore
  ) %>%
  filter(!is.na(TEN015) & !is.na(TOM008)) %>%
  mutate(
    diff = pmax(TEN015, TOM008) - pmin(TEN015, TOM008),
    outlier = diff > 3
  )

sp <- ggplot(pp_mean, aes(x = TEN015, y = TOM008)) +
  geom_point(size = 1, color = "steelblue4") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red3") +
  geom_text_repel(
    data = subset(pp_mean, outlier),
    aes(label = TE),
    size = 4,
    color = "black",
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  labs(
    x = "TEN015 ping-pong Z-score (pos10)",
    y = "TOM008 ping-pong Z-score (pos10)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )

sp
```

![](small_RNA_files/figure-gfm/pingpong_zscores-1.png)<!-- -->

``` r
ggsave("figs/pingpong_zscores.png", plot = sp, width = 8, height = 4, dpi = 400)
ggsave("figs/pingpong_zscores.svg", plot = sp, width = 8, height = 4)
ggsave("figs/pingpong_zscores.eps", plot = sp, width = 8, height = 4)
```

# piRNA abundance

``` r
library(ggplot2)
library(dplyr)
library(readr)

pa <- read.table("/Volumes/Temp/almoRNA/sRNA/results/pirna_comp/all/all_pirna_abu_fixed.txt", header = FALSE, sep = "\t")
colnames(pa) <- c("SampleID", "TE", "RNA", "Count", "Normalised")

# "p-element_dna/p" "hobo_dna/hat" "transib_riccardo"

pa_fil <- pa %>%
  filter(TE == "hobo_dna/hat") %>%
  filter(RNA == "pirna")

pab <- ggplot(pa_fil, aes(x = SampleID, y = Normalised)) +
      geom_bar(stat = "identity", color = "black", fill = "grey") +
      theme_bw() +
      labs(
        title = "Hobo",
        x = NULL,
        y = "piRNA abundance (per million piRNAs)"
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

pab
```

![](small_RNA_files/figure-gfm/pirna_abundance-1.png)<!-- -->

``` r
ggsave("figs/hobo_pirna_abundance.png", plot = pab, width = 8, height = 4, dpi = 400)
ggsave("figs/hobo_pirna_abundance.svg", plot = pab, width = 8, height = 4)
ggsave("figs/hobo_pirna_abundance.eps", plot = pab, width = 8, height = 4)
```

# TE specific piRNA comparison

``` r
library(tidyverse)
library(ggrepel)

pa <- read.table("/Volumes/Temp/almoRNA/sRNA/results/pirna_comp/all/all_pirna_abu_fixed.txt", header = FALSE, sep = "\t")
colnames(pa) <- c("SampleID", "TE", "RNA", "Count", "Normalised")

# "p-element_dna/p" "hobo_dna/hat" "transib_riccardo"

pa_fil <- pa %>%
  mutate(TE = sub("_.*", "", TE)) %>%
  filter(RNA == "pirna", SampleID %in% c("TEN015", "TOM008")) %>%
  distinct(SampleID, TE, .keep_all = TRUE) %>%
  pivot_wider(
    id_cols = TE,
    names_from = SampleID,
    values_from = Normalised,
    values_fill = 0
  ) %>%
  filter(TEN015 > 0 & TOM008 > 0) %>%
  mutate(
    TE = case_when(
      TE == "transib" ~ "transib1 (consensus)",
      TE == "transib1" ~ "transib1 (degraded)",
      TE == "p-element" ~ "P-element",
      TE == "hobo" ~ "Hobo",
      TRUE ~ TE
    )
  ) %>%
  mutate(
    ratio = pmax(TEN015, TOM008) / pmin(TEN015, TOM008),
    outlier = ratio > 10,
    label_this = outlier | TE %in% c("P-element", "Hobo")
  )

label_data <- pa_fil %>% filter(label_this)

sp <- ggplot(pa_fil, aes(x = `TEN015`, y = `TOM008`)) +
  geom_point(aes(color = label_this), size = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick3")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_text_repel(
    data = label_data,
    aes(label = TE, color = label_this),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_x_log10(
    limits = c(1, max(pa_fil$TENfxTOMm) * 1.1),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_number()
  ) +
  scale_y_log10(
    limits = c(1, max(pa_fil$TOMfxTENm) * 1.1),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_number()
  ) +
  labs(
    x = "TOMmTENf piRNAs (per million piRNAs)",
    y = "TENmTOMf piRNAs (per million piRNAs)",
    color = "Labeled TE"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    legend.position = "none"  # remove legend if not needed
  )

sp
```

![](small_RNA_files/figure-gfm/pirna_comp-1.png)<!-- -->

``` r
ggsave("figs/pirna_comp.png", plot = sp, width = 8, height = 4, dpi = 400)
ggsave("figs/pirna_comp.svg", plot = sp, width = 8, height = 4)
ggsave("figs/pirna_comp.eps", plot = sp, width = 8, height = 4)
```

``` r
library(tidyverse)
library(ggrepel)
library(scales) 

pa <- read.table("/Volumes/Temp/almoRNA/sRNA/results/pirna_comp/all/all_pirna_abu_fixed.txt",
                 header = FALSE, sep = "\t")
colnames(pa) <- c("SampleID", "TE", "RNA", "Count", "Normalised")

# "p-element_dna/p" "hobo_dna/hat" "transib_riccardo"

pa <- pa %>%
  mutate(TE = sub("_.*", "", TE)) %>%
  mutate(TE = case_when(
    TE == "transib" ~ "transib1 (consensus)",
    TE == "transib1" ~ "transib1 (degraded)",
    TE == "p-element" ~ "P-element",
    TE == "hobo" ~ "Hobo",
    TRUE ~ TE
  ))

pa_fil <- pa %>%
  filter(RNA == "pirna") 

pa_fil <- pa_fil %>%
  mutate(SampleGroup = sub("R[123]$", "", SampleID))

pa_mean <- pa_fil %>%
  group_by(SampleGroup, TE) %>%
  summarise(Normalised = mean(Normalised), .groups = "drop") %>%
  pivot_wider(
    id_cols = TE,
    names_from = SampleGroup,
    values_from = Normalised,
    values_fill = 0
  ) %>%
  filter(`TENfxTOMm` > 0 & `TOMfxTENm` > 0) %>%
  mutate(
    ratio = pmax(`TENfxTOMm`, `TOMfxTENm`) / pmin(`TENfxTOMm`, `TOMfxTENm`),
    outlier = ratio > 5,
    label_this = outlier | TE %in% c("P-element", "Hobo")
  )
label_data <- pa_mean %>% filter(label_this)

sp <- ggplot(pa_mean, aes(x = `TENfxTOMm`, y = `TOMfxTENm`)) +
  geom_point(aes(color = label_this), size = 1) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick3")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_text_repel(
    data = label_data,
    aes(label = TE, color = label_this),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_x_log10(
    limits = c(1, max(pa_mean$TENfxTOMm) * 1.1),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_number()
  ) +
  scale_y_log10(
    limits = c(1, max(pa_mean$TOMfxTENm) * 1.1),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_number()
  ) +
  labs(
    x = "TOMmTENf piRNAs (per million piRNAs)",
    y = "TENmTOMf piRNAs (per million piRNAs)",
    color = "Labeled TE"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    legend.position = "none"  # remove legend if not needed
  )

sp
```

![](small_RNA_files/figure-gfm/pirna_comp_rc-1.png)<!-- -->

``` r
ggsave("figs/pirna_comp_rc.png", plot = sp, width = 8, height = 4, dpi = 400)
ggsave("figs/pirna_comp_rc.svg", plot = sp, width = 8, height = 4)
ggsave("figs/pirna_comp_rc.eps", plot = sp, width = 8, height = 4)

#knitr::include_graphics("figs/pirna_comp_rc.png")
```
