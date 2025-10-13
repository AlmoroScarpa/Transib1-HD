1_RNA_prep_expression
================
Matthew Beaumont
2025-10-13

``` bash
knitr::opts_chunk$set(echo = TRUE)
```

# Visualisation

We first need to extract the expression data from our GSNAP aligned
reads.

``` bash
#!/bin/bash

fai="/Volumes/Data/Tools/RefGenomes/dmel/dna/dmel_consensus_TEs/dmel_allchrom_TEs_2021_newTEs.fasta.fai"
pyscript="/Volumes/Data/Projects/dmelR2_p-ele/scripts/mRNA-coverage-senseantisense.py"
indir="/Volumes/Temp/almoRNA/RNA/mapped/bam/100MM"
outdir="/Volumes/Temp/almoRNA/RNA/results/expression"
seqs="Transib_Riccardo, hobo, P-element"

mkdir -p "$outdir"

for sample in TEN1 TEN2 TEN3 TOM1 TOM2 TOM3 TENfTOMm1 TENfTOMm2 TENfTOMm3 TOMfTENm1 TOMfTENm2; do
    echo "Processing $sample"
    samtools view "$indir/${sample}.sort.bam" | \
        python "$pyscript" --sam - --sample-id "$sample" --seqs "$seqs" --fai "$fai" -a \
        > "$outdir/${sample}.txt"
done

cat "$outdir"/*.txt | perl -pe 's/-/mrna_data/g' > "$outdir/expr.forr"

# Tissue column for 'wf' needed adding back.
#awk 'BEGIN{OFS="\mrna_data"}{$2 = $2 "\mrna_data" "wf"; print}' expr.forr > expr_wf.forr
```

Then again for splicing.

``` bash
#!/bin/bash

pyscript="/Volumes/Data/Projects/dmelR2_p-ele/scripts/mRNA-splicing-senseantisense.py"
indir="/Volumes/Temp/almoRNA/RNA/mapped/bam/100MM"
outdir="/Volumes/Temp/almoRNA/RNA/results/splicing"
seqs="Transib_Riccardo, hobo, P-element"

# All transcripts are 'RA' variants unless stated otherwise.

mkdir -p "$outdir"

for sample in TEN1 TEN2 TEN3 TOM1 TOM2 TOM3 TENfTOMm1 TENfTOMm2 TENfTOMm3 TOMfTENm1 TOMfTENm2; do
    echo "Processing $sample"
    samtools view "$indir/${sample}.sort.bam" | \
        python "$pyscript" --sam - --sample-id "$sample" --seqs "$seqs" \
        > "$outdir/${sample}.txt"
done

cat "$outdir"/*.txt | perl -pe 's/-/mrna_data/g' > "$outdir/spli.forr"
```

And then we visualised both the P-element expression and the splicing of
its’ three introns together in ggplot2.

``` r
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite)

theme_set(theme_bw())

#target <- "Transib_Riccardo"
sminfreq <- 0.1

h <- read.table("/Volumes/Temp/almoRNA/RNA/results/expr/expr.forr")
names(h) <- c("sample_id","strand","te","pos","cov")

h[,3] <- gsub("mrna_data", "-", h[,3])

#a <- subset(h, te %in% c("Transib_Riccardo", "hobo", "P-element"))
a <- subset(h, te %in% c("P-element"))

a <- a %>%
  mutate(sample_id = as.character(sample_id)) %>%
  separate(sample_id, into = c("sample_type","rep"), sep = "(?<=[A-Za-z])(?=[0-9])")

a$key<-paste0(a$sample_type, "_",a$te, "_", a$rep,"_",a$pos,"_",a$strand)
 
s  <- subset(a, strand == "se")
as <- subset(a, strand == "ase")

spli <- read.table("/Volumes/Temp/almoRNA/RNA/results/spli/spli.forr")
names(spli) <- c("sample_id","strand","te","skey","start","end","rawcount","freq")

spli[,4] <- gsub("mrna_data", "-", spli[,4])
spli[,3] <- gsub("mrna_data", "-", spli[,3])

aspli <- subset(spli, te %in% c("Transib_Riccardo", "hobo", "P-element") & freq > sminfreq)

aspli <- aspli %>%
  mutate(sample_id = as.character(sample_id)) %>%
  separate(sample_id, into = c("sample_type","rep"), sep = "(?<=[A-Za-z])(?=[0-9])")

aspli$start <- aspli$start - 1

aspli$keystart<-paste0(aspli$sample_type, "_", aspli$te, "_", aspli$rep,"_",aspli$start,"_",aspli$strand)
aspli$keyend<-paste0(aspli$sample_type, "_", aspli$te, "_", aspli$rep,"_",aspli$end,"_",aspli$strand)

aspli<-merge(x=aspli,y=a[,c("key","cov")],by.x="keystart",by.y="key") 
aspli<-merge(x=aspli,y=a[,c("key","cov")],by.x="keyend",by.y="key")

aspli$size <- log(aspli$freq + 1)

a_s  <- subset(aspli, strand == "se")
a_as <- subset(aspli, strand == "ase")

a_s$size_scaled  <- a_s$size - 1
a_as$size_scaled <- a_as$size - 1

expr_spli_plot <- ggplot() +
  geom_polygon(data = s,  mapping = aes(x = pos, y = cov),     fill = 'grey',      color = 'grey') +
  geom_polygon(data = as, mapping = aes(x = pos, y = -cov),    fill = 'lightgrey', color = 'lightgrey') +
  geom_curve(data = a_s,  mapping = aes(x = start, y = cov.x,  xend = end, yend = cov.y, linewidth = size_scaled), curvature = -0.15, ncp = 10, show.legend = FALSE) +
  geom_curve(data = a_as, mapping = aes(x = start, y = -cov.x, xend = end, yend = -cov.y, linewidth = size_scaled), curvature = 0.15,  ncp = 10, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(~ te + sample_type + rep, ncol = 3) +
  scale_linewidth(range = c(0.1, 2)) +
  xlab("position") +
  ylab("expression [rpm]") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black"),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.4, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

expr_spli_plot
```

![](RNA_expression_files/figure-gfm/rna_expr-1.png)<!-- -->

``` r
ggsave("figs/pele_expr_spli.eps", expr_spli_plot, width = 20, height = 8)
ggsave("figs/pele_expr_spli.svg", expr_spli_plot, width = 20, height = 8)
ggsave("figs/pele_expr_spli.png", expr_spli_plot, width = 20, height = 8, dpi = 600)
```

``` r
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite)

theme_set(theme_bw())

sminfreq <- 0.1
region_start <- 5900861
region_end   <- 5905399

h <- read.table("/Volumes/Temp/almoRNA/RNA/results/expr_x/expr_act5C.forr")
names(h) <- c("sample_id","strand","te","pos","cov")

a  <- h  %>% filter(pos >= region_start & pos <= region_end)

a[,3] <- gsub("mrna_data", "-", a[,3])

a <- subset(h, te %in% c("X"))

a <- a %>%
  mutate(sample_id = as.character(sample_id)) %>%
  separate(sample_id, into = c("sample_type","rep"), sep = "(?<=[A-Za-z])(?=[0-9])")

a$key<-paste0(a$sample_type, "_",a$te, "_", a$rep,"_",a$pos,"_",a$strand)

s  <- subset(a, strand == "se")
as <- subset(a, strand == "ase")

spli <- read.table("/Volumes/Temp/almoRNA/RNA/results/spli_x/spli_act5C.forr")
names(spli) <- c("sample_id","strand","te","skey","start","end","rawcount","freq")

aspli <- aspli %>%
  filter((start >= region_start & start <= region_end) |
         (end   >= region_start & end   <= region_end))

spli[,4] <- gsub("mrna_data", "-", spli[,4])
spli[,3] <- gsub("mrna_data", "-", spli[,3])

aspli <- subset(spli, te %in% c("X") & freq > sminfreq)

aspli <- aspli %>%
  mutate(sample_id = as.character(sample_id)) %>%
  separate(sample_id, into = c("sample_type","rep"), sep = "(?<=[A-Za-z])(?=[0-9])")

aspli$start <- aspli$start - 1

aspli$keystart<-paste0(aspli$sample_type, "_", aspli$te, "_", aspli$rep,"_",aspli$start,"_",aspli$strand)
aspli$keyend<-paste0(aspli$sample_type, "_", aspli$te, "_", aspli$rep,"_",aspli$end,"_",aspli$strand)

aspli<-merge(x=aspli,y=a[,c("key","cov")],by.x="keystart",by.y="key") 
aspli<-merge(x=aspli,y=a[,c("key","cov")],by.x="keyend",by.y="key")

aspli$size <- log(aspli$freq + 1)

a_s  <- subset(aspli, strand == "se")
a_as <- subset(aspli, strand == "ase")

a_s$size_scaled  <- a_s$size - 1
a_as$size_scaled <- a_as$size - 1

expr_spli_plot <- ggplot() +
                  geom_polygon(data = s,  mapping = aes(x = pos, y = cov),     fill = 'grey',      color = 'grey') +
                  geom_polygon(data = as, mapping = aes(x = pos, y = -cov),    fill = 'lightgrey', color = 'lightgrey') +
                  geom_curve(data = a_s,  mapping = aes(x = start, y = cov.x,  xend = end, yend = cov.y, linewidth = size_scaled), curvature = -0.15, ncp = 10, show.legend = FALSE) +
                  geom_curve(data = a_as, mapping = aes(x = start, y = -cov.x, xend = end, yend = -cov.y, linewidth = size_scaled), curvature = 0.15,  ncp = 10, show.legend = FALSE) +
                  geom_hline(yintercept = 0, color = "black") +
                  facet_wrap(~ te + sample_type + rep, ncol = 3) +
                  scale_linewidth(range = c(0.1, 2)) +
                  xlab("position") +
                  ylab("expression [rpm]") +
                  xlim(5900861, 5905399) +
                  theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(color = "black", fill = NA, size = 1),
                    strip.background = element_rect(fill = "grey80", color = "black"),
                    strip.text = element_text(size = 10),
                    panel.spacing = unit(0.4, "lines"),
                    legend.position = "bottom",
                    legend.direction = "horizontal"
  )

expr_spli_plot
```

![](RNA_expression_files/figure-gfm/act5c_expr-1.png)<!-- -->

``` r
ggsave("figs/ACT5C_expr_spli.svg", expr_spli_plot, width = 20, height = 8)
ggsave("figs/ACT5C_expr_spli.png", expr_spli_plot, width = 20, height = 8, dpi = 600)
ggsave("figs/ACT5C_expr_spli.eps", expr_spli_plot, width = 20, height = 8)
```
