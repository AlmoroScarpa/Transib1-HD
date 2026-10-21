Transcriptome Analysis
================
2026-10-20

\#RNA analysis

Index fo STAR:

    STAR --runThreadN 10 --runMode genomeGenerate --genomeDir Genome_indices --genomeFastaFiles TE_lib.fa

Mapping (STAR):

    STAR --runThreadN 30 --genomeDir Genome_indices --readFilesIn /path/to/sample1.gz 
    /path/to/sample2.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix --outSAMstrandField intronMotif 

Running the extraction scripts:

    for i in *.bam;do base=$(basename "$i" .sorted.bam);
    samtools view "$i"|python mRNA-coverage-senseantisense.py --sam - --sample-id "$base" --seqs Transib_Riccardo --fai TE_lib.fa.fai -a > "$base"_cov.txt;
    done

    cat *._cov.txt > all_cov.txt

    for i in *.bam;do base=$(basename "$i" .sorted.bam);
    samtools view "$i"|python mRNA-splicing-senseantisense.py --sam - --sample-id "$base" --seqs Transib_Riccardo  > "$base"_as.txt;
    done

    cat *_as.txt > all_as.txt

\#Visualization of transcript abundance and splice sites (Transib 1 as
example - same code for hobo and p-element)

``` r
library(ggpubr)
```

    ## Loading required package: ggplot2

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ lubridate 1.9.5     ✔ tibble    3.3.1
    ## ✔ purrr     1.2.2     ✔ tidyr     1.3.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
all_cov<-read_tsv(file="/home/anna/transib_ana/Transib1-HD/Documentation_v2.0/Data/transib_merged_cov.txt",
                   col_names =c("sample_id","strand","te","pos","cov"))%>%
  mutate(sample_id=gsub("m", "♂", sample_id),
         sample_id=gsub("f","♀",sample_id))%>%
  separate(sample_id,sep="_",into = c("sample","replicate"),remove = F)%>%
  mutate(key=paste0(sample,"_",te,"_",replicate,"_",pos,"_",strand))
```

    ## Rows: 72768 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): sample_id, strand, te
    ## dbl (2): pos, cov
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
sense<-all_cov%>%filter(strand=="se")
antisense<-all_cov%>%filter(strand=="ase")

all_spli<-read_tsv(file = "/home/anna/transib_ana/Transib1-HD/Documentation_v2.0/Data/transib_merged_as.txt",
                    col_names=c("sample_id","strand","te","skey","start","end","rawcount","freq"))%>%
  mutate(sample_id=gsub("m", "♂", sample_id),
         sample_id=gsub("f","♀",sample_id))%>%
  separate(sample_id,sep="_",into = c("sample","replicate"),remove = F)%>%
  mutate(start=start-1,
         keystart=paste0(sample,"_",te,"_",replicate,"_",start,"_",strand),
         keyend=paste0(sample,"_",te,"_",replicate,"_",end,"_",strand))
```

    ## Rows: 36 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): sample_id, strand, te, skey
    ## dbl (4): start, end, rawcount, freq
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
spli_merg1<-merge(x=all_spli,y=all_cov[,c("key","cov")],by.x="keystart",by.y="key")
spli_merg2<-merge(x=spli_merg1,y=all_cov[,c("key","cov")],by.x="keyend",by.y="key")
  
sense_spli<-spli_merg2%>%filter(strand=="se")
asense_spli<-spli_merg2%>%filter(strand=="ase")

RNA_all<-ggplot()+
  geom_polygon(data=sense,mapping=aes(x = pos, y = cov),fill="grey",color="grey")+
  geom_polygon(data=antisense,mapping=aes(x = pos, y = -cov),fill="lightgrey",color="lightgrey")+
  geom_curve(data=sense_spli,mapping=aes(x = start, y = cov.x,  xend = end, yend = cov.y,alpha = freq,linewidth = freq), curvature = -0.15, ncp = 10, show.legend = FALSE)+
  geom_curve(data=asense_spli,mapping=aes(x = start, y = -cov.x,  xend = end, yend = -cov.y,alpha=freq,linewidth = freq),curvature = 0.15, ncp = 10, show.legend = FALSE)+
  theme_bw()+
  facet_wrap(~ sample+replicate,ncol=3)+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size=13),
        title = element_text(size = 18),
        strip.text.x = element_text(size = 16))+
  labs(y="expression [rpm]", x="position")

#ggsave(RNA_all,file="/home/anna/transib_ana/RNA_reanalysis/pele_all_exp_plots.png",
#       width=16,height=10)
#ggsave(RNA_all,file="/home/anna/transib_ana/RNA_reanalysis/pele_all_exp_plots.pdf",
#       width=16,height=10,device = cairo_pdf)
RNA_all
```

![](RNA_anadoc_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Stringtie:

    for i in *.bam; do base=$(basename "$i" .sorted.bam);stringtie "$i" -o stringtie/"$base".gtf

Transdecoder:

    for i in stringtie/*.gtf;do ~/TransDecoder-TransDecoder-v6.0.0/TransDecoder --genome ../DmelTEs_Bergmann_2021_edit.fa --gtf "$i" -O trns_dcdr/"$base"/ ;done

Plotting the abundance of the different transcripts:

``` r
tframe<-read_tsv("/home/anna/transib_ana/Transib1-HD/Documentation_v2.0/Data/transcript_summary.tsv")%>%
  mutate(strain=gsub("f","♀",strain))%>%
  mutate(strain=gsub("m","♂",strain))
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

    ## Rows: 36 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): ID, strain
    ## dbl (2): FPKM, scptfrm
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
by_Tform<-ggplot(tframe,aes(x=strain,y=FPKM))+
  geom_point()+
  stat_compare_means(method = "t.test",
                                comparisons = list(c("TEN♀TOM♂","TOM♀TEN♂"),
                                                    c("TEN♀TOM♂","TOM")))+
  facet_wrap(.~scptfrm)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  labs(y="FPKM")+
  theme_bw()

by_Tform
```

![](RNA_anadoc_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#ggsave(by_Tform,file="/home/anna/transib_ana/RNA_reanalysis/star_new/stringtie/byTform_plot.pdf",
#       height = 5,width=8, ,device = cairo_pdf)
```
