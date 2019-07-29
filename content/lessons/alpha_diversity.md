Resources
---------

This draws from [phyloseq plot\_bar tutorial](https://joey711.github.io/phyloseq/plot_bar-examples.html).

Data
----

This tutorial uses the 10% Atacama subset data (note that for the demux and dada2 tutorial we used the 1% Atacama subset)

Getting ready
=============

First we load libraries.

``` r
library(readr)
library(phyloseq)
```

``` r
atacama.ps = read_rds(atacama.rds)
print(atacama.ps)
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3388 taxa and 68 samples ]
    ## sample_data() Sample Data:       [ 68 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3388 taxa by 7 taxonomic ranks ]

Visualize alpha-diversity
=========================

``` r
plot_richness(atacama.ps, x="TransectName", 
              measures=c("Observed", "Shannon", "Chao1"), color="TransectName") + theme_bw()
```

![](alpha_diversity_plots_files/figure-markdown_github/unnamed-chunk-2-1.png)

Alpha-Diversity Boxplots
------------------------

It is a bit hard to compare the two different transects because many points are overlapping, let's add a boxplot layer so we can compare the distribution of alpha-diversity values between the transects.

``` r
plot_richness(atacama.ps, x="TransectName", 
              measures=c("Observed", "Shannon", "Chao1"), color="TransectName") + 
              geom_boxplot() +
              geom_point() +
              theme_bw() 
```

![](alpha_diversity_plots_files/figure-markdown_github/unnamed-chunk-3-1.png)

Alpha-Diversity as a function of other parameters
-------------------------------------------------

It also might be interesting to explore whether other parameters have an effect on alpha-diversity

``` r
sample_variables(atacama.ps)
```

    ##  [1] "BarcodeSequence"                 "LinkerPrimerSequence"           
    ##  [3] "Elevation"                       "ExtractConcen"                  
    ##  [5] "AmpliconConcentration"           "ExtractGroupNo"                 
    ##  [7] "TransectName"                    "SiteName"                       
    ##  [9] "Depth"                           "pH"                             
    ## [11] "TOC"                             "EC"                             
    ## [13] "AverageSoilRelativeHumidity"     "RelativeHumiditySoilHigh"       
    ## [15] "RelativeHumiditySoilLow"         "PercentRelativeHumiditySoil_100"
    ## [17] "AverageSoilTemperature"          "TemperatureSoilHigh"            
    ## [19] "TemperatureSoilLow"              "Vegetation"                     
    ## [21] "PercentCover"                    "Description"

### Elevation

``` r
plot_richness(atacama.ps, x="Elevation", 
              measures=c("Observed", "Shannon", "Chao1"), color="TransectName") + theme_bw()
```

![](alpha_diversity_plots_files/figure-markdown_github/unnamed-chunk-5-1.png)

### Depth

``` r
plot_richness(atacama.ps, x="Depth", 
              measures=c("Observed", "Shannon", "Chao1"), color="TransectName") + theme_bw()
```

![](alpha_diversity_plots_files/figure-markdown_github/unnamed-chunk-6-1.png)

Session Info
============

Always print `sessionInfo` for reproducibility!

``` r
sessionInfo()
```

    ## R version 3.4.1 (2017-06-30)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/openblas-base/libblas.so.3
    ## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] datasets  utils     stats     grDevices graphics  methods   base     
    ## 
    ## other attached packages:
    ##  [1] phyloseq_1.22.3    readr_1.1.1        GGally_1.3.2      
    ##  [4] broom_0.4.3        openintro_1.7.1    rvest_0.3.2       
    ##  [7] xml2_1.1.1         stringr_1.2.0      lubridate_1.6.0   
    ## [10] googlesheets_0.2.2 ggplot2_2.2.1      rmarkdown_1.8.6   
    ## [13] knitr_1.18         downloader_0.4    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Biobase_2.38.0      httr_1.3.1          tidyr_0.7.2        
    ##  [4] jsonlite_1.5        splines_3.4.1       foreach_1.4.3      
    ##  [7] assertthat_0.2.0    stats4_3.4.1        cellranger_1.1.0   
    ## [10] yaml_2.1.16         pillar_1.0.1        backports_1.1.2    
    ## [13] lattice_0.20-35     glue_1.2.0          digest_0.6.13      
    ## [16] RColorBrewer_1.1-2  XVector_0.18.0      colorspace_1.3-2   
    ## [19] htmltools_0.3.6     Matrix_1.2-10       plyr_1.8.4         
    ## [22] psych_1.7.8         pkgconfig_2.0.1     zlibbioc_1.24.0    
    ## [25] purrr_0.2.4         scales_0.5.0        tibble_1.4.1       
    ## [28] mgcv_1.8-18         IRanges_2.12.0      BiocGenerics_0.24.0
    ## [31] lazyeval_0.2.0      mnormt_1.5-5        survival_2.41-3    
    ## [34] magrittr_1.5        evaluate_0.10.1     nlme_3.1-131       
    ## [37] MASS_7.3-47         foreign_0.8-69      vegan_2.5-3        
    ## [40] tools_3.4.1         data.table_1.10.4-2 hms_0.3            
    ## [43] S4Vectors_0.16.0    munsell_0.4.3       cluster_2.0.6      
    ## [46] bindrcpp_0.2        Biostrings_2.46.0   ade4_1.7-13        
    ## [49] compiler_3.4.1      rlang_0.1.6         rhdf5_2.22.0       
    ## [52] grid_3.4.1          RCurl_1.95-4.8      iterators_1.0.8    
    ## [55] biomformat_1.6.0    igraph_1.1.2        labeling_0.3       
    ## [58] bitops_1.0-6        gtable_0.2.0        codetools_0.2-15   
    ## [61] multtest_2.34.0     reshape_0.8.7       reshape2_1.4.3     
    ## [64] R6_2.2.2            dplyr_0.7.4         bindr_0.1          
    ## [67] rprojroot_1.3-2     permute_0.9-4       ape_5.2            
    ## [70] stringi_1.1.6       parallel_3.4.1      Rcpp_0.12.14
