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

Simple absolute abundance plot
------------------------------

Let's start by making a simple absolute abundance bar plot with samples on the x-axis and absolute abundance on the y-axis.

``` r
plot_bar(atacama.ps)
```

![](absolute_abundance_plots_files/figure-markdown_github/unnamed-chunk-2-1.png)

### Absolute abundance plot with kingdom

Now let's look add the kingdom distribution of the reads. This is a good sanity check for contamination - we expect the *vast* majority of our reads to be bacterial.

First let's look up the names of the taxonomic levels

``` r
rank_names(atacama.ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

``` r
plot_bar(atacama.ps, fill="Kingdom")
```

![](absolute_abundance_plots_files/figure-markdown_github/unnamed-chunk-4-1.png)

### Getting rid of annoying black lines

You may have noticed that the bars are divided by short horizontal black lines. This is because each bar actually consists of a stack of sub-bars, each representing a different OTU. Usually the sub-bars delineation are more confusing than helpful, so I prefer to get rid of it. This isn't that hard to do! by default `plot_bar` outlines each sub-bar in black, but if we outline sub-bars with the same color as the bar fill, we get rid of the delineation. Let's re-make the *Absolute abundance plot with kingdom* from above, but without the black sub-bar outlines.

``` r
plot_bar(atacama.ps, fill="Kingdom") +
  geom_bar(aes(color=Kingdom, fill=Kingdom), stat="identity", position="stack")
```

![](absolute_abundance_plots_files/figure-markdown_github/unnamed-chunk-5-1.png)

### Absolute abundance plot with phyllum

Now do the same as above, but with coloring by phyllum instead of by kingdom.

``` r
plot_bar(atacama.ps, fill="Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
```

![](absolute_abundance_plots_files/figure-markdown_github/unnamed-chunk-6-1.png)

Facet Abundance Plots
---------------------

From our plots so far, we see some variation in number of reads per sample. This is to be expected, but we want to be sure that there are no confounding problems, for example, categories of samples that are outliers. To this end, we can make facet plots to group samples by different metadata features. Let's remake the kingdom absolute abundance plot so that samples are grouped by a variable, this will give us an idea if there is a correlation between the variable and abundance.

First let's see what variables are available.

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

``` r
get_variable(atacama.ps,"SiteName") %>%
  unique()
```

    ##  [1] "BAQ1370" "BAQ1552" "BAQ2420" "BAQ2462" "BAQ2687" "BAQ2838" "BAQ3473"
    ##  [8] "BAQ4166" "BAQ4697" "BAQ895"  "YUN1005" "YUN1242" "YUN1609" "YUN2029"
    ## [15] "YUN3008" "YUN3153" "YUN3184" "YUN3259" "YUN3346" "YUN3428" "YUN3533"
    ## [22] "YUN3856"

That's a lot of different values, it's going to make an illegible plot. Let's try another one.

``` r
get_variable(atacama.ps,"TransectName")  %>%
  unique()
```

    ## [1] "Baquedano" "Yungay"

Two different values seems better

``` r
plot_bar(atacama.ps, facet_grid=~TransectName, fill="Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
```

![](absolute_abundance_plots_files/figure-markdown_github/unnamed-chunk-10-1.png)

You might notice something weird things about the plot: There are large spaces between samples without bars. Let's try to fix this.

### "Missing" Values in Facet Plots

By default, the facets of a facet plot use the same scale. Often this is the right thing to do - if the y-axis of each of the abundance plots had a different scale, it would be hard to compare bars between the facets. But sometimes it is not, such as our case. For the TransectName facet plot, since the x-axis is samples, each sample has a place in each facet, even though TransectName is mutually exclusive - each sample corresponds to a specific site and depth combination, and each site occurs in only one of the transects. Let's fix the TransectName facet plot so that each facet only shows the samples that correspond to that Transect.

``` r
plot_bar(atacama.ps, facet_grid=~TransectName, fill="Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  facet_grid(~TransectName, scale="free_x")
```

![](absolute_abundance_plots_files/figure-markdown_github/unnamed-chunk-11-1.png)

Once we fix that, you may notice that the facets default to being the same size, despite the number of bars in the facet. Let's adjust the facet plot so the bars are all the same width, by making the size of the facet proportional to the number bars in it by using the `space` argument to `facet_grid`

``` r
plot_bar(atacama.ps, facet_grid=~TransectName, fill="Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  facet_grid(~TransectName, scale="free_x", space = "free_x")
```

![](absolute_abundance_plots_files/figure-markdown_github/unnamed-chunk-12-1.png)

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
