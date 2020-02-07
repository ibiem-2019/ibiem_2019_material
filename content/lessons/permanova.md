    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-6

Data
----

This tutorial uses the 10% Atacama subset data (note that for the demux
and dada2 tutorial we used the 1% Atacama subset) There is a problem
when running `ordinate` on the 1% dataset. Not sure what it is!

``` r
# Directories

if(exists("params") && 
   !is.null(params[["atacama_ps_rds"]])){
  atacama.ps.rds=params[["atacama_ps_rds"]]
} else {
  atacama.ps.rds = "/data/tutorial_data/atacama_10pct.rds"
}
print(atacama.ps.rds)
```

    ## [1] "/data/tutorial_data/atacama_10pct.rds"

``` r
atacama.ps = read_rds(atacama.ps.rds)
print(atacama.ps)
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3508 taxa and 68 samples ]
    ## sample_data() Sample Data:       [ 68 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3508 taxa by 7 taxonomic ranks ]

Data Preprocessing
------------------

As with relative abundance plots and ordination, before performing a
permanova test we will want to prune rare taxa and transform the data.
We prune because we don’t want small differences in rare taxa to swamp
out major trends. The transformation is important because otherwise the
permanova could be affected by meaningless differences in total counts
between samples.

### Prune

The decision about how to prune is important, we need to think about
what we are throwing away, and how it might affect the analysis.

``` r
sample_min_count = 50

atacama.ps %>%
  prune_samples(sample_sums(.)>=sample_min_count, .) ->
  atacama.sample_prune

sample_sums(atacama.sample_prune) %>% sort
```

    ## YUN3259.1.1   YUN3346.2 YUN3259.1.3   BAQ2838.3   BAQ2838.2   BAQ2838.1 
    ##         149         950        1633        1942        2727        2945 
    ##   YUN3153.2 BAQ2420.1.1   BAQ2462.3   YUN1005.3 BAQ2420.1.2 BAQ2420.1.3 
    ##        3344        3353        3381        3564        4028        4042 
    ##   BAQ2420.3 YUN3533.1.3 YUN3259.1.2   BAQ2687.1   BAQ2420.2 BAQ4166.1.1 
    ##        4253        4267        4298        4325        4405        4505 
    ##   BAQ2462.2   YUN3153.3   YUN3259.2   BAQ3473.2   BAQ2687.2   YUN1242.1 
    ##        4574        4641        4648        4793        4946        4975 
    ##   YUN3346.1   YUN3346.3 YUN3533.1.2   YUN1609.1 YUN3533.1.1 YUN3856.1.1 
    ##        4986        4987        5015        5081        5190        5244 
    ##   BAQ2687.3 YUN3856.1.3   BAQ4697.2   YUN1242.3   BAQ3473.1   YUN3856.2 
    ##        5263        5270        5332        5422        5585        5638 
    ##   YUN3428.1 BAQ4166.1.3   YUN3428.3   BAQ4697.1   YUN3259.3   BAQ2462.1 
    ##        5649        5753        6064        6102        6198        6222 
    ## YUN3856.1.2   BAQ4166.2   BAQ4166.3 BAQ4166.1.2   YUN2029.2 YUN1005.1.1 
    ##        6354        6558        6838        6952        7232        7360 
    ##   BAQ3473.3   YUN3856.3   YUN3533.3   BAQ4697.3   YUN3533.2   YUN3428.2 
    ##        7617        7625        7743        7773        7971        8825

``` r
min_count = 3
min_sample = 2

prune.vec = filter_taxa(atacama.sample_prune, 
                       function(x) sum(x >= min_count) >= min_sample)
sum(prune.vec)
```

    ## [1] 1015

### Transform to even sampling depth.

Here we are performing the same fractional abundance transformation we
did before, then multiplying by 1 x 10^6 to convert those proprotions
back into whole numbers.

Pay attention to the y-axes in these plots of the raw counts, the pruned
counts, and the transformed counts.

``` r
atacama.even = transform_sample_counts(atacama.sample_prune, function(x) 1E6 * x/sum(x))

atacama.st_prune.even = prune_taxa(prune.vec, atacama.even)
ntaxa(atacama.st_prune.even)
```

    ## [1] 1015

``` r
plot_bar(atacama.st_prune.even)
```

![](permanova_files/figure-markdown_github/unnamed-chunk-4-1.png)

Permanova
---------

The ordination plot, suggest that there is a differences in microbial
communities between sites with and without vegetation. A permanova
analysis allows us to quantify this.

``` r
atacama.st_prune.even %>%
  sample_data %>%
  as("data.frame") ->
  atacama.st_prune.even.metadata

atacama.st_prune.even %>%
  distance(method="bray") ->
  atacama.st_prune.even.bray

adonis(atacama.st_prune.even.bray ~ Vegetation,
       data = atacama.st_prune.even.metadata) ->
  vegetation.adonis

print(vegetation.adonis)
```

    ## 
    ## Call:
    ## adonis(formula = atacama.st_prune.even.bray ~ Vegetation, data = atacama.st_prune.even.metadata) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##            Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
    ## Vegetation  1    1.8272 1.82720  4.5895 0.0811  0.001 ***
    ## Residuals  52   20.7024 0.39812         0.9189           
    ## Total      53   22.5296                 1.0000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The permanova comparing with and without vegetation finds a
**significant** difference with a p-value of 0.001. However, the
permanova test is sensitive to differences in dispersion between the two
groups, so we need to test for homogeneity of dispersion to examine
whether differences in dispersion could lead us to falsely conclude that
vegetation results in a signifcant differences between the groups.

Homogeneity of Dispersion
-------------------------

``` r
betadisper(atacama.st_prune.even.bray,
           atacama.st_prune.even.metadata[["Vegetation"]]) ->
  vegetation.beta

permutest(vegetation.beta) ->
  vegetation.beta.permute

print(vegetation.beta.permute)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
    ## Groups     1 0.006827 0.0068271 2.7681    999  0.088 .
    ## Residuals 52 0.128251 0.0024664                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The test for homogeneity of dispersions is not significant (p-value of
0.088). In other words, there is not a significant difference in
dispersion between the groups, so a difference in dispersion between the
the two groups (with and without vegetation) is **not** the cause of the
significant permanova result. This provides strong evidence that
vegetation is associated with differences in microbial community.

Session Info
============

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] vegan_2.5-6     lattice_0.20-38 permute_0.9-5   phyloseq_1.30.0
    ##  [5] forcats_0.4.0   stringr_1.4.0   dplyr_0.8.3     purrr_0.3.3    
    ##  [9] readr_1.3.1     tidyr_1.0.0     tibble_2.1.3    ggplot2_3.2.1  
    ## [13] tidyverse_1.3.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-143        fs_1.3.1            lubridate_1.7.4    
    ##  [4] httr_1.4.1          tools_3.6.2         backports_1.1.5    
    ##  [7] R6_2.4.1            DBI_1.1.0           lazyeval_0.2.2     
    ## [10] BiocGenerics_0.32.0 mgcv_1.8-31         colorspace_1.4-1   
    ## [13] ade4_1.7-13         withr_2.1.2         tidyselect_0.2.5   
    ## [16] compiler_3.6.2      cli_2.0.1           rvest_0.3.5        
    ## [19] Biobase_2.46.0      xml2_1.2.2          labeling_0.3       
    ## [22] scales_1.1.0        digest_0.6.23       rmarkdown_2.1      
    ## [25] XVector_0.26.0      pkgconfig_2.0.3     htmltools_0.4.0    
    ## [28] dbplyr_1.4.2        rlang_0.4.2         readxl_1.3.1       
    ## [31] rstudioapi_0.10     generics_0.0.2      farver_2.0.3       
    ## [34] jsonlite_1.6        magrittr_1.5        biomformat_1.14.0  
    ## [37] Matrix_1.2-18       Rcpp_1.0.3          munsell_0.5.0      
    ## [40] S4Vectors_0.24.3    Rhdf5lib_1.8.0      fansi_0.4.1        
    ## [43] ape_5.3             lifecycle_0.1.0     stringi_1.4.5      
    ## [46] yaml_2.2.0          MASS_7.3-51.5       zlibbioc_1.32.0    
    ## [49] rhdf5_2.30.1        plyr_1.8.5          grid_3.6.2         
    ## [52] parallel_3.6.2      crayon_1.3.4        Biostrings_2.54.0  
    ## [55] haven_2.2.0         splines_3.6.2       multtest_2.42.0    
    ## [58] hms_0.5.3           zeallot_0.1.0       knitr_1.27         
    ## [61] pillar_1.4.3        igraph_1.2.4.2      reshape2_1.4.3     
    ## [64] codetools_0.2-16    stats4_3.6.2        reprex_0.3.0       
    ## [67] glue_1.3.1          evaluate_0.14       data.table_1.12.8  
    ## [70] modelr_0.1.5        vctrs_0.2.1         foreach_1.4.7      
    ## [73] cellranger_1.1.0    gtable_0.3.0        assertthat_0.2.1   
    ## [76] xfun_0.12           broom_0.5.3         survival_3.1-8     
    ## [79] iterators_1.0.12    IRanges_2.20.2      cluster_2.1.0
