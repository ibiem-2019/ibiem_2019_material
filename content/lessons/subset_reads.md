``` r
library(tools)
library(tibble)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:GGally':
    ## 
    ##     nasa

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(readr)
```

    ## 
    ## Attaching package: 'readr'

    ## The following object is masked from 'package:rvest':
    ## 
    ##     guess_encoding

``` r
# Directories
data.dir = "/data/ibiem_2016_data"
subset.dir = file.path("/data/tutorial_data/ibiem2016_subset")

# make directory for output
dir.create(subset.dir, recursive = TRUE)
```

    ## Warning in dir.create(subset.dir, recursive = TRUE): '/data/tutorial_data/
    ## ibiem2016_subset' already exists

``` r
# Set variables for bash
Sys.setenv(DATA_DIR = data.dir)
Sys.setenv(SUBSET_DIR = subset.dir)
```

Check Data Integrity
--------------------

``` bash
cd $DATA_DIR
md5sum -c md5sum.txt
```

    ## ibiem_2017_map_v1.txt: OK
    ## ibiem_2017_map_v2.txt: OK
    ## ibiem_2017_map_v3_decoder.csv: OK
    ## ibiem_2017_map_v3.txt: OK
    ## Undetermined_S0_L001_I1_001.fastq.gz: OK
    ## Undetermined_S0_L001_R1_001.fastq.gz: OK
    ## Undetermined_S0_L001_R2_001.fastq.gz: OK

Generate Data Subset for Demo Purposes
--------------------------------------

``` bash
set -u
NUM_READS=20000
RANDSEED=1
for FASTQ_FULL in $DATA_DIR/*.gz ; do
  echo $FASTQ_FULL
  FASTQ_BASE=`basename $FASTQ_FULL`
  echo $FASTQ_BASE
  seqtk sample -s $RANDSEED $FASTQ_FULL $NUM_READS | gzip -c > $SUBSET_DIR/$FASTQ_BASE
  # zcat $SUBSET_DIR/$FASTQ_BASE | wc
done
```

    ## /data/ibiem_2016_data/Undetermined_S0_L001_I1_001.fastq.gz
    ## Undetermined_S0_L001_I1_001.fastq.gz
    ## /data/ibiem_2016_data/Undetermined_S0_L001_R1_001.fastq.gz
    ## Undetermined_S0_L001_R1_001.fastq.gz
    ## /data/ibiem_2016_data/Undetermined_S0_L001_R2_001.fastq.gz
    ## Undetermined_S0_L001_R2_001.fastq.gz

Copy metadata to subset directory
=================================

``` r
map.file = file.path(data.dir,"ibiem_2017_map_v3.txt")
decoder.file = file.path(data.dir,"ibiem_2017_map_v3_decoder.csv")

file.copy(map.file, subset.dir)
```

    ## [1] TRUE

``` r
file.copy(decoder.file, subset.dir)
```

    ## [1] TRUE

Generate md5sums for subset directory
=====================================

``` r
md5sum(files = list.files(subset.dir,full.names = TRUE)) %>%
  as.data.frame %>%
  rownames_to_column(var="path") %>%
  mutate(path=basename(path)) %>%
  select(".",path) %>%
  write_delim(file.path(subset.dir,"md5sum.txt"), col_names = FALSE)
```

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
    ## [1] tools     datasets  utils     stats     grDevices graphics  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2       readr_1.1.1        dplyr_0.7.4       
    ##  [4] tibble_1.4.1       GGally_1.3.2       broom_0.4.3       
    ##  [7] openintro_1.7.1    rvest_0.3.2        xml2_1.1.1        
    ## [10] stringr_1.2.0      lubridate_1.6.0    googlesheets_0.2.2
    ## [13] ggplot2_2.2.1      rmarkdown_1.8.6    knitr_1.18        
    ## [16] downloader_0.4    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.14       compiler_3.4.1     RColorBrewer_1.1-2
    ##  [4] cellranger_1.1.0   pillar_1.0.1       plyr_1.8.4        
    ##  [7] bindr_0.1          bitops_1.0-6       digest_0.6.13     
    ## [10] lattice_0.20-35    nlme_3.1-131       evaluate_0.10.1   
    ## [13] gtable_0.2.0       pkgconfig_2.0.1    rlang_0.1.6       
    ## [16] psych_1.7.8        yaml_2.1.16        parallel_3.4.1    
    ## [19] httr_1.3.1         hms_0.3            rprojroot_1.3-2   
    ## [22] grid_3.4.1         reshape_0.8.7      glue_1.2.0        
    ## [25] R6_2.2.2           foreign_0.8-69     reshape2_1.4.3    
    ## [28] tidyr_0.7.2        purrr_0.2.4        magrittr_1.5      
    ## [31] backports_1.1.2    scales_0.5.0       htmltools_0.3.6   
    ## [34] mnormt_1.5-5       assertthat_0.2.0   colorspace_1.3-2  
    ## [37] stringi_1.1.6      RCurl_1.95-4.8     lazyeval_0.2.0    
    ## [40] munsell_0.4.3
