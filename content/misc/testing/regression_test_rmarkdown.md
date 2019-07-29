``` r
library(here)
```

    ## here() starts at /home/guest/IBIEM_2018_2019

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(fs)
library(rmarkdown)
library(magrittr)

# Load the following for sessioninfo
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 2.2.1     ✔ readr   1.3.1
    ## ✔ tibble  2.0.1     ✔ purrr   0.3.1
    ## ✔ tidyr   0.8.0     ✔ stringr 1.3.0
    ## ✔ ggplot2 2.2.1     ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ tidyr::extract()   masks magrittr::extract()
    ## ✖ dplyr::filter()    masks stats::filter()
    ## ✖ dplyr::lag()       masks stats::lag()
    ## ✖ purrr::set_names() masks magrittr::set_names()

``` r
regression_output_dir = here("regression_scratch")
if (dir_exists(regression_output_dir)) {dir_delete(regression_output_dir)}
dir.create(regression_output_dir, recursive = TRUE)
# tmp_scratch = "/tmp/scratch/atacama_1pct"
# list.files(tmp_scratch)
# dir_delete(tmp_scratch)
Sys.setenv(REGRESSION_OUT = regression_output_dir)
```

``` r
regression_rmds = c(
                    "/home/guest/IBIEM_2018_2019/content/lessons/dada2_tutorial_1_6.Rmd")
```

Apt-get Packages
================

sra-toolkit
-----------

``` bash
fastq-dump -X 5 -Z SRR390728
```

    ## 2019-03-13T21:25:26 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:26 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:26 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:27 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:27 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:27 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:27 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:27 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:27 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## 2019-03-13T21:25:27 fastq-dump.2.8.2 warn: block-size in local file 32768 does not match requested value 131072
    ## Read 5 spots for SRR390728
    ## Written 5 spots for SRR390728
    ## @SRR390728.1 1 length=72
    ## CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
    ## +SRR390728.1 1 length=72
    ## ;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
    ## @SRR390728.2 2 length=72
    ## AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
    ## +SRR390728.2 2 length=72
    ## ;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
    ## @SRR390728.3 3 length=72
    ## CCAGCCTGGCCAACAGAGTGTTACCCCGTTTTTACTTATTTATTATTATTATTTTGAGACAGAGCATTGGTC
    ## +SRR390728.3 3 length=72
    ## -;;;8;;;;;;;,*;;';-4,44;,:&,1,4'./&19;;;;;;669;;99;;;;;-;3;2;0;+;7442&2/
    ## @SRR390728.4 4 length=72
    ## ATAAAATCAGGGGTGTTGGAGATGGGATGCCTATTTCTGCACACCTTGGCCTCCCAAATTGCTGGGATTACA
    ## +SRR390728.4 4 length=72
    ## 1;;;;;;,;;4;3;38;8%&,,;)*;1;;,)/%4+,;1;;);;;;;;;4;(;1;;;;24;;;;41-444//0
    ## @SRR390728.5 5 length=72
    ## TTAAGAAATTTTTGCTCAAACCATGCCCTAAAGGGTTCTGTAATAAATAGGGCTGGGAAAACTGGCAAGCCA
    ## +SRR390728.5 5 length=72
    ## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;9445552;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;446662

Pip Packages
============

\#\# QIIME 1
------------

title: "Demultiplex" output: md\_document: variant: markdown\_github html\_document: df\_print: paged ---

### Paths, Directories, and Shell Variables

To keep the code readable and portable, it is nice to assign paths to variables. We also need to use the R `Sys.setenv` command to make shell variables that can be used in the bash chunks below.

``` r
# Directories
data.dir = "/data/tutorial_data/atacama_1pct"
output.dir = regression_output_dir
demux.dir = file.path(output.dir, "demux")

# make directory for output
if (dir_exists(demux.dir)) {dir_delete(demux.dir)}
dir_create(demux.dir)

# Files
map.file = file.path(data.dir,"sample_metadata.tsv")
barcode.fastq = file.path(data.dir,"barcodes.fastq.gz")
r1.fastq = file.path(data.dir,"forward.fastq.gz")
r2.fastq = file.path(data.dir,"reverse.fastq.gz")

# Set variables for bash
Sys.setenv(MAP_FILE = map.file)
Sys.setenv(OUT_DIR = output.dir)
Sys.setenv(DEMUX_DIR = demux.dir)
Sys.setenv(RAW_FASTQ_DIR = data.dir)
Sys.setenv(BARCODE_FASTQ = barcode.fastq)
```

### Demux R1 and R2

``` bash
set -u
for CURREAD in "forward" "reverse"
do
   CURREAD_DIR=$DEMUX_DIR/${CURREAD}
   TAGDIR=$CURREAD_DIR/tagged
    split_libraries_fastq.py -r 999 -n 999 -q 0 -p 0.0001 \
        --sequence_read_fps $RAW_FASTQ_DIR/${CURREAD}.fastq.gz \
        --output_dir $TAGDIR \
        --barcode_read_fps $BARCODE_FASTQ \
        --mapping_fps $MAP_FILE \
        --phred_offset 33 \
        --barcode_type golay_12 \
        --rev_comp_mapping_barcodes \
        --store_demultiplexed_fastq 
        
    split_sequence_file_on_sample_ids.py --input_seqs_fp $TAGDIR/seqs.fastq \
                     --file_type fastq \
                     --output_dir $CURREAD_DIR
                     
    rm -rf $TAGDIR
done
```

### Rename and move split FASTQs

``` r
for (curread in c("forward","reverse")) {
  curpath = file.path(demux.dir,curread)
  for (fastq_path in list.files(curpath, full.names = TRUE,pattern = ".fastq")){
    new_path = path_ext_remove(fastq_path)
    new_path = path_file(new_path)
    new_path = path(demux.dir, new_path, ext=paste0(curread,".fastq"))
    file_move(fastq_path, new_path)
  }
}
```

Bioconductor Packages
=====================

DADA2
-----

----------------------------------
==================================

----------------------------------
==================================

\#----------------------------------
------------------------------------

``` r
library(dada2)
# library(readr)
# library(stringr)
# library(dplyr)
# library(tibble)
# library(magrittr)
# library(ggplot2)
# library(fs)
```

``` r
# output.dir = regression_output_dir
demux.dir = file.path(output.dir, "demux")
scratch.dir = file.path(output.dir, "dada2")

data.dir = "/data/tutorial_data/atacama_1pct"
map.file = file.path(data.dir,"sample_metadata.tsv")

if (dir_exists(scratch.dir)) {
  dir_delete(scratch.dir)
}
dir_create(scratch.dir)

ps.rds = file.path(scratch.dir, "atacama_1pct.rds")

silva.ref = "/data/references/dada/silva_nr_v128_train_set.fa.gz"
silva.species.ref = "/data/references/dada/silva_species_assignment_v128.fa.gz"
```

### Filter and Trim

First we read in the names of the fastq files, and perform some string manipulation to get lists of the forward and reverse fastq files in matched order:

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(demux.dir, pattern="forward.fastq", full.names = TRUE))
fnRs <- sort(list.files(demux.dir, pattern="reverse.fastq", full.names = TRUE))

sample.names = fnFs %>% 
  basename %>%
  str_replace(".forward.fastq","") 
```

### Perform filtering and trimming

``` r
filt_path <- file.path(scratch.dir, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
```

``` r
filt.out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=10, truncLen=c(145,140),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(filt.out)
```

    ##                           reads.in reads.out
    ## BAQ1370.1.3.forward.fastq        1         0
    ## BAQ1552.1.1.forward.fastq        1         1
    ## BAQ2420.1.1.forward.fastq      721       672
    ## BAQ2420.1.2.forward.fastq      640       601
    ## BAQ2420.1.3.forward.fastq      656       633
    ## BAQ2420.2.forward.fastq        611       567

### Learn the Error Rates

``` r
filtFs = filtFs[file_exists(filtFs)]
filtRs = filtRs[file_exists(filtRs)]
```

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## Initializing error rates to maximum possible estimate.
    ## Sample 1 - 1 reads in 1 unique sequences.
    ## Sample 2 - 672 reads in 566 unique sequences.
    ## Sample 3 - 601 reads in 469 unique sequences.
    ## Sample 4 - 633 reads in 468 unique sequences.
    ## Sample 5 - 567 reads in 394 unique sequences.
    ## Sample 6 - 607 reads in 397 unique sequences.
    ## Sample 7 - 805 reads in 519 unique sequences.
    ## Sample 8 - 547 reads in 360 unique sequences.
    ## Sample 9 - 434 reads in 280 unique sequences.
    ## Sample 10 - 820 reads in 662 unique sequences.
    ## Sample 11 - 607 reads in 365 unique sequences.
    ## Sample 12 - 785 reads in 555 unique sequences.
    ## Sample 13 - 542 reads in 432 unique sequences.
    ## Sample 14 - 449 reads in 380 unique sequences.
    ## Sample 15 - 331 reads in 254 unique sequences.
    ## Sample 16 - 1020 reads in 785 unique sequences.
    ## Sample 17 - 790 reads in 655 unique sequences.
    ## Sample 18 - 1121 reads in 740 unique sequences.
    ## Sample 19 - 965 reads in 744 unique sequences.
    ## Sample 20 - 1203 reads in 905 unique sequences.
    ## Sample 21 - 990 reads in 718 unique sequences.
    ## Sample 22 - 1195 reads in 877 unique sequences.
    ## Sample 23 - 1144 reads in 888 unique sequences.
    ## Sample 24 - 783 reads in 413 unique sequences.
    ## Sample 25 - 723 reads in 428 unique sequences.
    ## Sample 26 - 1018 reads in 552 unique sequences.
    ## Sample 27 - 819 reads in 397 unique sequences.
    ## Sample 28 - 387 reads in 169 unique sequences.
    ## Sample 29 - 503 reads in 229 unique sequences.
    ## Sample 30 - 1 reads in 1 unique sequences.
    ## Sample 31 - 688 reads in 306 unique sequences.
    ## Sample 32 - 656 reads in 277 unique sequences.
    ## Sample 33 - 6 reads in 6 unique sequences.
    ## Sample 34 - 792 reads in 394 unique sequences.
    ## Sample 35 - 1 reads in 1 unique sequences.
    ## Sample 36 - 1 reads in 1 unique sequences.
    ## Sample 37 - 1 reads in 1 unique sequences.
    ## Sample 38 - 415 reads in 270 unique sequences.
    ## Sample 39 - 533 reads in 290 unique sequences.
    ## Sample 40 - 1 reads in 1 unique sequences.
    ## Sample 41 - 53 reads in 52 unique sequences.
    ## Sample 42 - 645 reads in 434 unique sequences.
    ## Sample 43 - 251 reads in 180 unique sequences.
    ## Sample 44 - 901 reads in 751 unique sequences.
    ## Sample 45 - 986 reads in 714 unique sequences.
    ## Sample 46 - 606 reads in 379 unique sequences.
    ## Sample 47 - 166 reads in 136 unique sequences.
    ## Sample 48 - 720 reads in 495 unique sequences.
    ## Sample 49 - 970 reads in 723 unique sequences.
    ## Sample 50 - 1278 reads in 934 unique sequences.
    ## Sample 51 - 951 reads in 739 unique sequences.
    ## Sample 52 - 892 reads in 645 unique sequences.
    ## Sample 53 - 846 reads in 681 unique sequences.
    ## Sample 54 - 732 reads in 548 unique sequences.
    ## Sample 55 - 1093 reads in 674 unique sequences.
    ## Sample 56 - 1190 reads in 839 unique sequences.
    ## Sample 57 - 872 reads in 665 unique sequences.
    ## Sample 58 - 917 reads in 598 unique sequences.
    ## Sample 59 - 612 reads in 404 unique sequences.
    ## Sample 60 - 921 reads in 657 unique sequences.
    ## Sample 61 - 1166 reads in 829 unique sequences.
    ##    selfConsist step 2 
    ##    selfConsist step 3 
    ##    selfConsist step 4 
    ##    selfConsist step 5 
    ## Convergence after  5  rounds.
    ## Total reads used:  40925

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## Initializing error rates to maximum possible estimate.
    ## Sample 1 - 1 reads in 1 unique sequences.
    ## Sample 2 - 672 reads in 592 unique sequences.
    ## Sample 3 - 601 reads in 519 unique sequences.
    ## Sample 4 - 633 reads in 527 unique sequences.
    ## Sample 5 - 567 reads in 444 unique sequences.
    ## Sample 6 - 607 reads in 473 unique sequences.
    ## Sample 7 - 805 reads in 624 unique sequences.
    ## Sample 8 - 547 reads in 441 unique sequences.
    ## Sample 9 - 434 reads in 359 unique sequences.
    ## Sample 10 - 820 reads in 715 unique sequences.
    ## Sample 11 - 607 reads in 463 unique sequences.
    ## Sample 12 - 785 reads in 657 unique sequences.
    ## Sample 13 - 542 reads in 473 unique sequences.
    ## Sample 14 - 449 reads in 404 unique sequences.
    ## Sample 15 - 331 reads in 284 unique sequences.
    ## Sample 16 - 1020 reads in 893 unique sequences.
    ## Sample 17 - 790 reads in 711 unique sequences.
    ## Sample 18 - 1121 reads in 905 unique sequences.
    ## Sample 19 - 965 reads in 839 unique sequences.
    ## Sample 20 - 1203 reads in 1066 unique sequences.
    ## Sample 21 - 990 reads in 818 unique sequences.
    ## Sample 22 - 1195 reads in 1021 unique sequences.
    ## Sample 23 - 1144 reads in 985 unique sequences.
    ## Sample 24 - 783 reads in 526 unique sequences.
    ## Sample 25 - 723 reads in 528 unique sequences.
    ## Sample 26 - 1018 reads in 708 unique sequences.
    ## Sample 27 - 819 reads in 552 unique sequences.
    ## Sample 28 - 387 reads in 244 unique sequences.
    ## Sample 29 - 503 reads in 319 unique sequences.
    ## Sample 30 - 1 reads in 1 unique sequences.
    ## Sample 31 - 688 reads in 436 unique sequences.
    ## Sample 32 - 656 reads in 432 unique sequences.
    ## Sample 33 - 6 reads in 6 unique sequences.
    ## Sample 34 - 792 reads in 536 unique sequences.
    ## Sample 35 - 1 reads in 1 unique sequences.
    ## Sample 36 - 1 reads in 1 unique sequences.
    ## Sample 37 - 1 reads in 1 unique sequences.
    ## Sample 38 - 415 reads in 322 unique sequences.
    ## Sample 39 - 533 reads in 377 unique sequences.
    ## Sample 40 - 1 reads in 1 unique sequences.
    ## Sample 41 - 53 reads in 52 unique sequences.
    ## Sample 42 - 645 reads in 505 unique sequences.
    ## Sample 43 - 251 reads in 214 unique sequences.
    ## Sample 44 - 901 reads in 816 unique sequences.
    ## Sample 45 - 986 reads in 801 unique sequences.
    ## Sample 46 - 606 reads in 446 unique sequences.
    ## Sample 47 - 166 reads in 154 unique sequences.
    ## Sample 48 - 720 reads in 572 unique sequences.
    ## Sample 49 - 970 reads in 811 unique sequences.
    ## Sample 50 - 1278 reads in 1048 unique sequences.
    ## Sample 51 - 951 reads in 809 unique sequences.
    ## Sample 52 - 892 reads in 737 unique sequences.
    ## Sample 53 - 846 reads in 751 unique sequences.
    ## Sample 54 - 732 reads in 609 unique sequences.
    ## Sample 55 - 1093 reads in 762 unique sequences.
    ## Sample 56 - 1190 reads in 946 unique sequences.
    ## Sample 57 - 872 reads in 755 unique sequences.
    ## Sample 58 - 917 reads in 683 unique sequences.
    ## Sample 59 - 612 reads in 512 unique sequences.
    ## Sample 60 - 921 reads in 781 unique sequences.
    ## Sample 61 - 1166 reads in 948 unique sequences.
    ##    selfConsist step 2 
    ##    selfConsist step 3 
    ##    selfConsist step 4 
    ##    selfConsist step 5 
    ## Convergence after  5  rounds.
    ## Total reads used:  40925

It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:

``` r
plotErrors(errF, nominalQ=TRUE)
```

![](regression_test_rmarkdown_files/figure-markdown_github/plot-errors-1.png)

### Dereplication

``` r
filtFs %>% 
  basename %>%
  str_replace("_F_filt.fastq.gz","") ->
  sample.names
```

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

### Sample Inference

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 1 reads in 1 unique sequences.
    ## Sample 2 - 672 reads in 566 unique sequences.
    ## Sample 3 - 601 reads in 469 unique sequences.
    ## Sample 4 - 633 reads in 468 unique sequences.
    ## Sample 5 - 567 reads in 394 unique sequences.
    ## Sample 6 - 607 reads in 397 unique sequences.
    ## Sample 7 - 805 reads in 519 unique sequences.
    ## Sample 8 - 547 reads in 360 unique sequences.
    ## Sample 9 - 434 reads in 280 unique sequences.
    ## Sample 10 - 820 reads in 662 unique sequences.
    ## Sample 11 - 607 reads in 365 unique sequences.
    ## Sample 12 - 785 reads in 555 unique sequences.
    ## Sample 13 - 542 reads in 432 unique sequences.
    ## Sample 14 - 449 reads in 380 unique sequences.
    ## Sample 15 - 331 reads in 254 unique sequences.
    ## Sample 16 - 1020 reads in 785 unique sequences.
    ## Sample 17 - 790 reads in 655 unique sequences.
    ## Sample 18 - 1121 reads in 740 unique sequences.
    ## Sample 19 - 965 reads in 744 unique sequences.
    ## Sample 20 - 1203 reads in 905 unique sequences.
    ## Sample 21 - 990 reads in 718 unique sequences.
    ## Sample 22 - 1195 reads in 877 unique sequences.
    ## Sample 23 - 1144 reads in 888 unique sequences.
    ## Sample 24 - 783 reads in 413 unique sequences.
    ## Sample 25 - 723 reads in 428 unique sequences.
    ## Sample 26 - 1018 reads in 552 unique sequences.
    ## Sample 27 - 819 reads in 397 unique sequences.
    ## Sample 28 - 387 reads in 169 unique sequences.
    ## Sample 29 - 503 reads in 229 unique sequences.
    ## Sample 30 - 1 reads in 1 unique sequences.
    ## Sample 31 - 688 reads in 306 unique sequences.
    ## Sample 32 - 656 reads in 277 unique sequences.
    ## Sample 33 - 6 reads in 6 unique sequences.
    ## Sample 34 - 792 reads in 394 unique sequences.
    ## Sample 35 - 1 reads in 1 unique sequences.
    ## Sample 36 - 1 reads in 1 unique sequences.
    ## Sample 37 - 1 reads in 1 unique sequences.
    ## Sample 38 - 415 reads in 270 unique sequences.
    ## Sample 39 - 533 reads in 290 unique sequences.
    ## Sample 40 - 1 reads in 1 unique sequences.
    ## Sample 41 - 53 reads in 52 unique sequences.
    ## Sample 42 - 645 reads in 434 unique sequences.
    ## Sample 43 - 251 reads in 180 unique sequences.
    ## Sample 44 - 901 reads in 751 unique sequences.
    ## Sample 45 - 986 reads in 714 unique sequences.
    ## Sample 46 - 606 reads in 379 unique sequences.
    ## Sample 47 - 166 reads in 136 unique sequences.
    ## Sample 48 - 720 reads in 495 unique sequences.
    ## Sample 49 - 970 reads in 723 unique sequences.
    ## Sample 50 - 1278 reads in 934 unique sequences.
    ## Sample 51 - 951 reads in 739 unique sequences.
    ## Sample 52 - 892 reads in 645 unique sequences.
    ## Sample 53 - 846 reads in 681 unique sequences.
    ## Sample 54 - 732 reads in 548 unique sequences.
    ## Sample 55 - 1093 reads in 674 unique sequences.
    ## Sample 56 - 1190 reads in 839 unique sequences.
    ## Sample 57 - 872 reads in 665 unique sequences.
    ## Sample 58 - 917 reads in 598 unique sequences.
    ## Sample 59 - 612 reads in 404 unique sequences.
    ## Sample 60 - 921 reads in 657 unique sequences.
    ## Sample 61 - 1166 reads in 829 unique sequences.

``` r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 1 reads in 1 unique sequences.
    ## Sample 2 - 672 reads in 592 unique sequences.
    ## Sample 3 - 601 reads in 519 unique sequences.
    ## Sample 4 - 633 reads in 527 unique sequences.
    ## Sample 5 - 567 reads in 444 unique sequences.
    ## Sample 6 - 607 reads in 473 unique sequences.
    ## Sample 7 - 805 reads in 624 unique sequences.
    ## Sample 8 - 547 reads in 441 unique sequences.
    ## Sample 9 - 434 reads in 359 unique sequences.
    ## Sample 10 - 820 reads in 715 unique sequences.
    ## Sample 11 - 607 reads in 463 unique sequences.
    ## Sample 12 - 785 reads in 657 unique sequences.
    ## Sample 13 - 542 reads in 473 unique sequences.
    ## Sample 14 - 449 reads in 404 unique sequences.
    ## Sample 15 - 331 reads in 284 unique sequences.
    ## Sample 16 - 1020 reads in 893 unique sequences.
    ## Sample 17 - 790 reads in 711 unique sequences.
    ## Sample 18 - 1121 reads in 905 unique sequences.
    ## Sample 19 - 965 reads in 839 unique sequences.
    ## Sample 20 - 1203 reads in 1066 unique sequences.
    ## Sample 21 - 990 reads in 818 unique sequences.
    ## Sample 22 - 1195 reads in 1021 unique sequences.
    ## Sample 23 - 1144 reads in 985 unique sequences.
    ## Sample 24 - 783 reads in 526 unique sequences.
    ## Sample 25 - 723 reads in 528 unique sequences.
    ## Sample 26 - 1018 reads in 708 unique sequences.
    ## Sample 27 - 819 reads in 552 unique sequences.
    ## Sample 28 - 387 reads in 244 unique sequences.
    ## Sample 29 - 503 reads in 319 unique sequences.
    ## Sample 30 - 1 reads in 1 unique sequences.
    ## Sample 31 - 688 reads in 436 unique sequences.
    ## Sample 32 - 656 reads in 432 unique sequences.
    ## Sample 33 - 6 reads in 6 unique sequences.
    ## Sample 34 - 792 reads in 536 unique sequences.
    ## Sample 35 - 1 reads in 1 unique sequences.
    ## Sample 36 - 1 reads in 1 unique sequences.
    ## Sample 37 - 1 reads in 1 unique sequences.
    ## Sample 38 - 415 reads in 322 unique sequences.
    ## Sample 39 - 533 reads in 377 unique sequences.
    ## Sample 40 - 1 reads in 1 unique sequences.
    ## Sample 41 - 53 reads in 52 unique sequences.
    ## Sample 42 - 645 reads in 505 unique sequences.
    ## Sample 43 - 251 reads in 214 unique sequences.
    ## Sample 44 - 901 reads in 816 unique sequences.
    ## Sample 45 - 986 reads in 801 unique sequences.
    ## Sample 46 - 606 reads in 446 unique sequences.
    ## Sample 47 - 166 reads in 154 unique sequences.
    ## Sample 48 - 720 reads in 572 unique sequences.
    ## Sample 49 - 970 reads in 811 unique sequences.
    ## Sample 50 - 1278 reads in 1048 unique sequences.
    ## Sample 51 - 951 reads in 809 unique sequences.
    ## Sample 52 - 892 reads in 737 unique sequences.
    ## Sample 53 - 846 reads in 751 unique sequences.
    ## Sample 54 - 732 reads in 609 unique sequences.
    ## Sample 55 - 1093 reads in 762 unique sequences.
    ## Sample 56 - 1190 reads in 946 unique sequences.
    ## Sample 57 - 872 reads in 755 unique sequences.
    ## Sample 58 - 917 reads in 683 unique sequences.
    ## Sample 59 - 612 reads in 512 unique sequences.
    ## Sample 60 - 921 reads in 781 unique sequences.
    ## Sample 61 - 1166 reads in 948 unique sequences.

Inspecting the dada-class object returned by dada:

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1 sample sequences were inferred from 1 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, BAND_SIZE = 16, USE_QUALS = TRUE

### Merge paired reads

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[2]])
```

    ##                                                                                                                                                                                                                                     sequence
    ## 2  GCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCGTGTAGGCGGTTCGGTAAGTCTGCCGTGAAAACCTGGGGCTCAACCCCGGGCGTGCGGTGGATACTGCCGGGCTAGAGGATGGTAGAGGCGAGTGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACACCAGTAGCGAAGGCGGCTCGCTGGGCCATTCCTGACGCTGAGACGCGAAAGCTAGGGG
    ## 5  GCAAGCGTTGTCCGGAATCATTGGGCGTAAAGAGCGCGTAGGCGGCCCGACAAGTCCGCTGTGAAAGTCAGGGGCTTAACCCTTGAATGCCGGTGGATACTGTCGGGCTAGAGTCCGGAAGAGGCGAGTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGCTCGCTGGGACGGTACTGACGCTGAGGCGCGAAAGCGTGGGG
    ## 6  GCTAGCGTTGTTCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGCTTTGTAAGTCGGGGGTGAAAGCCTGTGGCTCAACCACAGAATTGCCTTCGATACTGCATGGCTTGAGACCGGAAGAGGTAAGTGGAACTGCGAGTGTAGAGGTGAAATGCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTTACTGGTCCGGATCTGACGCTGAGGCGCGAAAGCGTGGGG
    ## 11 GCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGAAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGG
    ## 13 GCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGCCTATTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGCCATTGGAAACTGGTAGGCTTGAGTGCAGGAGAGGAGAGCGGAATTCCCGGTGTAGCGGTGAAATGCGTAGATATCGGGAGGAACACCCGTGGCGAAGGCGGCTCTCTGGCCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGG
    ## 16 GCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGCCTGATAAGTAGGGGGTGAAATCCTGCGGCTTAACCGCAGGGCTGCCTTCTAAACTGTCAGGCTCGAGCACAGTAGAGGCAGGTGGAATTCCCGGTGTAGCGGTGGAATGCGTAGAGATCGGGAAGAACATCAGTGGCGAAGGCGGCCTGCTGGGCTGTTGCTGACGCTGAGGCGCGACAGCGTGGGG
    ##    abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 2         58       2       1     32         0      0      2   TRUE
    ## 5         28      12       6     32         0      0      2   TRUE
    ## 6         26       4      10     32         0      0      1   TRUE
    ## 11        20       1       2     32         0      0      2   TRUE
    ## 13        16      16       5     32         0      0      2   TRUE
    ## 16        15       3       3     32         0      0      1   TRUE

### Construct sequence table

``` r
seqtab <- makeSequenceTable(mergers)
```

    ## The sequences being tabled vary in length.

``` r
dim(seqtab)
```

    ## [1]  61 461

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 232 233 234 235 236 
    ##   8 419  32   1   1

### Remove chimeras

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

    ## [1]  61 450

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.9936737

### Track reads through the pipeline

As a final check of our progress, we'll look at the number of reads that made it through each step in the pipeline:

``` r
getN <- function(x) sum(getUniques(x))
track = filt.out %>%
  as.data.frame %>%
  rownames_to_column %>%
  mutate(rowname=str_replace(rowname, ".forward.fastq","")) %>%
  rename(sample=rowname, input=reads.in, filtered=reads.out)

sapply(dadaFs, getN) %>%
  as.tibble %>%
  rownames_to_column() %>%
  rename(sample=rowname, denoised=value) ->
  denoised
```

    ## Warning: `as.tibble()` is deprecated, use `as_tibble()` (but mind the new semantics).
    ## This warning is displayed once per session.

``` r
track %<>% full_join(denoised, by=c("sample"))

sapply(mergers, getN) %>%
  as.tibble %>%
  rownames_to_column() %>%
  rename(sample=rowname, merged=value) ->
  merged

track %<>% full_join(merged, by=c("sample"))


rowSums(seqtab) %>%
  as.tibble %>%
  rownames_to_column() %>%
  rename(sample=rowname, tabled=value) -> 
  tabled
#   denoised

track %<>% full_join(tabled, by=c("sample"))

rowSums(seqtab.nochim) %>%
  as.tibble %>%
  rownames_to_column() %>%
  rename(sample=rowname, nonchim=value) -> 
  nonchim

track %<>% full_join(nonchim, by=c("sample"))

track
```

    ##          sample input filtered denoised merged tabled nonchim
    ## 1   BAQ1370.1.3     1        0       NA     NA     NA      NA
    ## 2   BAQ1552.1.1     1        1       NA     NA     NA      NA
    ## 3   BAQ2420.1.1   721      672       NA     NA     NA      NA
    ## 4   BAQ2420.1.2   640      601       NA     NA     NA      NA
    ## 5   BAQ2420.1.3   656      633       NA     NA     NA      NA
    ## 6     BAQ2420.2   611      567       NA     NA     NA      NA
    ## 7     BAQ2420.3   647      607       NA     NA     NA      NA
    ## 8     BAQ2462.1   869      805       NA     NA     NA      NA
    ## 9     BAQ2462.2   598      547       NA     NA     NA      NA
    ## 10    BAQ2462.3   467      434       NA     NA     NA      NA
    ## 11    BAQ2687.1   869      820       NA     NA     NA      NA
    ## 12    BAQ2687.2   659      607       NA     NA     NA      NA
    ## 13    BAQ2687.3   880      785       NA     NA     NA      NA
    ## 14    BAQ2838.1   588      542       NA     NA     NA      NA
    ## 15    BAQ2838.2   489      449       NA     NA     NA      NA
    ## 16    BAQ2838.3   369      331       NA     NA     NA      NA
    ## 17    BAQ3473.1  1097     1020       NA     NA     NA      NA
    ## 18    BAQ3473.2   824      790       NA     NA     NA      NA
    ## 19    BAQ3473.3  1173     1121       NA     NA     NA      NA
    ## 20  BAQ4166.1.1  1038      965       NA     NA     NA      NA
    ## 21  BAQ4166.1.2  1304     1203       NA     NA     NA      NA
    ## 22  BAQ4166.1.3  1052      990       NA     NA     NA      NA
    ## 23    BAQ4166.2  1264     1195       NA     NA     NA      NA
    ## 24    BAQ4166.3  1236     1144       NA     NA     NA      NA
    ## 25    BAQ4697.1   820      783       NA     NA     NA      NA
    ## 26    BAQ4697.2   761      723       NA     NA     NA      NA
    ## 27    BAQ4697.3  1084     1018       NA     NA     NA      NA
    ## 28  YUN1005.1.1   873      819       NA     NA     NA      NA
    ## 29    YUN1005.3   404      387       NA     NA     NA      NA
    ## 30    YUN1242.1   536      503       NA     NA     NA      NA
    ## 31    YUN1242.2     1        1       NA     NA     NA      NA
    ## 32    YUN1242.3   735      688       NA     NA     NA      NA
    ## 33    YUN1609.1   729      656       NA     NA     NA      NA
    ## 34    YUN2029.1     6        6       NA     NA     NA      NA
    ## 35    YUN2029.2   834      792       NA     NA     NA      NA
    ## 36    YUN2029.3     1        1       NA     NA     NA      NA
    ## 37  YUN3008.1.3     1        1       NA     NA     NA      NA
    ## 38    YUN3008.3     1        1       NA     NA     NA      NA
    ## 39    YUN3153.2   444      415       NA     NA     NA      NA
    ## 40    YUN3153.3   570      533       NA     NA     NA      NA
    ## 41    YUN3184.2     1        1       NA     NA     NA      NA
    ## 42  YUN3259.1.1    56       53       NA     NA     NA      NA
    ## 43  YUN3259.1.2   683      645       NA     NA     NA      NA
    ## 44  YUN3259.1.3   272      251       NA     NA     NA      NA
    ## 45    YUN3259.2  1011      901       NA     NA     NA      NA
    ## 46    YUN3259.3  1081      986       NA     NA     NA      NA
    ## 47    YUN3346.1   652      606       NA     NA     NA      NA
    ## 48    YUN3346.2   171      166       NA     NA     NA      NA
    ## 49    YUN3346.3   793      720       NA     NA     NA      NA
    ## 50    YUN3428.1  1066      970       NA     NA     NA      NA
    ## 51    YUN3428.2  1379     1278       NA     NA     NA      NA
    ## 52    YUN3428.3  1014      951       NA     NA     NA      NA
    ## 53  YUN3533.1.1   937      892       NA     NA     NA      NA
    ## 54  YUN3533.1.2   886      846       NA     NA     NA      NA
    ## 55  YUN3533.1.3   765      732       NA     NA     NA      NA
    ## 56    YUN3533.2  1159     1093       NA     NA     NA      NA
    ## 57    YUN3533.3  1253     1190       NA     NA     NA      NA
    ## 58  YUN3856.1.1   952      872       NA     NA     NA      NA
    ## 59  YUN3856.1.2   971      917       NA     NA     NA      NA
    ## 60  YUN3856.1.3   650      612       NA     NA     NA      NA
    ## 61    YUN3856.2   983      921       NA     NA     NA      NA
    ## 62    YUN3856.3  1233     1166       NA     NA     NA      NA
    ## 63            1    NA       NA        1      1      1       1
    ## 64            2    NA       NA      672    185    185     185
    ## 65            3    NA       NA      601    310    310     310
    ## 66            4    NA       NA      633    322    322     322
    ## 67            5    NA       NA      567    207    207     207
    ## 68            6    NA       NA      607    264    264     264
    ## 69            7    NA       NA      805    458    458     458
    ## 70            8    NA       NA      547    382    382     382
    ## 71            9    NA       NA      434    184    184     184
    ## 72           10    NA       NA      820    286    286     244
    ## 73           11    NA       NA      607    394    394     394
    ## 74           12    NA       NA      785    314    314     314
    ## 75           13    NA       NA      542    105    105     105
    ## 76           14    NA       NA      449     87     87      87
    ## 77           15    NA       NA      331    105    105     105
    ## 78           16    NA       NA     1020    520    520     520
    ## 79           17    NA       NA      790    261    261     261
    ## 80           18    NA       NA     1121    424    424     424
    ## 81           19    NA       NA      965    252    252     252
    ## 82           20    NA       NA     1203    515    515     515
    ## 83           21    NA       NA      990    434    434     434
    ## 84           22    NA       NA     1195    391    391     391
    ## 85           23    NA       NA     1144    556    556     556
    ## 86           24    NA       NA      783    545    545     542
    ## 87           25    NA       NA      723    420    420     416
    ## 88           26    NA       NA     1018    755    755     731
    ## 89           27    NA       NA      819    613    613     613
    ## 90           28    NA       NA      387    302    302     302
    ## 91           29    NA       NA      503    410    410     399
    ## 92           30    NA       NA        1      1      1       1
    ## 93           31    NA       NA      688    551    551     548
    ## 94           32    NA       NA      656    493    493     493
    ## 95           33    NA       NA        6      0      0       0
    ## 96           34    NA       NA      792    578    578     572
    ## 97           35    NA       NA        1      1      1       1
    ## 98           36    NA       NA        1      0      0       0
    ## 99           37    NA       NA        1      1      1       1
    ## 100          38    NA       NA      415    265    265     258
    ## 101          39    NA       NA      533    341    341     341
    ## 102          40    NA       NA        1      0      0       0
    ## 103          41    NA       NA       53      0      0       0
    ## 104          42    NA       NA      645    344    344     344
    ## 105          43    NA       NA      251     89     89      89
    ## 106          44    NA       NA      901    207    207     207
    ## 107          45    NA       NA      986    446    446     446
    ## 108          46    NA       NA      606    427    427     399
    ## 109          47    NA       NA      166     79     79      79
    ## 110          48    NA       NA      720    455    455     455
    ## 111          49    NA       NA      970    350    350     350
    ## 112          50    NA       NA     1278    603    603     603
    ## 113          51    NA       NA      951    263    263     263
    ## 114          52    NA       NA      892    522    522     522
    ## 115          53    NA       NA      846    369    369     369
    ## 116          54    NA       NA      732    335    335     335
    ## 117          55    NA       NA     1093    730    730     730
    ## 118          56    NA       NA     1190    625    625     625
    ## 119          57    NA       NA      872    451    451     451
    ## 120          58    NA       NA      917    547    547     547
    ## 121          59    NA       NA      612    314    314     314
    ## 122          60    NA       NA      921    545    545     543
    ## 123          61    NA       NA     1166    615    615     615

### Assign taxonomy

``` r
taxa <- assignTaxonomy(seqtab.nochim, silva.ref, multithread=TRUE)
```

------------------------------------------------------------------------

Phyloseq
--------

### Make Phyloseq Object

``` r
library(phyloseq)
meta.df = read_tsv(map.file, comment= "#q2") %>%
  rename(Sample = "#SampleID") %>%
  column_to_rownames("Sample") %>%
  as.data.frame
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   `#SampleID` = col_character(),
    ##   BarcodeSequence = col_character(),
    ##   LinkerPrimerSequence = col_character(),
    ##   ExtractGroupNo = col_character(),
    ##   TransectName = col_character(),
    ##   SiteName = col_character(),
    ##   Vegetation = col_character(),
    ##   Description = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
meta.df
```

    ##             BarcodeSequence   LinkerPrimerSequence Elevation ExtractConcen
    ## BAQ1370.1.2    GCCCAAGTTCAC CCGGACTACHVGGGTWTCTAAT      1370         0.019
    ## BAQ1370.3      GCGCCGAATCTT CCGGACTACHVGGGTWTCTAAT      1370         0.124
    ## BAQ1370.1.3    ATAAAGAGGAGG CCGGACTACHVGGGTWTCTAAT      1370         1.200
    ## BAQ1552.1.1    ATCCCAGCATGC CCGGACTACHVGGGTWTCTAAT      1552         0.722
    ## BAQ1552.2      GCTTCCAGACAA CCGGACTACHVGGGTWTCTAAT      1552         0.017
    ## BAQ2420.1.1    ACACAGTCCTGA CCGGACTACHVGGGTWTCTAAT      2420         0.350
    ## BAQ2420.1.2    ATTATACGGCGC CCGGACTACHVGGGTWTCTAAT      2420         0.108
    ## BAQ2420.2      TAAACGCGACTC CCGGACTACHVGGGTWTCTAAT      2420         0.089
    ## BAQ2420.3      CCTCGGGTACTA CCGGACTACHVGGGTWTCTAAT      2420         0.083
    ## BAQ2420.1.3    ATTCAGATGGCA CCGGACTACHVGGGTWTCTAAT      2420         0.132
    ## BAQ2462.1      TTCACCTGTATC CCGGACTACHVGGGTWTCTAAT      2462         1.166
    ## BAQ2462.2      CTCCAGGTCATG CCGGACTACHVGGGTWTCTAAT      2462         0.049
    ## BAQ2462.3      CAGGATTCGTAC CCGGACTACHVGGGTWTCTAAT      2462         0.066
    ## BAQ2687.1      GTGGCCTACTAC CCGGACTACHVGGGTWTCTAAT      2687         1.223
    ## BAQ2687.2      TTCCCTTCTCCG CCGGACTACHVGGGTWTCTAAT      2687         0.251
    ## BAQ2687.3      CATTTGACGACG CCGGACTACHVGGGTWTCTAAT      2687         1.176
    ## BAQ2838.1      CGCATACGACCT CCGGACTACHVGGGTWTCTAAT      2838         0.467
    ## BAQ2838.2      GCCTCGTACTGA CCGGACTACHVGGGTWTCTAAT      2838         0.309
    ## BAQ2838.3      ACCAACAGATTG CCGGACTACHVGGGTWTCTAAT      2838         0.145
    ## BAQ3473.1      AAGTGAAGCGAG CCGGACTACHVGGGTWTCTAAT      3473        16.374
    ## BAQ3473.2      TGCCGCCGTAAT CCGGACTACHVGGGTWTCTAAT      3473         5.170
    ## BAQ3473.3      AACCTCGGATAA CCGGACTACHVGGGTWTCTAAT      3473         0.910
    ## BAQ4166.1.1    GTGCTTGTGTAG CCGGACTACHVGGGTWTCTAAT      4166        33.492
    ## BAQ4166.1.2    CAACTAGACTCG CCGGACTACHVGGGTWTCTAAT      4166        16.086
    ## BAQ4166.2      GGAACGACGTGA CCGGACTACHVGGGTWTCTAAT      4166         2.958
    ## BAQ4166.3      TGTCAGCTGTCG CCGGACTACHVGGGTWTCTAAT      4166        28.920
    ## BAQ4166.1.3    AGTGCCCTTGGT CCGGACTACHVGGGTWTCTAAT      4166        14.532
    ## BAQ4697.1      CTGGTGCTGAAT CCGGACTACHVGGGTWTCTAAT      4697        10.044
    ## BAQ4697.2      GACAGAGGTGCA CCGGACTACHVGGGTWTCTAAT      4697         5.864
    ## BAQ4697.3      TCAGACCAACTG CCGGACTACHVGGGTWTCTAAT      4697         8.264
    ## BAQ895.2       TCCTCACTATCA CCGGACTACHVGGGTWTCTAAT       895         0.017
    ## BAQ895.3       GCCTGCAGTACT CCGGACTACHVGGGTWTCTAAT       895         0.011
    ## YUN1005.1.1    ACGTAACCACGT CCGGACTACHVGGGTWTCTAAT      1005         0.103
    ## YUN1005.2      TCTAACGAGTGC CCGGACTACHVGGGTWTCTAAT      1005         0.030
    ## YUN1005.3      CATCTGGGCAAT CCGGACTACHVGGGTWTCTAAT      1005         0.028
    ## YUN1005.1.3    GTCGGAAATTGT CCGGACTACHVGGGTWTCTAAT      1005         0.034
    ## YUN1242.1      TGTCCGTGGATC CCGGACTACHVGGGTWTCTAAT      1242         0.228
    ## YUN1242.2      ACTCGGCCAACT CCGGACTACHVGGGTWTCTAAT      1242         0.023
    ## YUN1242.3      GTTGGTTGGCAT CCGGACTACHVGGGTWTCTAAT      1242         0.050
    ## YUN1609.1      TTCCACACGTGG CCGGACTACHVGGGTWTCTAAT      1609         0.031
    ## YUN1609.3      AACCCAGATGAT CCGGACTACHVGGGTWTCTAAT      1609         0.025
    ## YUN2029.1      GTAGTGTCAACA CCGGACTACHVGGGTWTCTAAT      2029         0.022
    ## YUN2029.2      TGGAGAGGAGAT CCGGACTACHVGGGTWTCTAAT      2029         0.232
    ## YUN2029.3      CGTATAAATGCG CCGGACTACHVGGGTWTCTAAT      2029         0.078
    ## YUN3008.1.1    ATCGATCCACAG CCGGACTACHVGGGTWTCTAAT      3008         0.021
    ## YUN3008.1.2    ACACCGCACAAT CCGGACTACHVGGGTWTCTAAT      3008         0.029
    ## YUN3008.2      GTAGCACTCATG CCGGACTACHVGGGTWTCTAAT      3008         0.036
    ## YUN3008.3      CACCTGTAGTAG CCGGACTACHVGGGTWTCTAAT      3008         0.021
    ## YUN3008.1.3    GTCTCCTCCCTT CCGGACTACHVGGGTWTCTAAT      3008         0.017
    ## YUN3153.1      AATACAGACCTG CCGGACTACHVGGGTWTCTAAT      3153         0.011
    ## YUN3153.2      GACTCAACCAGT CCGGACTACHVGGGTWTCTAAT      3153         0.493
    ## YUN3153.3      GGAAGAAGTAGC CCGGACTACHVGGGTWTCTAAT      3153         0.220
    ## YUN3184.1      CACGAGCTACTC CCGGACTACHVGGGTWTCTAAT      3184         0.022
    ## YUN3184.2      TCTCGATAAGCG CCGGACTACHVGGGTWTCTAAT      3184         0.010
    ## YUN3259.1.1    TAGACACCGTGT CCGGACTACHVGGGTWTCTAAT      3259         3.710
    ## YUN3259.1.2    AGACAAGCTTCC CCGGACTACHVGGGTWTCTAAT      3259         1.543
    ## YUN3259.2      TCACTTGGTGCG CCGGACTACHVGGGTWTCTAAT      3259         1.935
    ## YUN3259.3      TTATGTACGGCG CCGGACTACHVGGGTWTCTAAT      3259         3.530
    ## YUN3259.1.3    TCCGCAACCTGA CCGGACTACHVGGGTWTCTAAT      3259         0.329
    ## YUN3346.1      GAGATACAGTTC CCGGACTACHVGGGTWTCTAAT      3346         0.314
    ## YUN3346.2      GCATGCATCCCA CCGGACTACHVGGGTWTCTAAT      3346         0.318
    ## YUN3346.3      GATCTAATCGAG CCGGACTACHVGGGTWTCTAAT      3346         1.218
    ## YUN3428.1      TACGCCCATCAG CCGGACTACHVGGGTWTCTAAT      3428         1.764
    ## YUN3428.2      AAGATCGTACTG CCGGACTACHVGGGTWTCTAAT      3428         1.988
    ## YUN3428.3      ACTCATCTTCCA CCGGACTACHVGGGTWTCTAAT      3428         2.782
    ## YUN3533.1.1    TTGGACGTCCAC CCGGACTACHVGGGTWTCTAAT      3533         3.352
    ## YUN3533.1.2    TCCAGGGCTATA CCGGACTACHVGGGTWTCTAAT      3533         3.210
    ## YUN3533.2      GAAACTCCTAGA CCGGACTACHVGGGTWTCTAAT      3533         9.236
    ## YUN3533.3      ATCGGGCTTAAC CCGGACTACHVGGGTWTCTAAT      3533         3.642
    ## YUN3533.1.3    GCGTAGAGAGAC CCGGACTACHVGGGTWTCTAAT      3533         1.987
    ## YUN3856.1.1    AATCTTGCGCCG CCGGACTACHVGGGTWTCTAAT      3856         6.444
    ## YUN3856.1.2    GGAAATCCCATC CCGGACTACHVGGGTWTCTAAT      3856         6.230
    ## YUN3856.2      TTGGAACGGCTT CCGGACTACHVGGGTWTCTAAT      3856         4.178
    ## YUN3856.3      TCCTAGGTCCGA CCGGACTACHVGGGTWTCTAAT      3856         6.695
    ## YUN3856.1.3    GACCGTCAATAC CCGGACTACHVGGGTWTCTAAT      3856         3.134
    ##             AmpliconConcentration ExtractGroupNo TransectName SiteName
    ## BAQ1370.1.2                  0.95              B    Baquedano  BAQ1370
    ## BAQ1370.3                   17.46              E    Baquedano  BAQ1370
    ## BAQ1370.1.3                  0.96              J    Baquedano  BAQ1370
    ## BAQ1552.1.1                 18.83              J    Baquedano  BAQ1552
    ## BAQ1552.2                    2.00              B    Baquedano  BAQ1552
    ## BAQ2420.1.1                  7.40              H    Baquedano  BAQ2420
    ## BAQ2420.1.2                 11.42              H    Baquedano  BAQ2420
    ## BAQ2420.2                   10.06              A    Baquedano  BAQ2420
    ## BAQ2420.3                   15.92              C    Baquedano  BAQ2420
    ## BAQ2420.1.3                 10.00              H    Baquedano  BAQ2420
    ## BAQ2462.1                   15.29              H    Baquedano  BAQ2462
    ## BAQ2462.2                   10.39              H    Baquedano  BAQ2462
    ## BAQ2462.3                    8.46              A    Baquedano  BAQ2462
    ## BAQ2687.1                    1.24              A    Baquedano  BAQ2687
    ## BAQ2687.2                   17.53              C    Baquedano  BAQ2687
    ## BAQ2687.3                   12.10              E    Baquedano  BAQ2687
    ## BAQ2838.1                   10.91              G    Baquedano  BAQ2838
    ## BAQ2838.2                   14.38              H    Baquedano  BAQ2838
    ## BAQ2838.3                   14.06              A    Baquedano  BAQ2838
    ## BAQ3473.1                    1.29              C    Baquedano  BAQ3473
    ## BAQ3473.2                   11.05              B    Baquedano  BAQ3473
    ## BAQ3473.3                    4.41              C    Baquedano  BAQ3473
    ## BAQ4166.1.1                  3.16              C    Baquedano  BAQ4166
    ## BAQ4166.1.2                 18.25              B    Baquedano  BAQ4166
    ## BAQ4166.2                    0.12              C    Baquedano  BAQ4166
    ## BAQ4166.3                    0.26              C    Baquedano  BAQ4166
    ## BAQ4166.1.3                  0.55              C    Baquedano  BAQ4166
    ## BAQ4697.1                    1.88              C    Baquedano  BAQ4697
    ## BAQ4697.2                    1.22              B    Baquedano  BAQ4697
    ## BAQ4697.3                    1.39              C    Baquedano  BAQ4697
    ## BAQ895.2                     1.01              B    Baquedano   BAQ895
    ## BAQ895.3                     2.36              E    Baquedano   BAQ895
    ## YUN1005.1.1                 13.17              I       Yungay  YUN1005
    ## YUN1005.2                    0.94              I       Yungay  YUN1005
    ## YUN1005.3                    4.32              I       Yungay  YUN1005
    ## YUN1005.1.3                  1.12              I       Yungay  YUN1005
    ## YUN1242.1                   19.20              K       Yungay  YUN1242
    ## YUN1242.2                    0.54              B       Yungay  YUN1242
    ## YUN1242.3                   15.61              K       Yungay  YUN1242
    ## YUN1609.1                    6.42              B       Yungay  YUN1609
    ## YUN1609.3                   12.84              E       Yungay  YUN1609
    ## YUN2029.1                    9.59              K       Yungay  YUN2029
    ## YUN2029.2                   16.71              B       Yungay  YUN2029
    ## YUN2029.3                   16.14              K       Yungay  YUN2029
    ## YUN3008.1.1                  1.29              G       Yungay  YUN3008
    ## YUN3008.1.2                  1.23              B       Yungay  YUN3008
    ## YUN3008.2                    1.10              K       Yungay  YUN3008
    ## YUN3008.3                    1.00              G       Yungay  YUN3008
    ## YUN3008.1.3                  0.99              G       Yungay  YUN3008
    ## YUN3153.1                    1.05              E       Yungay  YUN3153
    ## YUN3153.2                   16.89              G       Yungay  YUN3153
    ## YUN3153.3                   18.41              B       Yungay  YUN3153
    ## YUN3184.1                    1.18              B       Yungay  YUN3184
    ## YUN3184.2                    1.29              E       Yungay  YUN3184
    ## YUN3259.1.1                  7.07              D       Yungay  YUN3259
    ## YUN3259.1.2                 16.92              B       Yungay  YUN3259
    ## YUN3259.2                   12.11              D       Yungay  YUN3259
    ## YUN3259.3                   16.74              D       Yungay  YUN3259
    ## YUN3259.1.3                 19.17              D       Yungay  YUN3259
    ## YUN3346.1                   16.88              D       Yungay  YUN3346
    ## YUN3346.2                   16.26              D       Yungay  YUN3346
    ## YUN3346.3                   17.14              D       Yungay  YUN3346
    ## YUN3428.1                   15.32              D       Yungay  YUN3428
    ## YUN3428.2                    6.61              D       Yungay  YUN3428
    ## YUN3428.3                   15.46              D       Yungay  YUN3428
    ## YUN3533.1.1                 14.29              D       Yungay  YUN3533
    ## YUN3533.1.2                 14.32              D       Yungay  YUN3533
    ## YUN3533.2                   14.33              D       Yungay  YUN3533
    ## YUN3533.3                   17.43              D       Yungay  YUN3533
    ## YUN3533.1.3                  7.55              D       Yungay  YUN3533
    ## YUN3856.1.1                 12.93              G       Yungay  YUN3856
    ## YUN3856.1.2                 16.35              G       Yungay  YUN3856
    ## YUN3856.2                   17.12              H       Yungay  YUN3856
    ## YUN3856.3                   15.86              H       Yungay  YUN3856
    ## YUN3856.1.3                 18.17              G       Yungay  YUN3856
    ##             Depth   pH   TOC    EC AverageSoilRelativeHumidity
    ## BAQ1370.1.2     2 7.98   525 6.080                       16.17
    ## BAQ1370.3       2   NA   771 6.080                       16.17
    ## BAQ1370.1.3     3 8.13    NA    NA                       16.17
    ## BAQ1552.1.1     1 7.87    NA    NA                       15.75
    ## BAQ1552.2       2   NA   223 1.839                       15.75
    ## BAQ2420.1.1     1 9.33    NA    NA                       82.54
    ## BAQ2420.1.2     2 9.36   166 0.075                       82.54
    ## BAQ2420.2       2 9.36   337 0.075                       82.54
    ## BAQ2420.3       2 9.36   574 0.075                       82.54
    ## BAQ2420.1.3     3 8.90    NA    NA                       82.54
    ## BAQ2462.1       2 8.35   173 0.506                       69.08
    ## BAQ2462.2       2 8.35   331 0.506                       69.08
    ## BAQ2462.3       2 8.35   634 0.506                       69.08
    ## BAQ2687.1       2 8.35   554 0.131                       73.21
    ## BAQ2687.2       2 8.35   350 0.131                       73.21
    ## BAQ2687.3       2 8.35   400 0.131                       73.21
    ## BAQ2838.1       2 8.36   353 0.212                       44.74
    ## BAQ2838.2       2 8.36   260 0.212                       44.74
    ## BAQ2838.3       2 8.36   369 0.212                       44.74
    ## BAQ3473.1       2 7.93  3675 0.420                       82.05
    ## BAQ3473.2       2 7.93   838 0.420                       82.05
    ## BAQ3473.3       2 7.93 16449 0.420                       82.05
    ## BAQ4166.1.1     1 5.64    NA    NA                      100.00
    ## BAQ4166.1.2     2 7.22  2085 0.084                      100.00
    ## BAQ4166.2       2 7.22  3692 0.084                      100.00
    ## BAQ4166.3       2 7.22  2271 0.084                      100.00
    ## BAQ4166.1.3     3 7.41    NA    NA                      100.00
    ## BAQ4697.1       2 7.44   424 0.055                          NA
    ## BAQ4697.2       2 7.44   531 0.055                          NA
    ## BAQ4697.3       2 7.44   631 0.055                          NA
    ## BAQ895.2        2   NA    61 2.945                       32.22
    ## BAQ895.3        2   NA   474 2.945                       32.22
    ## YUN1005.1.1     1 7.54    NA    NA                       20.70
    ## YUN1005.2       2   NA   160 2.270                       20.70
    ## YUN1005.3       2 7.60   223 2.270                       20.70
    ## YUN1005.1.3     3 7.60    NA    NA                       20.70
    ## YUN1242.1       2 9.00   226 1.845                       20.90
    ## YUN1242.2       2 9.00   194 1.845                       20.90
    ## YUN1242.3       2 9.00   452 1.845                       20.90
    ## YUN1609.1       2 8.00   361 0.427                       17.18
    ## YUN1609.3       2 8.00   424 0.427                       17.18
    ## YUN2029.1       2 7.89   184 0.067                       28.79
    ## YUN2029.2       2 7.89   241 0.067                       28.79
    ## YUN2029.3       2 7.89   817 0.067                       28.79
    ## YUN3008.1.1     1 7.70    NA    NA                       70.89
    ## YUN3008.1.2     2 7.71   250 2.205                       70.89
    ## YUN3008.2       2   NA   367 2.205                       70.89
    ## YUN3008.3       2   NA   417 2.205                       70.89
    ## YUN3008.1.3     3 7.50    NA    NA                       70.89
    ## YUN3153.1       2 7.60   580 2.260                       59.69
    ## YUN3153.2       2 7.60   185 2.260                       59.69
    ## YUN3153.3       2 7.60   391 2.260                       59.69
    ## YUN3184.1       2 7.70   254 2.235                       26.97
    ## YUN3184.2       2   NA   411 2.235                       26.97
    ## YUN3259.1.1     1 8.48    NA    NA                       93.57
    ## YUN3259.1.2     2 8.10   419 0.133                       93.57
    ## YUN3259.2       2 8.10   342 0.133                       93.57
    ## YUN3259.3       2 8.10   539 0.133                       93.57
    ## YUN3259.1.3     3 8.40    NA    NA                       93.57
    ## YUN3346.1       2 7.10   429 0.044                       87.32
    ## YUN3346.2       2 7.10   666 0.044                       87.32
    ## YUN3346.3       2 7.10   387 0.044                       87.32
    ## YUN3428.1       2 7.20   623 0.023                       99.99
    ## YUN3428.2       2 7.20   621 0.023                       99.99
    ## YUN3428.3       2 7.20   658 0.023                       99.99
    ## YUN3533.1.1     1 7.10    NA    NA                      100.00
    ## YUN3533.1.2     2 8.00   521 0.024                      100.00
    ## YUN3533.2       2 8.00   773 0.024                      100.00
    ## YUN3533.3       2 8.00   664 0.024                      100.00
    ## YUN3533.1.3     3 7.80    NA    NA                      100.00
    ## YUN3856.1.1     1 6.98    NA    NA                       99.44
    ## YUN3856.1.2     2 7.43   713 0.029                       99.44
    ## YUN3856.2       2 7.43   404 0.029                       99.44
    ## YUN3856.3       2 7.43   988 0.029                       99.44
    ## YUN3856.1.3     3 7.54    NA    NA                       99.44
    ##             RelativeHumiditySoilHigh RelativeHumiditySoilLow
    ## BAQ1370.1.2                   23.970                  11.420
    ## BAQ1370.3                     23.970                  11.420
    ## BAQ1370.1.3                   23.970                  11.420
    ## BAQ1552.1.1                   35.360                  11.100
    ## BAQ1552.2                     35.360                  11.100
    ## BAQ2420.1.1                  100.000                  41.135
    ## BAQ2420.1.2                  100.000                  41.135
    ## BAQ2420.2                    100.000                  41.135
    ## BAQ2420.3                    100.000                  41.135
    ## BAQ2420.1.3                  100.000                  41.135
    ## BAQ2462.1                     98.396                  32.558
    ## BAQ2462.2                     98.396                  32.558
    ## BAQ2462.3                     98.396                  32.558
    ## BAQ2687.1                    100.000                  20.870
    ## BAQ2687.2                    100.000                  20.870
    ## BAQ2687.3                    100.000                  20.870
    ## BAQ2838.1                    100.000                  16.647
    ## BAQ2838.2                    100.000                  16.647
    ## BAQ2838.3                    100.000                  16.647
    ## BAQ3473.1                    100.000                  42.090
    ## BAQ3473.2                    100.000                  42.090
    ## BAQ3473.3                    100.000                  42.090
    ## BAQ4166.1.1                  100.000                 100.000
    ## BAQ4166.1.2                  100.000                 100.000
    ## BAQ4166.2                    100.000                 100.000
    ## BAQ4166.3                    100.000                 100.000
    ## BAQ4166.1.3                  100.000                 100.000
    ## BAQ4697.1                         NA                      NA
    ## BAQ4697.2                         NA                      NA
    ## BAQ4697.3                         NA                      NA
    ## BAQ895.2                      43.140                  25.600
    ## BAQ895.3                      43.140                  25.600
    ## YUN1005.1.1                   31.011                  16.343
    ## YUN1005.2                     31.011                  16.343
    ## YUN1005.3                     31.011                  16.343
    ## YUN1005.1.3                   31.011                  16.343
    ## YUN1242.1                     36.460                  16.960
    ## YUN1242.2                     36.460                  16.960
    ## YUN1242.3                     36.460                  16.960
    ## YUN1609.1                     50.400                  12.530
    ## YUN1609.3                     50.400                  12.530
    ## YUN2029.1                     52.550                  21.760
    ## YUN2029.2                     52.550                  21.760
    ## YUN2029.3                     52.550                  21.760
    ## YUN3008.1.1                  100.000                  20.370
    ## YUN3008.1.2                  100.000                  20.370
    ## YUN3008.2                    100.000                  20.370
    ## YUN3008.3                    100.000                  20.370
    ## YUN3008.1.3                  100.000                  20.370
    ## YUN3153.1                     93.150                  32.120
    ## YUN3153.2                     93.150                  32.120
    ## YUN3153.3                     93.150                  32.120
    ## YUN3184.1                     46.570                  15.090
    ## YUN3184.2                     46.570                  15.090
    ## YUN3259.1.1                  100.000                  69.040
    ## YUN3259.1.2                  100.000                  69.040
    ## YUN3259.2                    100.000                  69.040
    ## YUN3259.3                    100.000                  69.040
    ## YUN3259.1.3                  100.000                  69.040
    ## YUN3346.1                    100.000                  57.800
    ## YUN3346.2                    100.000                  57.800
    ## YUN3346.3                    100.000                  57.800
    ## YUN3428.1                    100.000                  99.910
    ## YUN3428.2                    100.000                  99.910
    ## YUN3428.3                    100.000                  99.910
    ## YUN3533.1.1                  100.000                  99.930
    ## YUN3533.1.2                  100.000                  99.930
    ## YUN3533.2                    100.000                  99.930
    ## YUN3533.3                    100.000                  99.930
    ## YUN3533.1.3                  100.000                  99.930
    ## YUN3856.1.1                  100.000                  73.750
    ## YUN3856.1.2                  100.000                  73.750
    ## YUN3856.2                    100.000                  73.750
    ## YUN3856.3                    100.000                  73.750
    ## YUN3856.1.3                  100.000                  73.750
    ##             PercentRelativeHumiditySoil_100 AverageSoilTemperature
    ## BAQ1370.1.2                            0.00                  22.61
    ## BAQ1370.3                              0.00                  22.61
    ## BAQ1370.1.3                            0.00                  22.61
    ## BAQ1552.1.1                            0.00                  22.63
    ## BAQ1552.2                              0.00                  22.63
    ## BAQ2420.1.1                           46.77                  22.45
    ## BAQ2420.1.2                           46.77                  22.45
    ## BAQ2420.2                             46.77                  22.45
    ## BAQ2420.3                             46.77                  22.45
    ## BAQ2420.1.3                           46.77                  22.45
    ## BAQ2462.1                              0.21                  21.83
    ## BAQ2462.2                              0.21                  21.83
    ## BAQ2462.3                              0.21                  21.83
    ## BAQ2687.1                             43.48                  18.20
    ## BAQ2687.2                             43.48                  18.20
    ## BAQ2687.3                             43.48                  18.20
    ## BAQ2838.1                             11.34                  19.08
    ## BAQ2838.2                             11.34                  19.08
    ## BAQ2838.3                             11.34                  19.08
    ## BAQ3473.1                             35.37                  13.46
    ## BAQ3473.2                             35.37                  13.46
    ## BAQ3473.3                             35.37                  13.46
    ## BAQ4166.1.1                          100.00                   9.41
    ## BAQ4166.1.2                          100.00                   9.41
    ## BAQ4166.2                            100.00                   9.41
    ## BAQ4166.3                            100.00                   9.41
    ## BAQ4166.1.3                          100.00                   9.41
    ## BAQ4697.1                                NA                     NA
    ## BAQ4697.2                                NA                     NA
    ## BAQ4697.3                                NA                     NA
    ## BAQ895.2                               0.00                  19.17
    ## BAQ895.3                               0.00                  19.17
    ## YUN1005.1.1                            0.00                  22.59
    ## YUN1005.2                              0.00                  22.59
    ## YUN1005.3                              0.00                  22.59
    ## YUN1005.1.3                            0.00                  22.59
    ## YUN1242.1                              0.00                  23.61
    ## YUN1242.2                              0.00                  23.61
    ## YUN1242.3                              0.00                  23.61
    ## YUN1609.1                              0.00                  21.47
    ## YUN1609.3                              0.00                  21.47
    ## YUN2029.1                              0.00                  17.85
    ## YUN2029.2                              0.00                  17.85
    ## YUN2029.3                              0.00                  17.85
    ## YUN3008.1.1                           29.10                  11.57
    ## YUN3008.1.2                           29.10                  11.57
    ## YUN3008.2                             29.10                  11.57
    ## YUN3008.3                             29.10                  11.57
    ## YUN3008.1.3                           29.10                  11.57
    ## YUN3153.1                              0.00                  12.50
    ## YUN3153.2                              0.00                  12.50
    ## YUN3153.3                              0.00                  12.50
    ## YUN3184.1                              0.00                  14.12
    ## YUN3184.2                              0.00                  14.12
    ## YUN3259.1.1                           64.69                  15.28
    ## YUN3259.1.2                           64.69                  15.28
    ## YUN3259.2                             64.69                  15.28
    ## YUN3259.3                             64.69                  15.28
    ## YUN3259.1.3                           64.69                  15.28
    ## YUN3346.1                             33.13                  15.89
    ## YUN3346.2                             33.13                  15.89
    ## YUN3346.3                             33.13                  15.89
    ## YUN3428.1                             99.42                  13.61
    ## YUN3428.2                             99.42                  13.61
    ## YUN3428.3                             99.42                  13.61
    ## YUN3533.1.1                           99.94                  11.41
    ## YUN3533.1.2                           99.94                  11.41
    ## YUN3533.2                             99.94                  11.41
    ## YUN3533.3                             99.94                  11.41
    ## YUN3533.1.3                           99.94                  11.41
    ## YUN3856.1.1                           93.70                   9.51
    ## YUN3856.1.2                           93.70                   9.51
    ## YUN3856.2                             93.70                   9.51
    ## YUN3856.3                             93.70                   9.51
    ## YUN3856.1.3                           93.70                   9.51
    ##             TemperatureSoilHigh TemperatureSoilLow Vegetation PercentCover
    ## BAQ1370.1.2              35.210             12.460         no         0.00
    ## BAQ1370.3                35.210             12.460         no         0.00
    ## BAQ1370.1.3              35.210             12.460         no         0.00
    ## BAQ1552.1.1              30.650             10.960         no         0.00
    ## BAQ1552.2                30.650             10.960         no         0.00
    ## BAQ2420.1.1              28.297             13.294         no         0.00
    ## BAQ2420.1.2              28.297             13.294         no         0.00
    ## BAQ2420.2                28.297             13.294         no         0.00
    ## BAQ2420.3                28.297             13.294         no         0.00
    ## BAQ2420.1.3              28.297             13.294         no         0.00
    ## BAQ2462.1                29.024             12.971         no         0.00
    ## BAQ2462.2                29.024             12.971         no         0.00
    ## BAQ2462.3                29.024             12.971         no         0.00
    ## BAQ2687.1                25.000             10.670        yes         0.10
    ## BAQ2687.2                25.000             10.670        yes         0.10
    ## BAQ2687.3                25.000             10.670        yes         0.10
    ## BAQ2838.1                24.995             10.709         no         0.00
    ## BAQ2838.2                24.995             10.709         no         0.00
    ## BAQ2838.3                24.995             10.709         no         0.00
    ## BAQ3473.1                20.690              5.120        yes         1.20
    ## BAQ3473.2                20.690              5.120        yes         1.20
    ## BAQ3473.3                20.690              5.120        yes         1.20
    ## BAQ4166.1.1              17.540              1.420        yes         7.10
    ## BAQ4166.1.2              17.540              1.420        yes         7.10
    ## BAQ4166.2                17.540              1.420        yes         7.10
    ## BAQ4166.3                17.540              1.420        yes         7.10
    ## BAQ4166.1.3              17.540              1.420        yes         7.10
    ## BAQ4697.1                    NA                 NA        yes         0.10
    ## BAQ4697.2                    NA                 NA        yes         0.10
    ## BAQ4697.3                    NA                 NA        yes         0.10
    ## BAQ895.2                 31.820             10.540         no         0.00
    ## BAQ895.3                 31.820             10.540         no         0.00
    ## YUN1005.1.1              26.260             16.510         no         0.00
    ## YUN1005.2                26.260             16.510         no         0.00
    ## YUN1005.3                26.260             16.510         no         0.00
    ## YUN1005.1.3              26.260             16.510         no         0.00
    ## YUN1242.1                27.620             18.330         no         0.00
    ## YUN1242.2                27.620             18.330         no         0.00
    ## YUN1242.3                27.620             18.330         no         0.00
    ## YUN1609.1                25.390             16.100         no         0.00
    ## YUN1609.3                25.390             16.100         no         0.00
    ## YUN2029.1                24.160             10.210         no         0.00
    ## YUN2029.2                24.160             10.210         no         0.00
    ## YUN2029.3                24.160             10.210         no         0.00
    ## YUN3008.1.1              30.980             -2.570         no         0.00
    ## YUN3008.1.2              30.980             -2.570         no         0.00
    ## YUN3008.2                30.980             -2.570         no         0.00
    ## YUN3008.3                30.980             -2.570         no         0.00
    ## YUN3008.1.3              30.980             -2.570         no         0.00
    ## YUN3153.1                19.620              3.480         no         0.00
    ## YUN3153.2                19.620              3.480         no         0.00
    ## YUN3153.3                19.620              3.480         no         0.00
    ## YUN3184.1                25.210              2.450         no         0.00
    ## YUN3184.2                25.210              2.450         no         0.00
    ## YUN3259.1.1              23.650              6.310        yes         2.40
    ## YUN3259.1.2              23.650              6.310        yes         2.40
    ## YUN3259.2                23.650              6.310        yes         2.40
    ## YUN3259.3                23.650              6.310        yes         2.40
    ## YUN3259.1.3              23.650              6.310        yes         2.40
    ## YUN3346.1                23.960              6.890        yes         0.01
    ## YUN3346.2                23.960              6.890        yes         0.01
    ## YUN3346.3                23.960              6.890        yes         0.01
    ## YUN3428.1                21.090              4.130        yes         8.80
    ## YUN3428.2                21.090              4.130        yes         8.80
    ## YUN3428.3                21.090              4.130        yes         8.80
    ## YUN3533.1.1              19.210              3.470        yes         8.60
    ## YUN3533.1.2              19.210              3.470        yes         8.60
    ## YUN3533.2                19.210              3.470        yes         8.60
    ## YUN3533.3                19.210              3.470        yes         8.60
    ## YUN3533.1.3              19.210              3.470        yes         8.60
    ## YUN3856.1.1              17.450              0.450        yes         3.10
    ## YUN3856.1.2              17.450              0.450        yes         3.10
    ## YUN3856.2                17.450              0.450        yes         3.10
    ## YUN3856.3                17.450              0.450        yes         3.10
    ## YUN3856.1.3              17.450              0.450        yes         3.10
    ##             Description
    ## BAQ1370.1.2 BAQ1370.1.2
    ## BAQ1370.3     BAQ1370.3
    ## BAQ1370.1.3 BAQ1370.1.3
    ## BAQ1552.1.1 BAQ1552.1.1
    ## BAQ1552.2     BAQ1552.2
    ## BAQ2420.1.1 BAQ2420.1.1
    ## BAQ2420.1.2 BAQ2420.1.2
    ## BAQ2420.2     BAQ2420.2
    ## BAQ2420.3     BAQ2420.3
    ## BAQ2420.1.3 BAQ2420.1.3
    ## BAQ2462.1     BAQ2462.1
    ## BAQ2462.2     BAQ2462.2
    ## BAQ2462.3     BAQ2462.3
    ## BAQ2687.1     BAQ2687.1
    ## BAQ2687.2     BAQ2687.2
    ## BAQ2687.3     BAQ2687.3
    ## BAQ2838.1     BAQ2838.1
    ## BAQ2838.2     BAQ2838.2
    ## BAQ2838.3     BAQ2838.3
    ## BAQ3473.1     BAQ3473.1
    ## BAQ3473.2     BAQ3473.2
    ## BAQ3473.3     BAQ3473.3
    ## BAQ4166.1.1 BAQ4166.1.1
    ## BAQ4166.1.2 BAQ4166.1.2
    ## BAQ4166.2     BAQ4166.2
    ## BAQ4166.3     BAQ4166.3
    ## BAQ4166.1.3 BAQ4166.1.3
    ## BAQ4697.1     BAQ4697.1
    ## BAQ4697.2     BAQ4697.2
    ## BAQ4697.3     BAQ4697.3
    ## BAQ895.2       BAQ895.2
    ## BAQ895.3       BAQ895.3
    ## YUN1005.1.1 YUN1005.1.1
    ## YUN1005.2     YUN1005.2
    ## YUN1005.3     YUN1005.3
    ## YUN1005.1.3 YUN1005.1.3
    ## YUN1242.1     YUN1242.1
    ## YUN1242.2     YUN1242.2
    ## YUN1242.3     YUN1242.3
    ## YUN1609.1     YUN1609.1
    ## YUN1609.3     YUN1609.3
    ## YUN2029.1     YUN2029.1
    ## YUN2029.2     YUN2029.2
    ## YUN2029.3     YUN2029.3
    ## YUN3008.1.1 YUN3008.1.1
    ## YUN3008.1.2 YUN3008.1.2
    ## YUN3008.2     YUN3008.2
    ## YUN3008.3     YUN3008.3
    ## YUN3008.1.3 YUN3008.1.3
    ## YUN3153.1     YUN3153.1
    ## YUN3153.2     YUN3153.2
    ## YUN3153.3     YUN3153.3
    ## YUN3184.1     YUN3184.1
    ## YUN3184.2     YUN3184.2
    ## YUN3259.1.1 YUN3259.1.1
    ## YUN3259.1.2 YUN3259.1.2
    ## YUN3259.2     YUN3259.2
    ## YUN3259.3     YUN3259.3
    ## YUN3259.1.3 YUN3259.1.3
    ## YUN3346.1     YUN3346.1
    ## YUN3346.2     YUN3346.2
    ## YUN3346.3     YUN3346.3
    ## YUN3428.1     YUN3428.1
    ## YUN3428.2     YUN3428.2
    ## YUN3428.3     YUN3428.3
    ## YUN3533.1.1 YUN3533.1.1
    ## YUN3533.1.2 YUN3533.1.2
    ## YUN3533.2     YUN3533.2
    ## YUN3533.3     YUN3533.3
    ## YUN3533.1.3 YUN3533.1.3
    ## YUN3856.1.1 YUN3856.1.1
    ## YUN3856.1.2 YUN3856.1.2
    ## YUN3856.2     YUN3856.2
    ## YUN3856.3     YUN3856.3
    ## YUN3856.1.3 YUN3856.1.3

``` r
otus = otu_table(seqtab.nochim, taxa_are_rows=FALSE)

sd = sample_data(meta.df)
ps <- phyloseq(otus,
               sd,
               tax_table(taxa))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 450 taxa and 61 samples ]
    ## sample_data() Sample Data:       [ 61 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 450 taxa by 6 taxonomic ranks ]

### Visualize alpha-diversity

``` r
plot_richness(ps, x="Elevation", measures=c("Shannon", "Simpson"), color="TransectName") + theme_bw()
```

![](regression_test_rmarkdown_files/figure-markdown_github/richness-1.png)

CRAN Packages
=============

multcomp
--------

``` r
### R code from vignette source 'multcomp-examples.Rnw'
library("multcomp")
```

    ## Loading required package: mvtnorm

    ## Loading required package: survival

    ## Loading required package: TH.data

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 
    ## Attaching package: 'TH.data'

    ## The following object is masked from 'package:MASS':
    ## 
    ##     geyser

``` r
library("mvtnorm")

dig <- 4
options(width = 65, digits = dig)
set.seed(290875)

lm.cars <- lm(dist ~ speed, data = cars)
summary(lm.cars)
```

    ## 
    ## Call:
    ## lm(formula = dist ~ speed, data = cars)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -29.07  -9.53  -2.27   9.21  43.20 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -17.579      6.758   -2.60    0.012 *  
    ## speed          3.932      0.416    9.46  1.5e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 15.4 on 48 degrees of freedom
    ## Multiple R-squared:  0.651,  Adjusted R-squared:  0.644 
    ## F-statistic: 89.6 on 1 and 48 DF,  p-value: 1.49e-12

``` r
betahat <- coef(lm.cars)
Vbetahat <- vcov(lm.cars)

K <- diag(2)
Sigma <- diag(1 / sqrt(diag(K %*% Vbetahat %*% t(K)))) 
z <- Sigma %*% K %*% betahat
Cor <- Sigma %*% (K %*% Vbetahat %*% t(K)) %*% t(Sigma)                  

df.cars <- nrow(cars) - length(betahat)
sapply(abs(z), function(x) 1 - pmvt(-rep(x, 2), rep(x, 2), corr = Cor, df = df.cars))
```

    ## [1] 1.661e-02 2.458e-12

``` r
rownames(K) <- names(betahat)

cars.ht <- glht(lm.cars, linfct = K)
summary(cars.ht)
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Fit: lm(formula = dist ~ speed, data = cars)
    ## 
    ## Linear Hypotheses:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) == 0  -17.579      6.758   -2.60    0.017 *  
    ## speed == 0          3.932      0.416    9.46   <1e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

agricolae
---------

``` r
### R code from vignette source 'tutorial.Rnw'
library(agricolae)

A<-as.data.frame(data(package="agricolae")$results[,3:4])
A[,2]<-paste(substr(A[,2],1,35),"..",sep=".")

weight<-c( 68, 53, 69.5, 55, 71, 63, 76.5, 65.5, 69, 75, 76, 57, 70.5, 71.5, 56, 81.5,
           69, 59, 67.5, 61, 68, 59.5, 56.5, 73, 61, 72.5, 71.5, 59.5, 74.5, 63)

par(mfrow=c(1,2),mar=c(4,4,0,1),cex=0.6)
h1<- graph.freq(weight,col=colors()[84],frequency=1,las=2,density=20,ylim=c(0,12),ylab="Frequency") 
x<-h1$breaks
h2<- plot(h1, frequency =2, axes= FALSE,ylim=c(0,0.4),xlab="weight",ylab="Relative (%)")
polygon.freq(h2, col=colors()[84], lwd=2, frequency =2)
axis(1,x,cex=0.6,las=2)
y<-seq(0,0.4,0.1)
axis(2, y,y*100,cex=0.6,las=1) 
```

![](regression_test_rmarkdown_files/figure-markdown_github/agricolae_test-1.png)

Other
=====

Transabyss
----------

The following chunk tests the transabyss installation to be sure that it runs correctly

``` bash
cd $REGRESSION_OUT
cp -r /opt/share/transabyss_sample_dataset .
chmod -R u+w transabyss_sample_dataset/
bash transabyss_sample_dataset/assemble.sh
```

    ## Trans-ABySS 2.0.1
    ## CMD: /opt/bin/transabyss -k 25 --se ./reads/rnaseq_1.fq.gz ./reads/rnaseq_2.fq.gz --outdir ./test.k25 --name test.k25 --threads 2 --island 0 -c 1
    ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ## Found Trans-ABySS directory at: /opt/bin
    ## Found Trans-ABySS `bin` directory at: /opt/bin/bin
    ## Found script at: /opt/bin/bin/skip_psl_self.awk
    ## Found script at: /opt/bin/bin/skip_psl_self_ss.awk
    ## Found `blat' at /opt/bin/blat
    ## Found `MergeContigs' at /usr/lib/abyss/MergeContigs
    ## Found `abyss-filtergraph' at /usr/lib/abyss/abyss-filtergraph
    ## Found `abyss-map' at /usr/lib/abyss/abyss-map
    ## Found `abyss-pe' at /usr/bin/abyss-pe
    ## Found `abyss-junction' at /usr/lib/abyss/abyss-junction
    ## # CPU(s) available:  16
    ## # thread(s) requested:   2
    ## # thread(s) to use:  2
    ## Creating output directory: /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25
    ## CMD: bash -euo pipefail -c 'abyss-pe graph=adj --directory=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25 k=25 name=test.k25 E=0 e=1 c=1 j=2 test.k25-1.fa test.k25-1.adj q=3 se="/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_1.fq.gz /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_2.fq.gz"'
    ## make: Entering directory '/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25'
    ## ABYSS -k25 -q3 -e1 -E0 -c1    --coverage-hist=coverage.hist -s test.k25-bubbles.fa  -o test.k25-1.fa  /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_1.fq.gz /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_2.fq.gz
    ## ABySS 2.0.2
    ## ABYSS -k25 -q3 -e1 -E0 -c1 --coverage-hist=coverage.hist -s test.k25-bubbles.fa -o test.k25-1.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_1.fq.gz /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_2.fq.gz
    ## Reading `/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_1.fq.gz'...
    ## Reading `/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_2.fq.gz'...
    ## Loaded 6114 k-mer
    ## Minimum k-mer coverage is 5
    ## Using a coverage threshold of 16...
    ## The median k-mer coverage is 254
    ## The reconstruction is 3508
    ## The k-mer coverage threshold is 15.9374
    ## Generating adjacency
    ## Added 12338 edges.
    ## Eroding tips
    ## Eroded 0 tips.
    ## Eroded 0 tips.
    ## Pruning tips shorter than 1 bp...
    ## Pruned 4 k-mer in 4 tips.
    ## Pruning tips shorter than 2 bp...
    ## Pruned 2 k-mer in 1 tips.
    ## Pruning tips shorter than 4 bp...
    ## Pruned 37 k-mer in 10 tips.
    ## Pruning tips shorter than 8 bp...
    ## Pruned 105 k-mer in 15 tips.
    ## Pruning tips shorter than 16 bp...
    ## Pruned 345 k-mer in 28 tips.
    ## Pruning tips shorter than 25 bp...
    ## Pruned 769 k-mer in 39 tips.
    ## Pruning tips shorter than 25 bp...
    ## Pruned 97 tips in 6 rounds.
    ## Marked 234 edges of 116 ambiguous vertices.
    ## Removing low-coverage contigs (mean k-mer coverage < 1)
    ## Found 4852 k-mer in 161 contigs before removing low-coverage contigs.
    ## Removed 0 k-mer in 0 low-coverage contigs.
    ## Split 0 ambigiuous branches.
    ## Eroding tips
    ## Eroded 0 tips.
    ## Eroded 0 tips.
    ## Pruning tips shorter than 1 bp...
    ## Pruning tips shorter than 2 bp...
    ## Pruning tips shorter than 4 bp...
    ## Pruning tips shorter than 8 bp...
    ## Pruning tips shorter than 16 bp...
    ## Pruning tips shorter than 25 bp...
    ## Pruned 0 tips in 5 rounds.
    ## Popping bubbles
    ## Removed 32 bubbles.
    ## Removed 32 bubbles
    ## Marked 70 edges of 34 ambiguous vertices.
    ## Assembled 3827 k-mer in 38 contigs.
    ## Removed 1262 k-mer.
    ## The signal-to-noise ratio (SNR) is 5.84861 dB.
    ## AdjList    -k25 -m50 --adj test.k25-1.fa >test.k25-1.adj
    ## make: Leaving directory '/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25'
    ## CHECKPOINT: De Bruijn graph assembly completed.
    ## Iteration 1 of graph simplification ...
    ## ADJ: 76 vertices, 104 edges
    ## Walked 3 paths and marked 12 vertices for removal.
    ## CMD: bash -euo pipefail -c 'abyss-filtergraph --shim --assemble --kmer=25 --island=0 --remove=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.braid.cids --graph=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.filtered.adj1 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-1.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-1.fa > /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.path'
    ## CMD: bash -euo pipefail -c 'MergeContigs --kmer=25 --out=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.filtered.fa --adj --graph=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.filtered.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-1.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.filtered.adj1 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.path'
    ## The minimum coverage of single-end contigs is 1.
    ## The minimum coverage of merged contigs is 210.747.
    ## Consider increasing the coverage threshold parameter, c, to 210.747.
    ## Completed iteration 1 of graph simplification.
    ## Iteration 2 of graph simplification ...
    ## ADJ: 36 vertices, 50 edges
    ## Walked 2 paths and marked 2 vertices for removal.
    ## CMD: bash -euo pipefail -c 'abyss-filtergraph --shim --assemble --kmer=25 --island=0 --remove=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r2.braid.cids --graph=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r2.filtered.adj1 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.filtered.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.filtered.fa > /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r2.path'
    ## CMD: bash -euo pipefail -c 'MergeContigs --kmer=25 --out=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r2.filtered.fa --adj --graph=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r2.filtered.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r1.filtered.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r2.filtered.adj1 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-unitigs.r2.path'
    ## Completed iteration 2 of graph simplification.
    ## Graph simplification stopped at iteration 3
    ## CHECKPOINT: Unitig assembly completed.
    ## ADJ: 34 vertices, 44 edges
    ## Walked 2 paths.
    ## CMD: bash -euo pipefail -c 'MergeContigs --kmer=25 --out=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-ref.fa --merged /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-3.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-3.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-ref.path'
    ## The minimum coverage of single-end contigs is 1.
    ## The minimum coverage of merged contigs is 223.167.
    ## Consider increasing the coverage threshold parameter, c, to 223.167.
    ## CHECKPOINT: Reference path assembly completed.
    ## CMD: bash -euo pipefail -c 'abyss-junction /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-3.adj  >/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-jn.path'
    ## CMD: bash -euo pipefail -c 'MergeContigs --kmer=25 --out=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-jn.fa --merged /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-3.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-3.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-jn.path'
    ## The minimum coverage of single-end contigs is 1.
    ## The minimum coverage of merged contigs is 1.
    ## CHECKPOINT: Junction extension completed.
    ## Using abyss-map to remove redundancy ...
    ## CMD: bash -euo pipefail -c 'abyss-map --dup --threads=2 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-concat.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-concat.fa > /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-concat.fa.dup_ids'
    ## Building the suffix array...
    ## Building the Burrows-Wheeler transform...
    ## Building the character occurrence table...
    ## Removing sequences shorter than 100 ...
    ## CHECKPOINT: Final assembly completed.
    ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ## Assembly generated with Trans-ABySS 2.0.1 :)
    ## Final assembly: /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-final.fa
    ## Total wallclock run time: 0 h 0 m 0 s
    ## Trans-ABySS 2.0.1
    ## CMD: /opt/bin/transabyss -k 32 --se ./reads/rnaseq_1.fq.gz ./reads/rnaseq_2.fq.gz --outdir ./test.k32 --name test.k32 --threads 2 --island 0 -c 1
    ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ## Found Trans-ABySS directory at: /opt/bin
    ## Found Trans-ABySS `bin` directory at: /opt/bin/bin
    ## Found script at: /opt/bin/bin/skip_psl_self.awk
    ## Found script at: /opt/bin/bin/skip_psl_self_ss.awk
    ## Found `blat' at /opt/bin/blat
    ## Found `MergeContigs' at /usr/lib/abyss/MergeContigs
    ## Found `abyss-filtergraph' at /usr/lib/abyss/abyss-filtergraph
    ## Found `abyss-map' at /usr/lib/abyss/abyss-map
    ## Found `abyss-pe' at /usr/bin/abyss-pe
    ## Found `abyss-junction' at /usr/lib/abyss/abyss-junction
    ## # CPU(s) available:  16
    ## # thread(s) requested:   2
    ## # thread(s) to use:  2
    ## Creating output directory: /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32
    ## CMD: bash -euo pipefail -c 'abyss-pe graph=adj --directory=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32 k=32 name=test.k32 E=0 e=1 c=1 j=2 test.k32-1.fa test.k32-1.adj q=3 se="/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_1.fq.gz /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_2.fq.gz"'
    ## make: Entering directory '/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32'
    ## ABYSS -k32 -q3 -e1 -E0 -c1    --coverage-hist=coverage.hist -s test.k32-bubbles.fa  -o test.k32-1.fa  /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_1.fq.gz /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_2.fq.gz
    ## ABySS 2.0.2
    ## ABYSS -k32 -q3 -e1 -E0 -c1 --coverage-hist=coverage.hist -s test.k32-bubbles.fa -o test.k32-1.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_1.fq.gz /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_2.fq.gz
    ## Reading `/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_1.fq.gz'...
    ## Reading `/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/reads/rnaseq_2.fq.gz'...
    ## Loaded 6361 k-mer
    ## Minimum k-mer coverage is 5
    ## Using a coverage threshold of 15...
    ## The median k-mer coverage is 219
    ## The reconstruction is 3499
    ## The k-mer coverage threshold is 14.7986
    ## Generating adjacency
    ## Added 12766 edges.
    ## Eroding tips
    ## Eroded 0 tips.
    ## Eroded 0 tips.
    ## Pruning tips shorter than 1 bp...
    ## Pruned 4 k-mer in 4 tips.
    ## Pruning tips shorter than 2 bp...
    ## Pruned 2 k-mer in 1 tips.
    ## Pruning tips shorter than 4 bp...
    ## Pruned 37 k-mer in 10 tips.
    ## Pruning tips shorter than 8 bp...
    ## Pruned 98 k-mer in 14 tips.
    ## Pruning tips shorter than 16 bp...
    ## Pruned 344 k-mer in 28 tips.
    ## Pruning tips shorter than 32 bp...
    ## Pruned 1606 k-mer in 68 tips.
    ## Pruning tips shorter than 32 bp...
    ## Pruned 125 tips in 6 rounds.
    ## Marked 108 edges of 54 ambiguous vertices.
    ## Removing low-coverage contigs (mean k-mer coverage < 1)
    ## Found 4270 k-mer in 78 contigs before removing low-coverage contigs.
    ## Removed 0 k-mer in 0 low-coverage contigs.
    ## Split 0 ambigiuous branches.
    ## Eroding tips
    ## Eroded 0 tips.
    ## Eroded 0 tips.
    ## Pruning tips shorter than 1 bp...
    ## Pruning tips shorter than 2 bp...
    ## Pruning tips shorter than 4 bp...
    ## Pruning tips shorter than 8 bp...
    ## Pruning tips shorter than 16 bp...
    ## Pruning tips shorter than 32 bp...
    ## Pruned 0 tips in 5 rounds.
    ## Popping bubbles
    ## Removed 17 bubbles.
    ## Removed 17 bubbles
    ## Marked 32 edges of 16 ambiguous vertices.
    ## Assembled 3662 k-mer in 21 contigs.
    ## Removed 2091 k-mer.
    ## The signal-to-noise ratio (SNR) is 3.10074 dB.
    ## AdjList    -k32 -m50 --adj test.k32-1.fa >test.k32-1.adj
    ## make: Leaving directory '/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32'
    ## CHECKPOINT: De Bruijn graph assembly completed.
    ## Iteration 1 of graph simplification ...
    ## ADJ: 42 vertices, 48 edges
    ## Walked 3 paths and marked 8 vertices for removal.
    ## CMD: bash -euo pipefail -c 'abyss-filtergraph --shim --assemble --kmer=32 --island=0 --remove=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.braid.cids --graph=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.filtered.adj1 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-1.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-1.fa > /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.path'
    ## CMD: bash -euo pipefail -c 'MergeContigs --kmer=32 --out=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.filtered.fa --adj --graph=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.filtered.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-1.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.filtered.adj1 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.path'
    ## The minimum coverage of single-end contigs is 1.
    ## The minimum coverage of merged contigs is 194.48.
    ## Consider increasing the coverage threshold parameter, c, to 194.48.
    ## Completed iteration 1 of graph simplification.
    ## Iteration 2 of graph simplification ...
    ## ADJ: 22 vertices, 24 edges
    ## Walked 2 paths and marked 0 vertices for removal.
    ## CMD: bash -euo pipefail -c 'abyss-filtergraph --shim --assemble --kmer=32 --island=0 --remove=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r2.braid.cids --graph=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r2.filtered.adj1 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.filtered.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.filtered.fa > /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r2.path'
    ## CMD: bash -euo pipefail -c 'MergeContigs --kmer=32 --out=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r2.filtered.fa --adj --graph=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r2.filtered.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r1.filtered.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r2.filtered.adj1 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-unitigs.r2.path'
    ## Completed iteration 2 of graph simplification.
    ## Graph simplification stopped at iteration 3
    ## CHECKPOINT: Unitig assembly completed.
    ## ADJ: 22 vertices, 24 edges
    ## Walked 2 paths.
    ## CMD: bash -euo pipefail -c 'MergeContigs --kmer=32 --out=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-ref.fa --merged /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-3.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-3.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-ref.path'
    ## The minimum coverage of single-end contigs is 1.
    ## The minimum coverage of merged contigs is 176.
    ## Consider increasing the coverage threshold parameter, c, to 176.
    ## CHECKPOINT: Reference path assembly completed.
    ## CMD: bash -euo pipefail -c 'abyss-junction /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-3.adj  >/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-jn.path'
    ## CMD: bash -euo pipefail -c 'MergeContigs --kmer=32 --out=/home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-jn.fa --merged /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-3.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-3.adj /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-jn.path'
    ## The minimum coverage of single-end contigs is 1.
    ## The minimum coverage of merged contigs is 1.
    ## CHECKPOINT: Junction extension completed.
    ## Using abyss-map to remove redundancy ...
    ## CMD: bash -euo pipefail -c 'abyss-map --dup --threads=2 /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-concat.fa /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-concat.fa > /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-concat.fa.dup_ids'
    ## Building the suffix array...
    ## Building the Burrows-Wheeler transform...
    ## Building the character occurrence table...
    ## Removing sequences shorter than 100 ...
    ## CHECKPOINT: Final assembly completed.
    ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ## Assembly generated with Trans-ABySS 2.0.1 :)
    ## Final assembly: /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-final.fa
    ## Total wallclock run time: 0 h 0 m 0 s
    ## Trans-ABySS 2.0.1
    ## CMD: /opt/bin/transabyss-merge --mink 25 --maxk 32 --prefixes k25. k32. --out ./merged.fa ./test.k25/test.k25-final.fa ./test.k32/test.k32-final.fa
    ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ## Found Trans-ABySS directory at: /opt/bin
    ## Found Trans-ABySS `bin` directory at: /opt/bin/bin
    ## Found script at: /opt/bin/bin/skip_psl_self.awk
    ## Found script at: /opt/bin/bin/skip_psl_self_ss.awk
    ## Found `abyss-map' at /usr/lib/abyss/abyss-map
    ## Found `blat' at /opt/bin/blat
    ## # CPU(s) available:  16
    ## # thread(s) requested:   1
    ## # thread(s) to use:  1
    ## Assembly Prefixes:
    ## k25. /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k25/test.k25-final.fa
    ## k32. /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/test.k32/test.k32-final.fa
    ## Loaded 1531 letters in 2 sequences
    ## Searched 1531 bases in 2 sequences
    ## Loaded 6800 letters in 4 sequences
    ## Searched 1531 bases in 2 sequences
    ## Loaded 2066 letters in 2 sequences
    ## Searched 1531 bases in 2 sequences
    ## Loaded 2332 letters in 2 sequences
    ## Searched 1531 bases in 2 sequences
    ## Loaded 1531 letters in 2 sequences
    ## Searched 1319 bases in 3 sequences
    ## Loaded 1531 letters in 2 sequences
    ## Searched 1372 bases in 2 sequences
    ## Loaded 1151 letters in 1 sequences
    ## Searched 1531 bases in 2 sequences
    ## Loaded 6800 letters in 4 sequences
    ## Searched 6800 bases in 4 sequences
    ## Loaded 6800 letters in 4 sequences
    ## Searched 2066 bases in 2 sequences
    ## Loaded 6800 letters in 4 sequences
    ## Searched 2332 bases in 2 sequences
    ## Loaded 6800 letters in 4 sequences
    ## Searched 1319 bases in 3 sequences
    ## Loaded 6800 letters in 4 sequences
    ## Searched 1372 bases in 2 sequences
    ## Loaded 6800 letters in 4 sequences
    ## Searched 1151 bases in 1 sequences
    ## Loaded 2066 letters in 2 sequences
    ## Searched 2066 bases in 2 sequences
    ## Loaded 2332 letters in 2 sequences
    ## Searched 2066 bases in 2 sequences
    ## Loaded 2066 letters in 2 sequences
    ## Searched 1319 bases in 3 sequences
    ## Loaded 2066 letters in 2 sequences
    ## Searched 1372 bases in 2 sequences
    ## Loaded 1151 letters in 1 sequences
    ## Searched 2066 bases in 2 sequences
    ## Loaded 2332 letters in 2 sequences
    ## Searched 2332 bases in 2 sequences
    ## Loaded 2332 letters in 2 sequences
    ## Searched 1319 bases in 3 sequences
    ## Loaded 2332 letters in 2 sequences
    ## Searched 1372 bases in 2 sequences
    ## Loaded 2332 letters in 2 sequences
    ## Searched 1151 bases in 1 sequences
    ## Loaded 1319 letters in 3 sequences
    ## Searched 1319 bases in 3 sequences
    ## Loaded 1372 letters in 2 sequences
    ## Searched 1319 bases in 3 sequences
    ## Loaded 1151 letters in 1 sequences
    ## Searched 1319 bases in 3 sequences
    ## Loaded 1372 letters in 2 sequences
    ## Searched 1372 bases in 2 sequences
    ## Loaded 1151 letters in 1 sequences
    ## Searched 1372 bases in 2 sequences
    ## Loaded 1151 letters in 1 sequences
    ## Searched 1151 bases in 1 sequences
    ## =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ## Merged assembly: /home/guest/IBIEM_2018_2019/regression_scratch/transabyss_sample_dataset/merged.fa
    ## Total wallclock run time: 0 h 0 m 1 s

``` r
dir_delete(regression_output_dir)
```

Session Info
============

``` r
sessionInfo()
```

    ## R version 3.4.4 (2018-03-15)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
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
    ## [1] stats     graphics  grDevices utils     datasets  methods  
    ## [7] base     
    ## 
    ## other attached packages:
    ##  [1] agricolae_1.3-0 multcomp_1.4-8  TH.data_1.0-10 
    ##  [4] MASS_7.3-49     survival_2.41-3 mvtnorm_1.0-8  
    ##  [7] phyloseq_1.22.3 dada2_1.6.0     Rcpp_1.0.0     
    ## [10] forcats_0.4.0   stringr_1.3.0   purrr_0.3.1    
    ## [13] readr_1.3.1     tidyr_0.8.0     tibble_2.0.1   
    ## [16] ggplot2_2.2.1   tidyverse_1.2.1 magrittr_1.5   
    ## [19] rmarkdown_1.11  fs_1.2.6        dplyr_0.8.0.1  
    ## [22] here_0.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.0               backports_1.1.1           
    ##   [3] plyr_1.8.4                 igraph_1.2.4              
    ##   [5] lazyeval_0.2.1             sp_1.3-1                  
    ##   [7] splines_3.4.4              BiocParallel_1.12.0       
    ##   [9] AlgDesign_1.1-7.3          GenomeInfoDb_1.14.0       
    ##  [11] digest_0.6.15              foreach_1.4.4             
    ##  [13] htmltools_0.3.6            gdata_2.18.0              
    ##  [15] cluster_2.0.6              Biostrings_2.46.0         
    ##  [17] modelr_0.1.4               RcppParallel_4.4.2        
    ##  [19] matrixStats_0.54.0         gmodels_2.18.1            
    ##  [21] sandwich_2.5-0             colorspace_1.3-2          
    ##  [23] rvest_0.3.2                haven_2.1.0               
    ##  [25] xfun_0.5                   crayon_1.3.4              
    ##  [27] RCurl_1.95-4.12            jsonlite_1.5              
    ##  [29] zoo_1.8-4                  iterators_1.0.10          
    ##  [31] ape_5.2                    glue_1.3.1                
    ##  [33] gtable_0.2.0               zlibbioc_1.24.0           
    ##  [35] XVector_0.18.0             DelayedArray_0.4.1        
    ##  [37] questionr_0.7.0            BiocGenerics_0.24.0       
    ##  [39] scales_1.0.0               DBI_1.0.0                 
    ##  [41] miniUI_0.1.1.1             xtable_1.8-3              
    ##  [43] spData_0.3.0               units_0.6-2               
    ##  [45] spdep_1.0-2                stats4_3.4.4              
    ##  [47] httr_1.4.0                 RColorBrewer_1.1-2        
    ##  [49] pkgconfig_2.0.2            deldir_0.1-16             
    ##  [51] tidyselect_0.2.5           labeling_0.3              
    ##  [53] rlang_0.3.1                reshape2_1.4.2            
    ##  [55] later_0.8.0                munsell_0.5.0             
    ##  [57] cellranger_1.1.0           tools_3.4.4               
    ##  [59] cli_1.0.1                  generics_0.0.2            
    ##  [61] ade4_1.7-13                broom_0.5.1               
    ##  [63] evaluate_0.10.1            biomformat_1.6.0          
    ##  [65] yaml_2.2.0                 knitr_1.22                
    ##  [67] nlme_3.1-131               mime_0.5                  
    ##  [69] xml2_1.2.0                 compiler_3.4.4            
    ##  [71] rstudioapi_0.7             e1071_1.7-0.1             
    ##  [73] klaR_0.6-14                stringi_1.1.6             
    ##  [75] highr_0.6                  lattice_0.20-35           
    ##  [77] Matrix_1.2-12              classInt_0.3-1            
    ##  [79] vegan_2.5-4                permute_0.9-5             
    ##  [81] multtest_2.34.0            pillar_1.3.1              
    ##  [83] LearnBayes_2.15.1          combinat_0.0-8            
    ##  [85] data.table_1.10.4-3        bitops_1.0-6              
    ##  [87] httpuv_1.4.5.1             GenomicRanges_1.30.3      
    ##  [89] R6_2.2.2                   latticeExtra_0.6-28       
    ##  [91] hwriter_1.3.2              promises_1.0.1            
    ##  [93] ShortRead_1.36.1           IRanges_2.12.0            
    ##  [95] codetools_0.2-15           boot_1.3-20               
    ##  [97] gtools_3.8.1               assertthat_0.2.0          
    ##  [99] rhdf5_2.22.0               SummarizedExperiment_1.8.1
    ## [101] rprojroot_1.3-2            GenomicAlignments_1.14.2  
    ## [103] Rsamtools_1.30.0           S4Vectors_0.16.0          
    ## [105] GenomeInfoDbData_1.0.0     mgcv_1.8-23               
    ## [107] expm_0.999-3               parallel_3.4.4            
    ## [109] hms_0.4.2                  grid_3.4.4                
    ## [111] coda_0.19-2                class_7.3-14              
    ## [113] sf_0.7-3                   Biobase_2.38.0            
    ## [115] shiny_1.2.0                lubridate_1.7.4
