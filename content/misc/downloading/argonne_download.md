<!-- 
R -e "rmarkdown::render('misc/argonne_download.Rmd', output_format=c('html_document', 'md_document'))" 
R -e "rmarkdown::render('misc/argonne_download.Rmd')" 
-->
> NOTE: This file will not work as-is because the URLs for data download
> are for one-time use.

Load Libraries
==============

``` r
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
```

Set up directory to receive data
================================

``` r
raw.data="/sharedspace/argonne_data"
uncompressed.md5 = file.path(raw.data,"md5_checksum_uncompressed_fastqs.txt")
compressed.md5 = file.path(raw.data,"md5_checksum_compressed_fastqs.txt")

Sys.setenv(RAW_DATA=raw.data)
Sys.setenv(COMPRESSED_MD5=compressed.md5)
Sys.setenv(UNCOMPRESSED_MD5=uncompressed.md5)
dir.create(raw.data, recursive = TRUE, showWarnings = FALSE)
```

Download Data from Argonne
==========================

Download with wget
------------------

These chunks will not work as-is because the one-time URLs have already
been used. New one-time URLs can be obtained by logging in to
<a href="https://sequencing.bio.anl.gov/index.html" class="uri">https://sequencing.bio.anl.gov/index.html</a>

-   `--no-verbose` tells `wget` to *not* print a lot of status
    information during the download.
-   `--content-disposition` tells `wget` to use the file name from the
    server (otherwise it will name files with the random charracters in
    the one-time URLs)
-   `--directory-prefix` tells `wget` what directory to download the
    files into

``` bash
wget --content-disposition --no-verbose \
  --directory-prefix $RAW_DATA \
  http://shock.metagenomics.anl.gov/preauth/pBiK1ZDIkDqGLEibz88c
```

``` bash
wget --content-disposition --no-verbose \
  --directory-prefix $RAW_DATA \
  http://shock.metagenomics.anl.gov/preauth/IUUT80Lu4aOKKkWLJsRl
```

``` bash
wget --content-disposition --no-verbose \
  --directory-prefix $RAW_DATA \
  http://shock.metagenomics.anl.gov/preauth/UKdpX2Qk7nszU790SBDn
```

Make MD5 sum file
-----------------

Annoyingly, the md5sum hashes are provided in the ANL web interface, not
as a downloadable file, we have to manually construct the md5sum file by
copying the md5 hash for each FASTQ from the web interface.

``` r
tribble(
  ~md5sum, ~file,
  "8d1d121829a8519f404026654c26ba3d", "Undetermined_S0_L001_I1_001.fastq",
  "0a0e4475e81fa0cd2ea3d891064c521f", "Undetermined_S0_L001_R1_001.fastq",
  "b1c1503dd7a05e92ad0eb296e66deadb", "Undetermined_S0_L001_R2_001.fastq",
  "d68f457d872e1aaba215f0b7b7585ec0", "200114_McCumber_16SFW_AS_200110.txt"
) %>%
  write_tsv(uncompressed.md5, col_names = FALSE)
```

Confirm Downloaded files
------------------------

Annoyingly, the FASTQs are supplied unzipped, so we need to:

1.  Run `md5sum -c` on the downloaded FASTQs using the manually
    constructed `unzipped_md5sum.txt`
2.  gzip the FASTQs
3.  Generate a new md5sum file for the gzipped FASTQs
4.  It seems like a good idea to chain these together with `&&` so that
    the whole thing fails if one step fails, instead of giving the false
    impression that it worked if an intermediate step fails.

``` bash
cd $RAW_DATA && \
  md5sum -c $UNCOMPRESSED_MD5 && \
  gzip Undetermined_S0_L001_??_001.fastq && \
  md5sum Undetermined_S0_L001_??_001.fastq.gz 200114_McCumber_16SFW_AS_200110.txt > $COMPRESSED_MD5
```

Make directory read-only
------------------------

``` bash
chmod -R a-w $RAW_DATA
```

Archive Data
============

This is only necessary if data isn’t already archived on DDS

Upload data to DDS
------------------

``` bash
ddsclient upload -p IBIEM_2019 $RAW_DATA
```

Confirm DDS upload
------------------

``` bash
ddsclient upload --dry-run -p IBIEM_2019 $RAW_DATA
```
