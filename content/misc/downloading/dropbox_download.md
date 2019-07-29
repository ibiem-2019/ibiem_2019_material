Set up directory to receive data
================================

Replace "/sharedspace/platypus\_poop/raw\_data" with the path for your raw data

``` r
raw.data="/sharedspace/platypus_poop/raw_data"
Sys.setenv(RAW_DATA=raw.data)
dir.create(raw.data, recursive = TRUE)
```

Get Dropbox Link
================

Follow the instructions to [Share a link from dropbox.com](https://www.dropbox.com/help/files-folders/view-only-access#link)

Download Data
=============

The link you copied from dropbox should look something like: <https://www.dropbox.com/s/xv3mjz64bypkgbh/2015Miseq.zip?dl=0>

You will need to use `wget` to download within RStudio, either in an R Notebook (preferred) or from the terminal. Here is an example `wget` command.

``` bash
wget --content-disposition --no-verbose \
  --directory-prefix $RAW_DATA \
  https://www.dropbox.com/s/xv3mjz64bypkgbh/2015Miseq.zip?dl=1
```

Note that I changed the "0" to "1" at the end of the dropbox link (see <https://superuser.com/a/486351>).

As to the command line options: - `--content-disposition`: Use filename information from the server. Otherwise the downloaded file will end up with a funny name - `--no-verbose`: Don't output details. Otherwise you will get a lot of output, which is particularly annoying when run in a notebook - `--directory-prefix $RAW_DATA`: Save the file to the `$RAW_DATA` directory

Check to see that the file downloaded
=====================================

``` bash
ls -ltr $RAW_DATA
```

Unzip the downloaded file
=========================

``` bash
unzip -q -d $RAW_DATA $RAW_DATA/2015Miseq.zip
```

Check zip archive contents
==========================

``` bash
ls -ltr $RAW_DATA/2015Miseq/
```

Clean up
========

No need to keep the original zip file around anymore

``` bash
rm $RAW_DATA/2015Miseq.zip
```

Gzip the fastqs
===============

Saves space! Never keep unzipped fastqs around

``` bash
gzip $RAW_DATA/2015Miseq/101515_MiSeq/*.fastq
```

Make the data directory read-only
=================================

``` bash
chmod -R a-w $RAW_DATA/2015Miseq
```
