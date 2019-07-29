Replace "/sharedspace/platypus\_poop/raw\_data" with the path for your raw data Replace "SRR390728" with your accession number

``` r
raw.data="/sharedspace/platypus_poop/raw_data"
Sys.setenv(RAW_DATA=raw.data)
Sys.setenv(SRA_ACC="SRR390728")
dir.create(raw.data, recursive = TRUE)
```

Quick test of fastq-dump. Should output the first 3 reads from SRR390728

``` bash
fastq-dump --maxSpotId 3 --stdout  $SRA_ACC
```

Check if `--origfmt` gives more useful header Compare this to above

``` bash
fastq-dump --maxSpotId 3 --stdout --origfmt $SRA_ACC
```

Do a test run: download the first five reads and save to compressed FASTQ

``` bash
fastq-dump --maxSpotId 5 --split-files --gzip --outdir $RAW_DATA $SRA_ACC
ls -ltr $RAW_DATA
```

Download the full dataset

``` bash
fastq-dump --split-files --gzip --outdir $RAW_DATA $SRA_ACC
ls -ltr $RAW_DATA
```

Don't forget to set to read-only

``` bash
chomd -R a-w $RAW_DATA
ls -ltrd $RAW_DATA
```
