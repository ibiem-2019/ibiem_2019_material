Download references for taxonomic assignment in DADA2: <https://benjjneb.github.io/dada2/training.html>

Set up directory to receive data
================================

``` r
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
```

Set up directory to receive data
================================

``` r
ref.dir = "/data/references/dada"
md5.file = file.path(ref.dir, "md5sum.txt")

Sys.setenv(REF_DIR=ref.dir)
Sys.setenv(MD5_FILE=md5.file)

dir.create(ref.dir, recursive = TRUE, showWarnings = FALSE)
```

Generate MD5 file from zenodo
=============================

``` r
tribble(
  ~path, ~md5sum,
  "silva_nr_v132_train_set.fa.gz", "2c4e2d8d9a936cdf24a91c0d2c113a43", 
  "silva_species_assignment_v132.fa.gz", "059fa72391d6aa2f17bf69f2cd63b8ea", 
  "silva_nr_v128_train_set.fa.gz", "26b23e13ed310174ae5563e8f7258ecc",
  "silva_species_assignment_v128.fa.gz", "f67d7d9443312ca38a9fbf61fc744ae4",
  "silva_nr_v123_train_set.fa.gz", "005e0f0fd4c8478c8840e7a35647c36d", 
  "silva_species_assignment_v123.fa.gz", "ec82e11498bb32656bd96d3fa8c79c76"
) %>%
  select(md5sum, path) %>%
  write_tsv(md5.file, col_names=FALSE)
```

Download data
=============

``` bash
set -u
wget --no-verbose --directory-prefix $REF_DIR \
  https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz \
  https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz \
  https://zenodo.org/record/824551/files/silva_nr_v128_train_set.fa.gz \
  https://zenodo.org/record/824551/files/silva_species_assignment_v128.fa.gz \
  https://zenodo.org/record/158958/files/silva_nr_v123_train_set.fa.gz \
  https://zenodo.org/record/158958/files/silva_species_assignment_v123.fa.gz
```

Check md5sums
=============

``` bash
cd $REF_DIR
md5sum -c $MD5_FILE
```

Make the data directory read-only
=================================

``` bash
chmod -R a-w $REF_DIR
```
