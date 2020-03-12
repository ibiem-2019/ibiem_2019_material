Setup
=====

Load Libraries
--------------

``` r
library(readr)
library(phyloseq)
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)
library(fs)
library(knitr)
```

Setup Paths
-----------

``` r
# Directories

if(exists("params") && 
   !is.null(params[["atacama_ps_rds"]])){
  atacama.ps.rds=params[["atacama_ps_rds"]]
} else {
  atacama.ps.rds = "/data/tutorial_data/atacama_1pct.rds"
}
print(atacama.ps.rds)
```

    ## [1] "/data/tutorial_data/atacama_1pct.rds"

``` r
# atacama.ps.rds = "/data/tutorial_data/atacama_1pct.rds"
outdir=path.expand("~/scratch/lefse")
# lefse.input.file=file.path(outdir,"data_for_lefse.tsv")

if(fig_path()!="-1"){
  fig_path() %>%
    dirname() -> outdir
}

atacama.ps.rds %>%
  path_file %>%
  path_ext_remove %>%
  path(outdir, .) ->
  output_basename 

lefse.input.file = path_ext_set(output_basename, ".tsv")



if(dir_exists(outdir)){
  dir_delete(outdir)
}
dir.create(outdir)


Sys.setenv(LEFSE_INPUT_FILE=lefse.input.file)
Sys.setenv(OUTDIR=outdir)
# Sys.setenv(FIGURE_DIR=lefse_figure_dir)
Sys.setenv(BASENAME=output_basename)

# NORMALIZATION=""
Sys.setenv(NORMALIZATION="-o 1000000")
Sys.setenv(PLOT_FORMAT="png")

atacama.ps = read_rds(atacama.ps.rds)
print(atacama.ps)
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 459 taxa and 61 samples ]
    ## sample_data() Sample Data:       [ 61 samples by 22 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 459 taxa by 7 taxonomic ranks ]

Run Lefse
=========

Define Functions
----------------

``` r
#' RepseqToTaxa:
#' Convert Repseq column names to Taxa column names in a spread data frame
#' The big problem here is that this needs to be done after all other 
#' manipulations to the dataframe, otherwise most functions will balk if there
#' are dataframe columns with identical names
#'
#' @param spread.df The dataframe generated from phyloseq object.
#' @param source.ps phyloseq object from which spread.df was derived.
RepseqToTaxa <- function(spread.df, source.ps) {
  tax.df = as.data.frame(tax_table(source.ps)) %>%
    rownames_to_column("repseq") %>%
    mutate(taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus, sep="|")) %>%
    select(repseq, taxonomy)
  
  # need to preserve non-OTU column names (otherwise they get lost)
  colname_match = match(names(spread.df), tax.df$repseq)
  cols_to_keep = which(is.na(colname_match))
  colnames_to_keep = names(spread.df)[cols_to_keep]
  
  # replace repseqs with taxonomies
  names(spread.df) = tax.df$taxonomy[match(names(spread.df), tax.df$repseq)]
  # now reset the non-OTU column names
  names(spread.df)[cols_to_keep] = colnames_to_keep
  return(spread.df)
}

GenerateLefseOutput = function(ps,output_columns,outfile){
  if(length(output_columns)==1){
    format.flags="-f c -u 1 -c 2"
  }else if(length(output_columns)==2){
    format.flags="-f c -u 1 -c 2 -s 3"
  }else{
    stop("output_columns must be a vector of length 1 or 2")
  }
  base_columns = c("SampleID","OTU","Abundance")
  spread.df = psmelt(ps) %>% 
    mutate(SampleID = str_replace(Sample, pattern="\\.", replacement="_")) %>%
    mutate(taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus, sep="|")) %>%
    select(one_of(c(base_columns,output_columns))) %>%
    spread(OTU, Abundance)
  
  RepseqToTaxa(spread.df, ps) %>%
    write.table(file=outfile, 
                sep="\t", quote = FALSE,
                row.names = FALSE)
  
  return(format.flags)
}
```

Select variables
----------------

List variables to select parameters of interest

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
# format.flags = GenerateLefseOutput(atacama.ps, c("Vegetation","TransectName"), lefse.input.file)
grouping_parameter="Vegetation"
format.flags = GenerateLefseOutput(atacama.ps, grouping_parameter, lefse.input.file)
Sys.setenv(FORMAT_FLAGS=format.flags)
```

Format data for Lefse
---------------------

``` bash
format_input.py $LEFSE_INPUT_FILE  "${BASENAME}.in" $FORMAT_FLAGS $NORMALIZATION  --output_table ${BASENAME}.tab
```

Run Lefse itself
----------------

``` bash
run_lefse.py "${BASENAME}.in" ${BASENAME}.res
```

    ## Number of significantly discriminative features: 15 ( 15 ) before internal wilcoxon
    ## Number of discriminative features with abs LDA score > 2.0 : 15

Generate bar plot of significant taxa
-------------------------------------

``` bash
plot_res.py --format ${PLOT_FORMAT} ${BASENAME}.res ${BASENAME}.${PLOT_FORMAT}
```

    ## /usr/lib/python2.7/dist-packages/matplotlib/artist.py:896: MatplotlibDeprecationWarning: The set_axis_bgcolor function was deprecated in version 2.0. Use set_facecolor instead.
    ##   return func(v)

``` r
lefse.bar = path_ext_set(output_basename, ".png")
cat(paste0("![Barplot of taxa associated with ", grouping_parameter, "](", lefse.bar, ")"), fill = FALSE)
```

![Barplot of taxa associated with
Vegetation](lefse_files/figure-markdown_github/atacama_1pct.png)

Generate cladogram plot of significant taxa
-------------------------------------------

``` bash
plot_cladogram.py ${BASENAME}.res ${BASENAME}.cladogram.${PLOT_FORMAT} --format ${PLOT_FORMAT}
```

    ## /usr/lib/python2.7/dist-packages/matplotlib/artist.py:896: MatplotlibDeprecationWarning: The set_axis_bgcolor function was deprecated in version 2.0. Use set_facecolor instead.
    ##   return func(v)
    ## clade_sep parameter too large, lowered to 0.200225830078

``` r
lefse.cladogram = path_ext_set(output_basename, ".cladogram.png")

cat(paste0("![Cladogram of taxa associated with ", grouping_parameter, "](", lefse.cladogram, ")"), fill = FALSE)
```

![Cladogram of taxa associated with
Vegetation](lefse_files/figure-markdown_github/atacama_1pct.cladogram.png)

Generate separate plots for each significant taxon
--------------------------------------------------

``` bash
mkdir -p ${BASENAME}_individual_taxon_plots

plot_features.py ${BASENAME}.in ${BASENAME}.res ${BASENAME}_individual_taxon_plots/
```

    ## /usr/lib/python2.7/dist-packages/matplotlib/artist.py:896: MatplotlibDeprecationWarning: The set_axis_bgcolor function was deprecated in version 2.0. Use set_facecolor instead.
    ##   return func(v)
    ## Exporting  Bacteria.Bacteroidetes.Bacteroidia.Chitinophagales
    ## Exporting  Bacteria.Proteobacteria.Gammaproteobacteria.Steroidobacterales.Steroidobacteraceae
    ## Exporting  Bacteria.Proteobacteria.Gammaproteobacteria.Steroidobacterales
    ## Exporting  Bacteria.Actinobacteria.Actinobacteria.Micrococcales
    ## Exporting  Bacteria.Bacteroidetes.Bacteroidia.Chitinophagales.Chitinophagaceae
    ## Exporting  Bacteria.Proteobacteria.Gammaproteobacteria.Methylococcales.Methylococcaceae
    ## Exporting  Bacteria.Actinobacteria.Actinobacteria.Micrococcales.Micrococcaceae
    ## Exporting  Bacteria.Proteobacteria.Gammaproteobacteria.Methylococcales
    ## Exporting  Bacteria.Actinobacteria.Actinobacteria.Micrococcales.Micrococcaceae.Pseudarthrobacter
    ## Exporting  Bacteria.Proteobacteria.Gammaproteobacteria.Methylococcales.Methylococcaceae.NA
    ## Exporting  Bacteria.Gemmatimonadetes.Gemmatimonadetes.Gemmatimonadales.Gemmatimonadaceae.NA
    ## Exporting  Bacteria.Proteobacteria.Deltaproteobacteria.Myxococcales
    ## Exporting  Bacteria.Firmicutes.Bacilli.Bacillales.Paenibacillaceae.Ammoniphilus
    ## Exporting  Bacteria.Proteobacteria.Deltaproteobacteria
    ## Exporting  Bacteria.Proteobacteria.Gammaproteobacteria.Steroidobacterales.Steroidobacteraceae.Steroidobacter

``` r
taxon_plots = paste0(output_basename, "_individual_taxon_plots")

for (curplot in list.files(taxon_plots, full.names = TRUE)){
  cat(paste0("![Individual Taxon Associated With ", grouping_parameter, "](", curplot, ")"), fill = FALSE)
}
```

![Individual Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Actinobacteria-Actinobacteria-Micrococcales-Micrococcaceae-Pseudarthrobacter.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Actinobacteria-Actinobacteria-Micrococcales-Micrococcaceae.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Actinobacteria-Actinobacteria-Micrococcales.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Bacteroidetes-Bacteroidia-Chitinophagales-Chitinophagaceae.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Bacteroidetes-Bacteroidia-Chitinophagales.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Firmicutes-Bacilli-Bacillales-Paenibacillaceae-Ammoniphilus.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Gemmatimonadetes-Gemmatimonadetes-Gemmatimonadales-Gemmatimonadaceae-NA.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Proteobacteria-Deltaproteobacteria-Myxococcales.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Proteobacteria-Deltaproteobacteria.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Proteobacteria-Gammaproteobacteria-Methylococcales-Methylococcaceae-NA.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Proteobacteria-Gammaproteobacteria-Methylococcales-Methylococcaceae.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Proteobacteria-Gammaproteobacteria-Methylococcales.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Proteobacteria-Gammaproteobacteria-Steroidobacterales-Steroidobacteraceae-Steroidobacter.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Proteobacteria-Gammaproteobacteria-Steroidobacterales-Steroidobacteraceae.png)![Individual
Taxon Associated With
Vegetation](lefse_files/figure-markdown_github/atacama_1pct_individual_taxon_plots/1_Bacteria-Proteobacteria-Gammaproteobacteria-Steroidobacterales.png)
