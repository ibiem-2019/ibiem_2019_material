Atacama Data
============

Atacama data is from <https://docs.qiime2.org/2018.11/tutorials/atacama-soils/#subsample-data>

Setup
=====

First we load libraries.

``` r
library(readr)
library(fs)
```

Paths, Directories, and Shell Variables
---------------------------------------

To keep the code readable and portable, it is nice to assign paths to variables. We also need to use the R `Sys.setenv` command to make shell variables that can be used in the bash chunks below.

``` r
# Directories
data.dir = "/data/tutorial_data/atacama_1pct"
output.dir = path.expand("~/scratch/atacama_1pct")
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

Check Data Integrity
--------------------

``` bash
cd $RAW_FASTQ_DIR
md5sum -c md5sum.txt
```

    ## barcodes.fastq.gz: OK
    ## forward.fastq.gz: OK
    ## reverse.fastq.gz: OK
    ## sample_metadata.tsv: OK

Assemble Metadata Table (Map)
=============================

You are in luck, because you there is a metadata table already made for you. Let's check it out \#\# Examine Metadata Table (Map)

``` r
meta.df = read_tsv(map.file)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
head(meta.df)
```

    ## # A tibble: 6 x 23
    ##   `#SampleID` BarcodeSequence LinkerPrimerSeq… Elevation ExtractConcen
    ##   <chr>       <chr>           <chr>            <chr>     <chr>        
    ## 1 #q2:types   categorical     categorical      numeric   numeric      
    ## 2 BAQ1370.1.2 GCCCAAGTTCAC    CCGGACTACHVGGGT… 1370      0.019        
    ## 3 BAQ1370.3   GCGCCGAATCTT    CCGGACTACHVGGGT… 1370      0.124        
    ## 4 BAQ1370.1.3 ATAAAGAGGAGG    CCGGACTACHVGGGT… 1370      1.200        
    ## 5 BAQ1552.1.1 ATCCCAGCATGC    CCGGACTACHVGGGT… 1552      0.722        
    ## 6 BAQ1552.2   GCTTCCAGACAA    CCGGACTACHVGGGT… 1552      0.017        
    ## # … with 18 more variables: AmpliconConcentration <chr>,
    ## #   ExtractGroupNo <chr>, TransectName <chr>, SiteName <chr>, Depth <chr>,
    ## #   pH <chr>, TOC <chr>, EC <chr>, AverageSoilRelativeHumidity <chr>,
    ## #   RelativeHumiditySoilHigh <chr>, RelativeHumiditySoilLow <chr>,
    ## #   PercentRelativeHumiditySoil_100 <chr>, AverageSoilTemperature <chr>,
    ## #   TemperatureSoilHigh <chr>, TemperatureSoilLow <chr>, Vegetation <chr>,
    ## #   PercentCover <chr>, Description <chr>

Check Map File
--------------

QIIME is inflexible about map file formatting. Fortunately, QIIME includes the [validate\_mapping\_file.py](http://qiime.org/scripts/validate_mapping_file.html) script that checks your map file to see if the format meets its specifications. Unfortunately the script is not very robust, so incorrectly formated map files sometimes make it crash without giving a useful error message. Let's run it anyway . . .

``` bash
validate_mapping_file.py -m $MAP_FILE -o $OUT_DIR/validate_mapfile
```

    ## Errors and/or warnings detected in mapping file.  Please check the log and html file for details.

Once you have run `validate_mapping_file.py` you can view the report through RStudio:

1.  In the *Files* pane, Navigate to /home/guest/scratch/atacama\_1pct/validate\_mapfile
2.  Click on `sample_metadata.tsv.html` and select *View in Web Browser*
3.  Look around the output! How does it look?

Demultiplexing
==============

Remember from the map file that our dataset consists of reads from multiple samples. In order to tally the number of different bacteria corresponding to each sample, we need to separate reads according to the sample of origin.

We will be using `split_libraries_fastq.py` and `split_sequence_file_on_sample_ids.py` to demultiplex the data. The documentation for these programs is available here: - [split\_libraries\_fastq.py](http://qiime.org/scripts/split_libraries_fastq.html) - [split\_sequence\_file\_on\_sample\_ids.py](http://qiime.org/scripts/split_sequence_file_on_sample_ids.html)

Alternatively we can get some instructions if we run the programs with the `--help` command line option.

``` bash
split_libraries_fastq.py --help
```

    ## Usage: split_libraries_fastq.py [options] {-i/--sequence_read_fps SEQUENCE_READ_FPS -o/--output_dir OUTPUT_DIR}
    ## 
    ## [] indicates optional input (order unimportant)
    ## {} indicates required input (order unimportant)
    ## 
    ## 
    ## 
    ## Example usage: 
    ## Print help message and exit
    ##  split_libraries_fastq.py -h
    ## 
    ## Demultiplex and quality filter (at Phred >= Q20) one lane of Illumina fastq data and write results to ./slout_q20.: 
    ##  split_libraries_fastq.py -i lane1_read1.fastq.gz -b lane1_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt -q 19
    ## 
    ## Demultiplex and quality filter (at Phred >= Q20) one lane of Illumina fastq data and write results to ./slout_q20. Store trimmed quality scores in addition to sequence data.: 
    ##  split_libraries_fastq.py -i lane1_read1.fastq.gz -b lane1_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt --store_qual_scores -q 19
    ## 
    ## Demultiplex and quality filter (at Phred >= Q20) two lanes of Illumina fastq data and write results to ./slout_q20.: 
    ##  split_libraries_fastq.py -i lane1_read1.fastq.gz,lane2_read1.fastq.gz -b lane1_barcode.fastq.gz,lane2_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt,map.txt --store_qual_scores -q 19
    ## 
    ## Quality filter (at Phred >= Q20) one non-multiplexed lane of Illumina fastq data and write results to ./slout_single_sample_q20.: 
    ##  split_libraries_fastq.py -i lane1_read1.fastq.gz --sample_ids my.sample.1 -o slout_single_sample_q20/ -q 19 --barcode_type 'not-barcoded'
    ## 
    ## Quality filter (at Phred >= Q20) two non-multiplexed lanes of Illumina fastq data with different samples in each and write results to ./slout_not_multiplexed_q20.: 
    ##  split_libraries_fastq.py -i lane1_read1.fastq.gz,lane2_read1.fastq.gz --sample_ids my.sample.1,my.sample.2 -o slout_not_multiplexed_q20/ -q 19 --barcode_type 'not-barcoded'
    ## 
    ## Options:
    ##   --version             show program's version number and exit
    ##   -h, --help            show this help message and exit
    ##   -v, --verbose         Print information during execution -- useful for
    ##                         debugging [default: False]
    ##   -m MAPPING_FPS, --mapping_fps=MAPPING_FPS
    ##                         metadata mapping files (comma-separated if more than
    ##                         one) [default: none]
    ##   -b BARCODE_READ_FPS, --barcode_read_fps=BARCODE_READ_FPS
    ##                         the barcode read fastq files (comma-separated if more
    ##                         than one) [default: none]
    ##   --store_qual_scores   store qual strings in .qual files [default: False]
    ##   --sample_ids=SAMPLE_IDS
    ##                         comma-separated list of samples ids to be applied to
    ##                         all sequences, must be one per input file path (used
    ##                         when data is not multiplexed) [default: none]
    ##   --store_demultiplexed_fastq
    ##                         write demultiplexed fastq files [default: False]
    ##   --retain_unassigned_reads
    ##                         retain sequences which don't map to a barcode in the
    ##                         mapping file (sample ID will be "Unassigned")
    ##                         [default: False]
    ##   -r MAX_BAD_RUN_LENGTH, --max_bad_run_length=MAX_BAD_RUN_LENGTH
    ##                         max number of consecutive low quality base calls
    ##                         allowed before truncating a read [default: 3]
    ##   -p MIN_PER_READ_LENGTH_FRACTION, --min_per_read_length_fraction=MIN_PER_READ_LENGTH_FRACTION
    ##                         min number of consecutive high quality base calls to
    ##                         include a read (per single end read) as a fraction of
    ##                         the input read length [default: 0.75]
    ##   -n SEQUENCE_MAX_N, --sequence_max_n=SEQUENCE_MAX_N
    ##                         maximum number of N characters allowed in a sequence
    ##                         to retain it -- this is applied after quality
    ##                         trimming, and is total over combined paired end reads
    ##                         if applicable [default: 0]
    ##   -s START_SEQ_ID, --start_seq_id=START_SEQ_ID
    ##                         start seq_ids as ascending integers beginning with
    ##                         start_seq_id [default: 0]
    ##   --rev_comp_barcode    reverse complement barcode reads before lookup
    ##                         [default: False]
    ##   --rev_comp_mapping_barcodes
    ##                         reverse complement barcode in mapping before lookup
    ##                         (useful if barcodes in mapping file are reverse
    ##                         complements of golay codes) [default: False]
    ##   --rev_comp            reverse complement sequence before writing to output
    ##                         file (useful for reverse-orientation reads) [default:
    ##                         False]
    ##   -q PHRED_QUALITY_THRESHOLD, --phred_quality_threshold=PHRED_QUALITY_THRESHOLD
    ##                         the maximum unacceptable Phred quality score (e.g.,
    ##                         for Q20 and better, specify -q 19) [default: 3]
    ##   --last_bad_quality_char=LAST_BAD_QUALITY_CHAR
    ##                         DEPRECATED: use -q instead. This method of setting is
    ##                         not robust to different versions of CASAVA.
    ##   --barcode_type=BARCODE_TYPE
    ##                         The type of barcode used. This can be an integer, e.g.
    ##                         for length 6 barcodes, or "golay_12" for golay error-
    ##                         correcting barcodes. Error correction will only be
    ##                         applied for "golay_12" barcodes. If data is not
    ##                         barcoded, pass "not-barcoded". [default: golay_12]
    ##   --max_barcode_errors=MAX_BARCODE_ERRORS
    ##                         maximum number of errors in barcode [default: 1.5]
    ##   --phred_offset=PHRED_OFFSET
    ##                         the ascii offset to use when decoding phred scores
    ##                         (either 33 or 64). Warning: in most cases you don't
    ##                         need to pass this value [default: determined
    ##                         automatically]
    ## 
    ##   REQUIRED options:
    ##     The following options must be provided under all circumstances.
    ## 
    ##     -i SEQUENCE_READ_FPS, --sequence_read_fps=SEQUENCE_READ_FPS
    ##                         the sequence read fastq files (comma-separated if more
    ##                         than one) [REQUIRED]
    ##     -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
    ##                         directory to store output files [REQUIRED]

Running split\_libraries\_fastq.py
----------------------------------

We are now ready to run `split_libraries_fastq.py`. Here is an explanation for the command line options that we use

-   Quality Parameters: `split_libraries_fastq.py` can do quality filtering, but we don't want it to because we will be doing quality filtering in DADA2, so we will set these paramters to the least stringent values possible to be sure that nothing gets filtered
    -   -r, --max\_bad\_run\_length
    -   -p, --min\_per\_read\_length\_fraction
    -   -n, --sequence\_max\_n
    -   -q, --phred\_quality\_threshold
-   --sequence\_read\_fps : this is the FASTQ that is being demuxtiplexed, we have to run `split_libraries_fastq.py` separately on the R1 and R2 files
-   --output\_dir : the output directory; there is more than one output file
-   --barcode\_read\_fps : the barcode FASTQ
-   --mapping\_fps : the file mapping samples to barcodes (and other metadata)
-   --phred\_offset : are the FASTQs phred33 or phred66?
-   --barcode\_type : Are the barcodes EMP golay codes? How long are they
-   --store\_demultiplexed\_fastq : save demultiplexed FASTQs (default is to only generate a FASTA because that's what the qiime pipeline uses after this point)
-   --retain\_unassigned\_reads : be a packrat - don't throw out reads that don't match a barcode

``` bash
set -u
TAGDIR=$DEMUX_DIR/tagged_1
split_libraries_fastq.py -r 999 -n 999 -q 0 -p 0.0001 \
        --sequence_read_fps $RAW_FASTQ_DIR/forward.fastq.gz \
        --output_dir $TAGDIR \
        --barcode_read_fps $BARCODE_FASTQ \
        --mapping_fps $MAP_FILE \
        --phred_offset 33 \
        --barcode_type golay_12 \
        --store_demultiplexed_fastq \
        --retain_unassigned_reads
```

    ## Error in split_libraries_fastq.py: Some or all barcodes are not valid golay codes. Do they need to be reverse complemented? If these are not golay barcodes pass --barcode_type 12 to disable barcode error correction, or pass --barcode_type # if the barcodes are not 12 base pairs, where # is the size of the barcodes. Invalid codes:
    ##  CTGGTGCTGAAT AAGATCGTACTG CATCTGGGCAAT GCGCCGAATCTT CGCATACGACCT GACCGTCAATAC TTATGTACGGCG TCCTAGGTCCGA TCAGACCAACTG AACCCAGATGAT AGACAAGCTTCC TTGGAACGGCTT TTCCACACGTGG GGAAATCCCATC ACTCGGCCAACT GCATGCATCCCA TGGAGAGGAGAT GGAAGAAGTAGC TTGGACGTCCAC TACGCCCATCAG TTCACCTGTATC TCTCGATAAGCG AACCTCGGATAA CTCCAGGTCATG GTAGTGTCAACA GTAGCACTCATG CACGAGCTACTC GTGGCCTACTAC CACCTGTAGTAG GCCTCGTACTGA TGTCCGTGGATC ATAAAGAGGAGG ACTCATCTTCCA AATCTTGCGCCG CGTATAAATGCG TTCCCTTCTCCG ACGTAACCACGT ATTCAGATGGCA GTCTCCTCCCTT TCCTCACTATCA GTGCTTGTGTAG GCTTCCAGACAA GAGATACAGTTC ATTATACGGCGC ATCGATCCACAG TCACTTGGTGCG TCCGCAACCTGA TCCAGGGCTATA GCCTGCAGTACT ATCCCAGCATGC TGTCAGCTGTCG TAAACGCGACTC AAGTGAAGCGAG CAACTAGACTCG CCTCGGGTACTA ACACCGCACAAT GATCTAATCGAG ACCAACAGATTG GTCGGAAATTGT GAAACTCCTAGA AATACAGACCTG GGAACGACGTGA CATTTGACGACG AGTGCCCTTGGT ATCGGGCTTAAC GCGTAGAGAGAC GACAGAGGTGCA TCTAACGAGTGC CAGGATTCGTAC GACTCAACCAGT TGCCGCCGTAAT GCCCAAGTTCAC GTTGGTTGGCAT ACACAGTCCTGA TAGACACCGTGT
    ## 
    ## If you need help with QIIME, see:
    ## http://help.qiime.org

Hmm, an error message! Let's read it carefully to see what it says.

It is saying that some of the barcodes are not valid golay codes. We know that these barcodes are golay barcodes (as is common with 16s barcodes), so that can't be right.

It says that we could disable barcode error correction with the `--barcode_type 12` option, but that's a bit like disabling the brakes in our car because they squeal.

It asks if they "need to be reverse complemented". There are two options for doing that: - --rev\_comp\_barcode: Reverse complement barcode reads before lookup \[default: False\] - --rev\_comp\_mapping\_barcodes: Reverse complement barcode in mapping before lookup (useful if barcodes in mapping file are reverse complements of golay codes) \[default: False\]

Let's try one

``` bash
set -u
TAGDIR=$DEMUX_DIR/tagged_2
split_libraries_fastq.py -r 999 -n 999 -q 0 -p 0.0001 \
        --sequence_read_fps $RAW_FASTQ_DIR/forward.fastq.gz \
        --output_dir $TAGDIR \
        --barcode_read_fps $BARCODE_FASTQ \
        --mapping_fps $MAP_FILE \
        --phred_offset 33 \
        --barcode_type golay_12 \
        --rev_comp_barcode \
        --store_demultiplexed_fastq \
        --retain_unassigned_reads
```

    ## Error in split_libraries_fastq.py: Some or all barcodes are not valid golay codes. Do they need to be reverse complemented? If these are not golay barcodes pass --barcode_type 12 to disable barcode error correction, or pass --barcode_type # if the barcodes are not 12 base pairs, where # is the size of the barcodes. Invalid codes:
    ##  CTGGTGCTGAAT AAGATCGTACTG CATCTGGGCAAT GCGCCGAATCTT CGCATACGACCT GACCGTCAATAC TTATGTACGGCG TCCTAGGTCCGA TCAGACCAACTG AACCCAGATGAT AGACAAGCTTCC TTGGAACGGCTT TTCCACACGTGG GGAAATCCCATC ACTCGGCCAACT GCATGCATCCCA TGGAGAGGAGAT GGAAGAAGTAGC TTGGACGTCCAC TACGCCCATCAG TTCACCTGTATC TCTCGATAAGCG AACCTCGGATAA CTCCAGGTCATG GTAGTGTCAACA GTAGCACTCATG CACGAGCTACTC GTGGCCTACTAC CACCTGTAGTAG GCCTCGTACTGA TGTCCGTGGATC ATAAAGAGGAGG ACTCATCTTCCA AATCTTGCGCCG CGTATAAATGCG TTCCCTTCTCCG ACGTAACCACGT ATTCAGATGGCA GTCTCCTCCCTT TCCTCACTATCA GTGCTTGTGTAG GCTTCCAGACAA GAGATACAGTTC ATTATACGGCGC ATCGATCCACAG TCACTTGGTGCG TCCGCAACCTGA TCCAGGGCTATA GCCTGCAGTACT ATCCCAGCATGC TGTCAGCTGTCG TAAACGCGACTC AAGTGAAGCGAG CAACTAGACTCG CCTCGGGTACTA ACACCGCACAAT GATCTAATCGAG ACCAACAGATTG GTCGGAAATTGT GAAACTCCTAGA AATACAGACCTG GGAACGACGTGA CATTTGACGACG AGTGCCCTTGGT ATCGGGCTTAAC GCGTAGAGAGAC GACAGAGGTGCA TCTAACGAGTGC CAGGATTCGTAC GACTCAACCAGT TGCCGCCGTAAT GCCCAAGTTCAC GTTGGTTGGCAT ACACAGTCCTGA TAGACACCGTGT
    ## 
    ## If you need help with QIIME, see:
    ## http://help.qiime.org

That didn't work, let's try both

``` bash
set -u
TAGDIR=$DEMUX_DIR/tagged_3
split_libraries_fastq.py -r 999 -n 999 -q 0 -p 0.0001 \
        --sequence_read_fps $RAW_FASTQ_DIR/forward.fastq.gz \
        --output_dir $TAGDIR \
        --barcode_read_fps $BARCODE_FASTQ \
        --mapping_fps $MAP_FILE \
        --phred_offset 33 \
        --barcode_type golay_12 \
        --rev_comp_barcode \
        --rev_comp_mapping_barcodes \
        --store_demultiplexed_fastq \
        --retain_unassigned_reads
```

That's better, let's check the output:

``` bash
ls $DEMUX_DIR/tagged_3/
```

    ## histograms.txt
    ## seqs.fastq
    ## seqs.fna
    ## split_library_log.txt

``` bash
cat $DEMUX_DIR/tagged_3/histograms.txt
```

    ## Length   Count
    ## 151.0    2712
    ## --

``` bash
cat $DEMUX_DIR/tagged_3/split_library_log.txt
```

    ## Input file paths
    ## Mapping filepath: /data/tutorial_data/atacama_1pct/sample_metadata.tsv (md5: 3f7eda5424236e66300a7a503ac7926e)
    ## Sequence read filepath: /data/tutorial_data/atacama_1pct/forward.fastq.gz (md5: b5ee4718605158e5db91ef906fb61f62)
    ## Barcode read filepath: /data/tutorial_data/atacama_1pct/barcodes.fastq.gz (md5: df32df8dc627319d7bc1d9ce49b89849)
    ## 
    ## Quality filter results
    ## Total number of input sequences: 135487
    ## Barcode not in mapping file: 0
    ## Read too short after quality truncation: 0
    ## Count of N characters exceeds limit: 0
    ## Illumina quality digit = 0: 0
    ## Barcode errors exceed max: 132775
    ## 
    ## Result summary (after quality filtering)
    ## Median sequence length: 151.00
    ## Unassigned   2711
    ## YUN3533.1.1  1
    ## YUN3856.3    0
    ## YUN3856.2    0
    ## YUN3856.1.3  0
    ## YUN3856.1.2  0
    ## YUN3856.1.1  0
    ## YUN3533.3    0
    ## YUN3533.2    0
    ## YUN3533.1.3  0
    ## YUN3533.1.2  0
    ## YUN3428.3    0
    ## YUN3428.2    0
    ## YUN3428.1    0
    ## YUN3346.3    0
    ## YUN3346.2    0
    ## YUN3346.1    0
    ## YUN3259.3    0
    ## YUN3259.2    0
    ## YUN3259.1.3  0
    ## YUN3259.1.2  0
    ## YUN3259.1.1  0
    ## YUN3184.2    0
    ## YUN3184.1    0
    ## YUN3153.3    0
    ## YUN3153.2    0
    ## YUN3153.1    0
    ## YUN3008.3    0
    ## YUN3008.2    0
    ## YUN3008.1.3  0
    ## YUN3008.1.2  0
    ## YUN3008.1.1  0
    ## YUN2029.3    0
    ## YUN2029.2    0
    ## YUN2029.1    0
    ## YUN1609.3    0
    ## YUN1609.1    0
    ## YUN1242.3    0
    ## YUN1242.2    0
    ## YUN1242.1    0
    ## YUN1005.3    0
    ## YUN1005.2    0
    ## YUN1005.1.3  0
    ## YUN1005.1.1  0
    ## BAQ895.3 0
    ## BAQ895.2 0
    ## BAQ4697.3    0
    ## BAQ4697.2    0
    ## BAQ4697.1    0
    ## BAQ4166.3    0
    ## BAQ4166.2    0
    ## BAQ4166.1.3  0
    ## BAQ4166.1.2  0
    ## BAQ4166.1.1  0
    ## BAQ3473.3    0
    ## BAQ3473.2    0
    ## BAQ3473.1    0
    ## BAQ2838.3    0
    ## BAQ2838.2    0
    ## BAQ2838.1    0
    ## BAQ2687.3    0
    ## BAQ2687.2    0
    ## BAQ2687.1    0
    ## BAQ2462.3    0
    ## BAQ2462.2    0
    ## BAQ2462.1    0
    ## BAQ2420.3    0
    ## BAQ2420.2    0
    ## BAQ2420.1.3  0
    ## BAQ2420.1.2  0
    ## BAQ2420.1.1  0
    ## BAQ1552.2    0
    ## BAQ1552.1.1  0
    ## BAQ1370.3    0
    ## BAQ1370.1.3  0
    ## BAQ1370.1.2  0
    ## 
    ## Total number seqs written    2712
    ## ---

`split_library_log.txt` gives us some summary statistics: total number of reads processed, number of reads that fail various quality tests, number of reads assigned to each sample (based on their barcode), and total number of reads that were assigned to all barcodes.

This doesn't look so good!

"Total number of input sequences: 135487"

So we are inputting 135487 reads into `split_libraries_fastq.py` (Is this what you expected?), but only 2712 are assigned to a barcode (confusingly 2711 are assigned to the "Unassigned" barcode, so only one read is assigned to a sample). We are losing almost all of our reds!

The quality filtering looks good, except "Barcode errors exceed max: 132775"; that's where the most of the reads are going (98%). This could be a problem with our data, but it suggests that there might be a problem with our analysis, particularly how we are handling the barcodes. Perhaps we still haven't figured it out yet! Let's try this:

``` bash
set -u
TAGDIR=$DEMUX_DIR/tagged_4
split_libraries_fastq.py -r 999 -n 999 -q 0 -p 0.0001 \
        --sequence_read_fps $RAW_FASTQ_DIR/forward.fastq.gz \
        --output_dir $TAGDIR \
        --barcode_read_fps $BARCODE_FASTQ \
        --mapping_fps $MAP_FILE \
        --phred_offset 33 \
        --barcode_type golay_12 \
        --rev_comp_mapping_barcodes \
        --store_demultiplexed_fastq \
        --retain_unassigned_reads
```

``` bash
cat $DEMUX_DIR/tagged_4/split_library_log.txt
```

    ## Input file paths
    ## Mapping filepath: /data/tutorial_data/atacama_1pct/sample_metadata.tsv (md5: 3f7eda5424236e66300a7a503ac7926e)
    ## Sequence read filepath: /data/tutorial_data/atacama_1pct/forward.fastq.gz (md5: b5ee4718605158e5db91ef906fb61f62)
    ## Barcode read filepath: /data/tutorial_data/atacama_1pct/barcodes.fastq.gz (md5: df32df8dc627319d7bc1d9ce49b89849)
    ## 
    ## Quality filter results
    ## Total number of input sequences: 135487
    ## Barcode not in mapping file: 0
    ## Read too short after quality truncation: 0
    ## Count of N characters exceeds limit: 0
    ## Illumina quality digit = 0: 0
    ## Barcode errors exceed max: 28688
    ## 
    ## Result summary (after quality filtering)
    ## Median sequence length: 151.00
    ## Unassigned   62978
    ## YUN3428.2    1379
    ## BAQ4166.1.2  1304
    ## BAQ4166.2    1264
    ## YUN3533.3    1253
    ## BAQ4166.3    1236
    ## YUN3856.3    1233
    ## BAQ3473.3    1173
    ## YUN3533.2    1159
    ## BAQ3473.1    1097
    ## BAQ4697.3    1084
    ## YUN3259.3    1081
    ## YUN3428.1    1066
    ## BAQ4166.1.3  1052
    ## BAQ4166.1.1  1038
    ## YUN3428.3    1014
    ## YUN3259.2    1011
    ## YUN3856.2    983
    ## YUN3856.1.2  971
    ## YUN3856.1.1  952
    ## YUN3533.1.1  937
    ## YUN3533.1.2  886
    ## BAQ2687.3    880
    ## YUN1005.1.1  873
    ## BAQ2687.1    869
    ## BAQ2462.1    869
    ## YUN2029.2    834
    ## BAQ3473.2    824
    ## BAQ4697.1    820
    ## YUN3346.3    793
    ## YUN3533.1.3  765
    ## BAQ4697.2    761
    ## YUN1242.3    735
    ## YUN1609.1    729
    ## BAQ2420.1.1  721
    ## YUN3259.1.2  683
    ## BAQ2687.2    659
    ## BAQ2420.1.3  656
    ## YUN3346.1    652
    ## YUN3856.1.3  650
    ## BAQ2420.3    647
    ## BAQ2420.1.2  640
    ## BAQ2420.2    611
    ## BAQ2462.2    598
    ## BAQ2838.1    588
    ## YUN3153.3    570
    ## YUN1242.1    536
    ## BAQ2838.2    489
    ## BAQ2462.3    467
    ## YUN3153.2    444
    ## YUN1005.3    404
    ## BAQ2838.3    369
    ## YUN3259.1.3  272
    ## YUN3346.2    171
    ## YUN3259.1.1  56
    ## YUN2029.1    6
    ## YUN3184.2    1
    ## YUN3008.3    1
    ## YUN3008.1.3  1
    ## YUN2029.3    1
    ## YUN1242.2    1
    ## BAQ1552.1.1  1
    ## BAQ1370.1.3  1
    ## YUN3184.1    0
    ## YUN3153.1    0
    ## YUN3008.2    0
    ## YUN3008.1.2  0
    ## YUN3008.1.1  0
    ## YUN1609.3    0
    ## YUN1005.2    0
    ## YUN1005.1.3  0
    ## BAQ895.3 0
    ## BAQ895.2 0
    ## BAQ1552.2    0
    ## BAQ1370.3    0
    ## BAQ1370.1.2  0
    ## 
    ## Total number seqs written    106799
    ## ---

We still have a bunch of Barcode errors but many fewer (17% instead of 98%). Many are still Unassigned, but most of our samples have some reads, and many of our samples have a large number. It looks like we needed to reverse complement both the barcodes as supplied in the map file and the barcodes as sequenced.

Running `split_sequence_file_on_sample_ids.py`
----------------------------------------------

Despite its name `split_libraries_fastq.py` does not actually *spilt* the FASTQ, it just relabels or "tags" it. To actually do the demultiplexing we need another program: `split_sequence_file_on_sample_ids.py`. Fortunately the commands for `split_sequence_file_on_sample_ids.py` are a little simpler.

``` bash
split_sequence_file_on_sample_ids.py -h
```

    ## Usage: split_sequence_file_on_sample_ids.py [options] {-i/--input_seqs_fp INPUT_SEQS_FP -o/--output_dir OUTPUT_DIR}
    ## 
    ## [] indicates optional input (order unimportant)
    ## {} indicates required input (order unimportant)
    ## 
    ## Split a single post-split_libraries.py fasta (or post-split_libraries_fastq.py fastq) file into per-sample fasta files. This script requires that the sequences identitifers are in post-split_libraries.py format (i.e., SampleID_SeqID). A file will be created for each unique SampleID.
    ## 
    ## Example usage: 
    ## Print help message and exit
    ##  split_sequence_file_on_sample_ids.py -h
    ## 
    ## Split seqs.fna into one fasta file per sample and store the resulting fasta files in 'out'
    ##  split_sequence_file_on_sample_ids.py -i seqs.fna -o out/
    ## 
    ## Split seqs.fastq into one fastq file per sample and store the resulting fastq files in 'out_fastq'
    ##  split_sequence_file_on_sample_ids.py -i seqs.fastq --file_type fastq -o out_fastq/
    ## 
    ## Options:
    ##   --version             show program's version number and exit
    ##   -h, --help            show this help message and exit
    ##   -v, --verbose         Print information during execution -- useful for
    ##                         debugging [default: False]
    ##   --buffer_size=BUFFER_SIZE
    ##                         the number of sequences to read into memory before
    ##                         writing to file (you usually won't need to change
    ##                         this) [default: 500]
    ##   --file_type=FILE_TYPE
    ##                         Type of file. Either fasta or fastq
    ## 
    ##   REQUIRED options:
    ##     The following options must be provided under all circumstances.
    ## 
    ##     -i INPUT_SEQS_FP, --input_seqs_fp=INPUT_SEQS_FP
    ##                         the input fasta file to split [REQUIRED]
    ##     -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
    ##                         the output directory [default: none] [REQUIRED]

Here's what we will use:

-   --input\_seqs\_fp: the "tagged" fastq
-   --file\_type: FASTA or FASTQ?
-   --output\_dir: where to put the demuxed FASTQs

``` bash
TAGDIR=$DEMUX_DIR/tagged_4
SPLITDIR=$DEMUX_DIR/split_4
split_sequence_file_on_sample_ids.py --input_seqs_fp $TAGDIR/seqs.fastq \
                     --file_type fastq \
                     --output_dir $SPLITDIR
```

Now let's check that it worked

``` bash
ls -lSrh $DEMUX_DIR/split_4
```

    ## total 44M
    ## -rw-r--r-- 1 guest users  427 Apr  1 14:41 YUN1242.2.fastq
    ## -rw-r--r-- 1 guest users  428 Apr  1 14:41 YUN3184.2.fastq
    ## -rw-r--r-- 1 guest users  428 Apr  1 14:41 YUN3008.3.fastq
    ## -rw-r--r-- 1 guest users  428 Apr  1 14:41 YUN2029.3.fastq
    ## -rw-r--r-- 1 guest users  429 Apr  1 14:41 BAQ1370.1.3.fastq
    ## -rw-r--r-- 1 guest users  430 Apr  1 14:41 YUN3008.1.3.fastq
    ## -rw-r--r-- 1 guest users  430 Apr  1 14:41 BAQ1552.1.1.fastq
    ## -rw-r--r-- 1 guest users 2.6K Apr  1 14:41 YUN2029.1.fastq
    ## -rw-r--r-- 1 guest users  24K Apr  1 14:41 YUN3259.1.1.fastq
    ## -rw-r--r-- 1 guest users  72K Apr  1 14:41 YUN3346.2.fastq
    ## -rw-r--r-- 1 guest users 115K Apr  1 14:41 YUN3259.1.3.fastq
    ## -rw-r--r-- 1 guest users 155K Apr  1 14:41 BAQ2838.3.fastq
    ## -rw-r--r-- 1 guest users 169K Apr  1 14:41 YUN1005.3.fastq
    ## -rw-r--r-- 1 guest users 186K Apr  1 14:41 YUN3153.2.fastq
    ## -rw-r--r-- 1 guest users 195K Apr  1 14:41 BAQ2462.3.fastq
    ## -rw-r--r-- 1 guest users 205K Apr  1 14:41 BAQ2838.2.fastq
    ## -rw-r--r-- 1 guest users 224K Apr  1 14:41 YUN1242.1.fastq
    ## -rw-r--r-- 1 guest users 238K Apr  1 14:41 YUN3153.3.fastq
    ## -rw-r--r-- 1 guest users 246K Apr  1 14:41 BAQ2838.1.fastq
    ## -rw-r--r-- 1 guest users 250K Apr  1 14:41 BAQ2462.2.fastq
    ## -rw-r--r-- 1 guest users 256K Apr  1 14:41 BAQ2420.2.fastq
    ## -rw-r--r-- 1 guest users 269K Apr  1 14:41 BAQ2420.1.2.fastq
    ## -rw-r--r-- 1 guest users 271K Apr  1 14:41 BAQ2420.3.fastq
    ## -rw-r--r-- 1 guest users 273K Apr  1 14:41 YUN3346.1.fastq
    ## -rw-r--r-- 1 guest users 273K Apr  1 14:41 YUN3856.1.3.fastq
    ## -rw-r--r-- 1 guest users 276K Apr  1 14:41 BAQ2687.2.fastq
    ## -rw-r--r-- 1 guest users 276K Apr  1 14:41 BAQ2420.1.3.fastq
    ## -rw-r--r-- 1 guest users 287K Apr  1 14:41 YUN3259.1.2.fastq
    ## -rw-r--r-- 1 guest users 303K Apr  1 14:41 BAQ2420.1.1.fastq
    ## -rw-r--r-- 1 guest users 305K Apr  1 14:41 YUN1609.1.fastq
    ## -rw-r--r-- 1 guest users 307K Apr  1 14:41 YUN1242.3.fastq
    ## -rw-r--r-- 1 guest users 318K Apr  1 14:41 BAQ4697.2.fastq
    ## -rw-r--r-- 1 guest users 321K Apr  1 14:41 YUN3533.1.3.fastq
    ## -rw-r--r-- 1 guest users 332K Apr  1 14:41 YUN3346.3.fastq
    ## -rw-r--r-- 1 guest users 343K Apr  1 14:41 BAQ4697.1.fastq
    ## -rw-r--r-- 1 guest users 344K Apr  1 14:41 BAQ3473.2.fastq
    ## -rw-r--r-- 1 guest users 349K Apr  1 14:41 YUN2029.2.fastq
    ## -rw-r--r-- 1 guest users 363K Apr  1 14:41 BAQ2462.1.fastq
    ## -rw-r--r-- 1 guest users 363K Apr  1 14:41 BAQ2687.1.fastq
    ## -rw-r--r-- 1 guest users 367K Apr  1 14:41 YUN1005.1.1.fastq
    ## -rw-r--r-- 1 guest users 368K Apr  1 14:41 BAQ2687.3.fastq
    ## -rw-r--r-- 1 guest users 372K Apr  1 14:41 YUN3533.1.2.fastq
    ## -rw-r--r-- 1 guest users 393K Apr  1 14:41 YUN3533.1.1.fastq
    ## -rw-r--r-- 1 guest users 400K Apr  1 14:41 YUN3856.1.1.fastq
    ## -rw-r--r-- 1 guest users 408K Apr  1 14:41 YUN3856.1.2.fastq
    ## -rw-r--r-- 1 guest users 411K Apr  1 14:41 YUN3856.2.fastq
    ## -rw-r--r-- 1 guest users 423K Apr  1 14:41 YUN3259.2.fastq
    ## -rw-r--r-- 1 guest users 424K Apr  1 14:41 YUN3428.3.fastq
    ## -rw-r--r-- 1 guest users 436K Apr  1 14:41 BAQ4166.1.1.fastq
    ## -rw-r--r-- 1 guest users 442K Apr  1 14:41 BAQ4166.1.3.fastq
    ## -rw-r--r-- 1 guest users 445K Apr  1 14:41 YUN3428.1.fastq
    ## -rw-r--r-- 1 guest users 452K Apr  1 14:41 YUN3259.3.fastq
    ## -rw-r--r-- 1 guest users 453K Apr  1 14:41 BAQ4697.3.fastq
    ## -rw-r--r-- 1 guest users 458K Apr  1 14:41 BAQ3473.1.fastq
    ## -rw-r--r-- 1 guest users 484K Apr  1 14:41 YUN3533.2.fastq
    ## -rw-r--r-- 1 guest users 490K Apr  1 14:41 BAQ3473.3.fastq
    ## -rw-r--r-- 1 guest users 515K Apr  1 14:41 YUN3856.3.fastq
    ## -rw-r--r-- 1 guest users 516K Apr  1 14:41 BAQ4166.3.fastq
    ## -rw-r--r-- 1 guest users 524K Apr  1 14:41 YUN3533.3.fastq
    ## -rw-r--r-- 1 guest users 528K Apr  1 14:41 BAQ4166.2.fastq
    ## -rw-r--r-- 1 guest users 547K Apr  1 14:41 BAQ4166.1.2.fastq
    ## -rw-r--r-- 1 guest users 576K Apr  1 14:41 YUN3428.2.fastq
    ## -rw-r--r-- 1 guest users  26M Apr  1 14:41 Unassigned.fastq

Looks good like we generated a demultiplexed FASTQ for each sample!

### Putting it together for R1 and R2

This will run `split_libraries_fastq.py` and `split_sequence_file_on_sample_ids.py` on both R1 and R2, and do a little cleanup (get rid of the results of `split_libraries_fastq.py` once we have demuxed it. We can drop "--retain\_unassigned\_reads" since we have already reviewed the results.

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

``` bash
ls $DEMUX_DIR/
```

    ## forward
    ## reverse
    ## split_4
    ## tagged_1
    ## tagged_2
    ## tagged_3
    ## tagged_4

``` bash
ls $DEMUX_DIR/forward $DEMUX_DIR/reverse
```

    ## /home/guest/scratch/atacama_1pct/demux/forward:
    ## BAQ1370.1.3.fastq
    ## BAQ1552.1.1.fastq
    ## BAQ2420.1.1.fastq
    ## BAQ2420.1.2.fastq
    ## BAQ2420.1.3.fastq
    ## BAQ2420.2.fastq
    ## BAQ2420.3.fastq
    ## BAQ2462.1.fastq
    ## BAQ2462.2.fastq
    ## BAQ2462.3.fastq
    ## BAQ2687.1.fastq
    ## BAQ2687.2.fastq
    ## BAQ2687.3.fastq
    ## BAQ2838.1.fastq
    ## BAQ2838.2.fastq
    ## BAQ2838.3.fastq
    ## BAQ3473.1.fastq
    ## BAQ3473.2.fastq
    ## BAQ3473.3.fastq
    ## BAQ4166.1.1.fastq
    ## BAQ4166.1.2.fastq
    ## BAQ4166.1.3.fastq
    ## BAQ4166.2.fastq
    ## BAQ4166.3.fastq
    ## BAQ4697.1.fastq
    ## BAQ4697.2.fastq
    ## BAQ4697.3.fastq
    ## YUN1005.1.1.fastq
    ## YUN1005.3.fastq
    ## YUN1242.1.fastq
    ## YUN1242.2.fastq
    ## YUN1242.3.fastq
    ## YUN1609.1.fastq
    ## YUN2029.1.fastq
    ## YUN2029.2.fastq
    ## YUN2029.3.fastq
    ## YUN3008.1.3.fastq
    ## YUN3008.3.fastq
    ## YUN3153.2.fastq
    ## YUN3153.3.fastq
    ## YUN3184.2.fastq
    ## YUN3259.1.1.fastq
    ## YUN3259.1.2.fastq
    ## YUN3259.1.3.fastq
    ## YUN3259.2.fastq
    ## YUN3259.3.fastq
    ## YUN3346.1.fastq
    ## YUN3346.2.fastq
    ## YUN3346.3.fastq
    ## YUN3428.1.fastq
    ## YUN3428.2.fastq
    ## YUN3428.3.fastq
    ## YUN3533.1.1.fastq
    ## YUN3533.1.2.fastq
    ## YUN3533.1.3.fastq
    ## YUN3533.2.fastq
    ## YUN3533.3.fastq
    ## YUN3856.1.1.fastq
    ## YUN3856.1.2.fastq
    ## YUN3856.1.3.fastq
    ## YUN3856.2.fastq
    ## YUN3856.3.fastq
    ## 
    ## /home/guest/scratch/atacama_1pct/demux/reverse:
    ## BAQ1370.1.3.fastq
    ## BAQ1552.1.1.fastq
    ## BAQ2420.1.1.fastq
    ## BAQ2420.1.2.fastq
    ## BAQ2420.1.3.fastq
    ## BAQ2420.2.fastq
    ## BAQ2420.3.fastq
    ## BAQ2462.1.fastq
    ## BAQ2462.2.fastq
    ## BAQ2462.3.fastq
    ## BAQ2687.1.fastq
    ## BAQ2687.2.fastq
    ## BAQ2687.3.fastq
    ## BAQ2838.1.fastq
    ## BAQ2838.2.fastq
    ## BAQ2838.3.fastq
    ## BAQ3473.1.fastq
    ## BAQ3473.2.fastq
    ## BAQ3473.3.fastq
    ## BAQ4166.1.1.fastq
    ## BAQ4166.1.2.fastq
    ## BAQ4166.1.3.fastq
    ## BAQ4166.2.fastq
    ## BAQ4166.3.fastq
    ## BAQ4697.1.fastq
    ## BAQ4697.2.fastq
    ## BAQ4697.3.fastq
    ## YUN1005.1.1.fastq
    ## YUN1005.3.fastq
    ## YUN1242.1.fastq
    ## YUN1242.2.fastq
    ## YUN1242.3.fastq
    ## YUN1609.1.fastq
    ## YUN2029.1.fastq
    ## YUN2029.2.fastq
    ## YUN2029.3.fastq
    ## YUN3008.1.3.fastq
    ## YUN3008.3.fastq
    ## YUN3153.2.fastq
    ## YUN3153.3.fastq
    ## YUN3184.2.fastq
    ## YUN3259.1.1.fastq
    ## YUN3259.1.2.fastq
    ## YUN3259.1.3.fastq
    ## YUN3259.2.fastq
    ## YUN3259.3.fastq
    ## YUN3346.1.fastq
    ## YUN3346.2.fastq
    ## YUN3346.3.fastq
    ## YUN3428.1.fastq
    ## YUN3428.2.fastq
    ## YUN3428.3.fastq
    ## YUN3533.1.1.fastq
    ## YUN3533.1.2.fastq
    ## YUN3533.1.3.fastq
    ## YUN3533.2.fastq
    ## YUN3533.3.fastq
    ## YUN3856.1.1.fastq
    ## YUN3856.1.2.fastq
    ## YUN3856.1.3.fastq
    ## YUN3856.2.fastq
    ## YUN3856.3.fastq

So the demuxed forward reads are in the `forward` directory and the demuxed reverse reads are in the `reverse` directory. We are ready for DADA2!

Bonus: Rename and move split FASTQs
-----------------------------------

``` r
for (curread in c("forward","reverse")) {
  curpath = file.path(demux.dir,curread)
  print(curpath)
  # cur_fastqs = list.files(curpath, full.names = TRUE,pattern = ".fastq")
  # print(cur_fastqs)
  for (fastq_path in list.files(curpath, full.names = TRUE,pattern = ".fastq")){
    print(fastq_path)
    new_path = path_ext_remove(fastq_path)
    print(new_path)
    new_path = path_file(new_path)
    print(new_path)
    new_path = path(demux.dir, new_path, ext=paste0(curread,".fastq"))
    print(new_path)
    file_move(fastq_path, new_path)
  }
}
```

    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ1370.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ1370.1.3"
    ## BAQ1370.1.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ1370.1.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ1552.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ1552.1.1"
    ## BAQ1552.1.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ1552.1.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.1.1"
    ## BAQ2420.1.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.1.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.1.2"
    ## BAQ2420.1.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.1.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.1.3"
    ## BAQ2420.1.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.1.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.2"
    ## BAQ2420.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2420.3"
    ## BAQ2420.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2462.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2462.1"
    ## BAQ2462.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2462.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2462.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2462.2"
    ## BAQ2462.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2462.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2462.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2462.3"
    ## BAQ2462.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2462.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2687.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2687.1"
    ## BAQ2687.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2687.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2687.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2687.2"
    ## BAQ2687.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2687.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2687.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2687.3"
    ## BAQ2687.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2687.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2838.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2838.1"
    ## BAQ2838.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2838.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2838.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2838.2"
    ## BAQ2838.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2838.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2838.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ2838.3"
    ## BAQ2838.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2838.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ3473.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ3473.1"
    ## BAQ3473.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ3473.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ3473.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ3473.2"
    ## BAQ3473.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ3473.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ3473.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ3473.3"
    ## BAQ3473.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ3473.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.1.1"
    ## BAQ4166.1.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.1.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.1.2"
    ## BAQ4166.1.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.1.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.1.3"
    ## BAQ4166.1.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.1.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.2"
    ## BAQ4166.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4166.3"
    ## BAQ4166.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4697.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4697.1"
    ## BAQ4697.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4697.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4697.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4697.2"
    ## BAQ4697.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4697.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4697.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/BAQ4697.3"
    ## BAQ4697.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4697.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1005.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1005.1.1"
    ## YUN1005.1.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN1005.1.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1005.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1005.3"
    ## YUN1005.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN1005.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1242.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1242.1"
    ## YUN1242.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN1242.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1242.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1242.2"
    ## YUN1242.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN1242.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1242.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1242.3"
    ## YUN1242.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN1242.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1609.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN1609.1"
    ## YUN1609.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN1609.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN2029.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN2029.1"
    ## YUN2029.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN2029.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN2029.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN2029.2"
    ## YUN2029.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN2029.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN2029.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN2029.3"
    ## YUN2029.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN2029.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3008.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3008.1.3"
    ## YUN3008.1.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3008.1.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3008.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3008.3"
    ## YUN3008.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3008.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3153.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3153.2"
    ## YUN3153.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3153.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3153.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3153.3"
    ## YUN3153.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3153.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3184.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3184.2"
    ## YUN3184.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3184.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.1.1"
    ## YUN3259.1.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.1.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.1.2"
    ## YUN3259.1.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.1.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.1.3"
    ## YUN3259.1.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.1.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.2"
    ## YUN3259.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3259.3"
    ## YUN3259.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3346.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3346.1"
    ## YUN3346.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3346.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3346.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3346.2"
    ## YUN3346.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3346.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3346.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3346.3"
    ## YUN3346.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3346.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3428.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3428.1"
    ## YUN3428.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3428.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3428.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3428.2"
    ## YUN3428.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3428.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3428.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3428.3"
    ## YUN3428.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3428.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.1.1"
    ## YUN3533.1.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.1.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.1.2"
    ## YUN3533.1.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.1.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.1.3"
    ## YUN3533.1.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.1.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.2"
    ## YUN3533.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3533.3"
    ## YUN3533.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.1.1"
    ## YUN3856.1.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.1.1.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.1.2"
    ## YUN3856.1.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.1.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.1.3"
    ## YUN3856.1.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.1.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.2"
    ## YUN3856.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.2.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/forward/YUN3856.3"
    ## YUN3856.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.3.forward.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ1370.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ1370.1.3"
    ## BAQ1370.1.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ1370.1.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ1552.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ1552.1.1"
    ## BAQ1552.1.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ1552.1.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.1.1"
    ## BAQ2420.1.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.1.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.1.2"
    ## BAQ2420.1.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.1.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.1.3"
    ## BAQ2420.1.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.1.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.2"
    ## BAQ2420.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2420.3"
    ## BAQ2420.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2420.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2462.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2462.1"
    ## BAQ2462.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2462.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2462.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2462.2"
    ## BAQ2462.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2462.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2462.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2462.3"
    ## BAQ2462.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2462.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2687.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2687.1"
    ## BAQ2687.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2687.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2687.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2687.2"
    ## BAQ2687.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2687.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2687.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2687.3"
    ## BAQ2687.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2687.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2838.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2838.1"
    ## BAQ2838.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2838.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2838.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2838.2"
    ## BAQ2838.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2838.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2838.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ2838.3"
    ## BAQ2838.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ2838.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ3473.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ3473.1"
    ## BAQ3473.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ3473.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ3473.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ3473.2"
    ## BAQ3473.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ3473.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ3473.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ3473.3"
    ## BAQ3473.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ3473.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.1.1"
    ## BAQ4166.1.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.1.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.1.2"
    ## BAQ4166.1.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.1.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.1.3"
    ## BAQ4166.1.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.1.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.2"
    ## BAQ4166.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4166.3"
    ## BAQ4166.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4166.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4697.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4697.1"
    ## BAQ4697.1
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4697.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4697.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4697.2"
    ## BAQ4697.2
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4697.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4697.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/BAQ4697.3"
    ## BAQ4697.3
    ## /home/guest/scratch/atacama_1pct/demux/BAQ4697.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1005.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1005.1.1"
    ## YUN1005.1.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN1005.1.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1005.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1005.3"
    ## YUN1005.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN1005.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1242.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1242.1"
    ## YUN1242.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN1242.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1242.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1242.2"
    ## YUN1242.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN1242.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1242.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1242.3"
    ## YUN1242.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN1242.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1609.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN1609.1"
    ## YUN1609.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN1609.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN2029.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN2029.1"
    ## YUN2029.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN2029.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN2029.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN2029.2"
    ## YUN2029.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN2029.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN2029.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN2029.3"
    ## YUN2029.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN2029.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3008.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3008.1.3"
    ## YUN3008.1.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3008.1.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3008.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3008.3"
    ## YUN3008.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3008.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3153.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3153.2"
    ## YUN3153.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3153.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3153.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3153.3"
    ## YUN3153.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3153.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3184.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3184.2"
    ## YUN3184.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3184.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.1.1"
    ## YUN3259.1.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.1.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.1.2"
    ## YUN3259.1.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.1.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.1.3"
    ## YUN3259.1.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.1.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.2"
    ## YUN3259.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3259.3"
    ## YUN3259.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3259.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3346.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3346.1"
    ## YUN3346.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3346.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3346.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3346.2"
    ## YUN3346.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3346.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3346.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3346.3"
    ## YUN3346.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3346.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3428.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3428.1"
    ## YUN3428.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3428.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3428.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3428.2"
    ## YUN3428.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3428.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3428.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3428.3"
    ## YUN3428.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3428.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.1.1"
    ## YUN3533.1.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.1.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.1.2"
    ## YUN3533.1.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.1.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.1.3"
    ## YUN3533.1.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.1.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.2"
    ## YUN3533.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3533.3"
    ## YUN3533.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3533.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.1.1.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.1.1"
    ## YUN3856.1.1
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.1.1.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.1.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.1.2"
    ## YUN3856.1.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.1.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.1.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.1.3"
    ## YUN3856.1.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.1.3.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.2.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.2"
    ## YUN3856.2
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.2.reverse.fastq
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.3.fastq"
    ## [1] "/home/guest/scratch/atacama_1pct/demux/reverse/YUN3856.3"
    ## YUN3856.3
    ## /home/guest/scratch/atacama_1pct/demux/YUN3856.3.reverse.fastq

Session Info
============

Always print `sessionInfo` for reproducibility!

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] fs_1.2.7    readr_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.1       fansi_0.4.0      assertthat_0.2.0 utf8_1.1.3      
    ##  [5] crayon_1.3.4     digest_0.6.15    R6_2.2.2         magrittr_1.5    
    ##  [9] evaluate_0.13    pillar_1.3.1     cli_1.1.0        rlang_0.3.2     
    ## [13] stringi_1.1.6    rmarkdown_1.12   tools_3.4.4      stringr_1.3.0   
    ## [17] hms_0.4.2        xfun_0.5         yaml_2.2.0       compiler_3.4.4  
    ## [21] pkgconfig_2.0.2  htmltools_0.3.6  knitr_1.22       tibble_2.1.1
