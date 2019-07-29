Start Docker Demo
=================

Overview
========

What is Reproducible Analysis?
------------------------------

"Reproducible Analysis is an important part of reproducible research. Reproducible analysis requires that all components of the analysis be archived so that *anyone* can independently repeat the analysis and arrive at exactly the same results" - Josh Granek

### Who is "anyone"

> -   Labmates
> -   Collaborators
> -   Competitors
> -   Everyone else interested in your work
> -   Someone who wants to apply your analysis to their own data
> -   You, 6 months from now

Not Covered Here
----------------

-   Reproducibility in the Wet Lab
-   Lots of Other Stuff
-   Many Other Valid Approaches
-   Important Details

Soapbox
-------

1.  Setup
2.  Stand on

Reproducible Analysis
---------------------

. . . is like eating your vegetables

Three Pillars of Reproducible Analysis
--------------------------------------

-   Raw Data
-   Compute Environment
-   Analysis Process

### Executive Summary

-   Archive the Necessary Components!

TLDR
====

Raw Data
========

Care and Handling of Raw Data
-----------------------------

1.  Download
2.  Provenance
    -   Verify checksum
    -   Generate checksum if it didn't come with the data
3.  Protection
    -   Set *READ-ONLY* permissions on data files and data directory
    -   Archive ASAP

Switch to RStudio
-----------------

Checksums
---------

``` bash
cd $RAWDATA_DIR
head *.csv
```

    ## ==> data1.csv <==
    ## value
    ## 1
    ## 2
    ## 3
    ## 4
    ## 
    ## ==> data2.csv <==
    ## value
    ## 5
    ## 6
    ## 7
    ## 8

Checksums
---------

``` bash
cd $RAWDATA_DIR
md5sum *.csv
```

    ## 25ef3a11e49406d5f4be5163c0d1ce06  data1.csv
    ## 7e9e9bf71d2aaa2486a6536f3043311e  data2.csv

Checksums File
--------------

``` bash
cd $RAWDATA_DIR
md5sum *.csv > mydata_md5.txt
md5sum -c mydata_md5.txt
```

    ## data1.csv: OK
    ## data2.csv: OK

Catching Changes
----------------

``` bash
sed -i s/7/3/ $DATA2
cat $DATA2
```

    ## value
    ## 5
    ## 6
    ## 3
    ## 8

``` bash
cd $RAWDATA_DIR
md5sum -c mydata_md5.txt
```

    ## data1.csv: OK
    ## data2.csv: FAILED
    ## md5sum: WARNING: 1 computed checksum did NOT match

READ-ONLY
---------

``` bash
ls -ld $RAWDATA_DIR
ls -ltr $RAWDATA_DIR
```

    ## drwxr-xr-x 2 guest users 4096 Aug 10 23:04 /tmp/Rtmpgogiql/rr_raw_data
    ## total 12
    ## -rw-r--r-- 1 guest users 14 Aug 10 23:04 data1.csv
    ## -rw-r--r-- 1 guest users 88 Aug 10 23:04 mydata_md5.txt
    ## -rw-r--r-- 1 guest users 14 Aug 10 23:04 data2.csv

``` bash
chmod -R a-w $RAWDATA_DIR
ls -ld $RAWDATA_DIR
ls -ltr $RAWDATA_DIR
```

    ## dr-xr-xr-x 2 guest users 4096 Aug 10 23:04 /tmp/Rtmpgogiql/rr_raw_data
    ## total 12
    ## -r--r--r-- 1 guest users 14 Aug 10 23:04 data1.csv
    ## -r--r--r-- 1 guest users 88 Aug 10 23:04 mydata_md5.txt
    ## -r--r--r-- 1 guest users 14 Aug 10 23:04 data2.csv

Preventing Modification
-----------------------

``` bash
echo "This is not the data you are looking for" > $DATA2
```

    ## bash: /tmp/Rtmpgogiql/rr_raw_data/data2.csv: Permission denied

``` bash
sed s/7/3/ -i $DATA2
```

    ## sed: couldn't open temporary file /tmp/Rtmpgogiql/rr_raw_data/sedzW9P7T: Permission denied

Preventing Deletion
-------------------

``` bash
rm -rf $RAWDATA_DIR
```

    ## rm: cannot remove '/tmp/Rtmpgogiql/rr_raw_data/mydata_md5.txt': Permission denied
    ## rm: cannot remove '/tmp/Rtmpgogiql/rr_raw_data/data1.csv': Permission denied
    ## rm: cannot remove '/tmp/Rtmpgogiql/rr_raw_data/data2.csv': Permission denied

Archiving Raw Sequence Data @ NCBI
----------------------------------

-   [GEO: Gene "Expression" Data](https://www.ncbi.nlm.nih.gov/geo/)
    -   GEO: Gene Expression Omnibus
        -   RNA-Seq
        -   ChIP-Seq
-   [SRA: Everything else](https://www.ncbi.nlm.nih.gov/sra)
    -   SRA: Sequence Read Archive

"Archive First" Movement
------------------------

-   Free Backup
-   Organize metadata

Alternatives to NCBI
--------------------

-   [International Nucleotide Sequence Database Collaboration](http://www.insdc.org)
    -   [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena)
    -   [DNA Database of Japan](https://www.ddbj.nig.ac.jp/index-e.html)

Other Data Types
----------------

???????????????

Compute Environment
===================

Reproducible Computing Environment
----------------------------------

### Containerization for Reproducible Research

-   *Versioning*: Lock down the specific computing environment used for an analysis
-   Portability: Runs on Linux, Mac, and Windows
-   Sharebility: [Docker Hub](https://hub.docker.com)/[Singularity Hub](https://singularity-hub.org)
-   Scalability: Runs on a laptop, massive server, and everything in between

### Container Platforms

-   [Docker](https://docs.docker.com/get-started/#docker-concepts)
-   [Singularity](http://singularity.lbl.gov/)
-   Etc

Docker Demo
-----------

Analysis Process
================

Analysis Process
----------------

-   Script Everything
-   Share Everything
-   Version Control

Script Everything
-----------------

-   R, Python, Shell Scripts, Rmarkdown, Makefiles, Jupyter, etc
-   [No Excel](https://www.bloomberg.com/news/articles/2013-04-18/faq-reinhart-rogoff-and-the-excel-error-that-changed-history)

Share Everything
----------------

-   Scripts
    -   run parameters are embedded
-   Metadata
-   Documentation
-   Manuscript (optional)

Version Control
---------------

### Heard of Version Control?

### Heard of Github?

<img src="./github_news.png" height="600px" />

What is Version Control?
------------------------

<img src="http://www.phdcomics.com/comics/archive/phd101212s.gif" width="400px" />

<https://stackoverflow.com/a/1408464>

What is Version Control?
------------------------

-   Track
-   Backup
-   Rewind
-   Branch
-   Collaborate
-   Publish

What is Version Control?
------------------------

<img src="./git_diff.png" height="400px" />

Version Control Software
------------------------

-   [Git](https://git-scm.com/)
-   [Mercurial](https://www.mercurial-scm.org/)
-   etc

Git-repository Hosts
--------------------

-   Github ([Education Discount](https://help.github.com/categories/teaching-and-learning-with-github-education/))
-   [Bitbucket](https://bitbucket.org)
-   [Gitlab](https://about.gitlab.com)
    -   [Duke Gitlab](https://gitlab.oit.duke.edu)

Organization
============

My strategy
-----------

1.  Raw Data directory: must be *read-only*
2.  Output directory: everything generated by a script
3.  Git repository: Code and metadata

Alternatives
------------

1.  [ProjectTemplate](https://swcarpentry.github.io/r-novice-gapminder/02-project-intro/#tip-projecttemplate---a-possible-solution)
2.  [Resources on Project Directory Organization](https://discuss.ropensci.org/t/resources-on-project-directory-organization/340)
3.  [A Quick Guide to Organizing Computational Biology Projects](https://doi.org/10.1371/journal.pcbi.1000424)
4.  [Designing projects](https://nicercode.github.io/blog/2013-04-05-projects/)

Resources
=========

Version Control
---------------

-   [Why should I use version control?](https://stackoverflow.com/questions/1408450/why-should-i-use-version-control)

### Git for Version Control

-   Hands On: [Introduction to Version Control with Git](https://swcarpentry.github.io/git-novice/)
-   [Installing Git](https://swcarpentry.github.io/workshop-template/#git)
-   [Sourcetree: a free GUI for git](https://www.sourcetreeapp.com)
-   [Git in RStudio](https://gitlab.oit.duke.edu/IBIEM/IBIEM_2017_2018/blob/master/git_material/git_overview.md)
-   [Happy Git and GitHub for the useR](http://happygitwithr.com)
-   [Excuse me, do you have a moment to talk about version control?](https://doi.org/10.7287/peerj.preprints.3159v2)

Scientific Computing Advice
---------------------------

-   [Good enough practices in scientific computing](https://doi.org/10.1371/journal.pcbi.1005510)
-   [Best Practices for Scientific Computing](https://doi.org/10.1371/journal.pbio.1001745)
