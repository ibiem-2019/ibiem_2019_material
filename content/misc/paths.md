Unix
====

List Files and Directories
--------------------------

``` bash
ls
```

    ## bash_test_dir
    ## downloading
    ## downloading.md
    ## downloading.Rmd
    ## github_notes.md
    ## github_notes.Rmd
    ## ibiem_computing.md
    ## ibiem_computing.Rmd
    ## paths.md
    ## paths.nb.html
    ## paths.Rmd
    ## r_test_dir

``` bash
ls -a
```

    ## .
    ## ..
    ## bash_test_dir
    ## downloading
    ## downloading.md
    ## downloading.Rmd
    ## github_notes.md
    ## github_notes.Rmd
    ## .gitignore
    ## ibiem_computing.md
    ## ibiem_computing.Rmd
    ## paths.md
    ## paths.nb.html
    ## paths.Rmd
    ## r_test_dir

``` bash
ls -F
```

    ## bash_test_dir/
    ## downloading/
    ## downloading.md
    ## downloading.Rmd
    ## github_notes.md
    ## github_notes.Rmd
    ## ibiem_computing.md
    ## ibiem_computing.Rmd
    ## paths.md
    ## paths.nb.html
    ## paths.Rmd
    ## r_test_dir/

Working Directory
-----------------

``` bash
pwd
```

    ## /home/guest/IBIEM_2018_2019/content/misc

``` bash
cd ..
pwd
```

    ## /home/guest/IBIEM_2018_2019/content

``` bash
pwd
```

    ## /home/guest/IBIEM_2018_2019/content/misc

Absolute vs Relative Paths
--------------------------

``` bash
head paths.Rmd
```

    ## ---
    ## title: "Paths in R, RStudio, and Unix"
    ## output:
    ##   md_document:
    ##     variant: markdown_github
    ## ---
    ## 
    ## # Unix
    ## ## List Files and Directories
    ## ```{bash}

``` bash
head /home/guest/IBIEM_2018_2019/content/misc/paths.Rmd
```

    ## ---
    ## title: "Paths in R, RStudio, and Unix"
    ## output:
    ##   md_document:
    ##     variant: markdown_github
    ## ---
    ## 
    ## # Unix
    ## ## List Files and Directories
    ## ```{bash}

Environment Variables
---------------------

``` bash
echo "Hello world"
```

    ## Hello world

``` bash
echo $HOME
```

    ## /home/guest

``` bash
echo $PWD
```

    ## /home/guest/IBIEM_2018_2019/content/misc

``` bash
ls $HOME/IBIEM_2018_2019
```

    ## content
    ## IBIEM_2018_2019.Rproj
    ## misc
    ## README.md
    ## README.Rmd

Tilde Expansion
---------------

``` bash
echo ~
```

    ## /home/guest

``` bash
ls ~
```

    ## challenge_2
    ## IBIEM_2018_2019
    ## R

``` bash
cd ~
pwd
```

    ## /home/guest

Make a New Directory
--------------------

``` bash
mkdir bash_test_dir
ls -F
```

    ## mkdir: cannot create directory ‘bash_test_dir’: File exists
    ## bash_test_dir/
    ## downloading/
    ## downloading.md
    ## downloading.Rmd
    ## github_notes.md
    ## github_notes.Rmd
    ## ibiem_computing.md
    ## ibiem_computing.Rmd
    ## paths.md
    ## paths.nb.html
    ## paths.Rmd
    ## r_test_dir/

R
=

List Files and Directories
--------------------------

``` r
list.files()
```

    ##  [1] "bash_test_dir"       "downloading"         "downloading.md"     
    ##  [4] "downloading.Rmd"     "github_notes.md"     "github_notes.Rmd"   
    ##  [7] "ibiem_computing.md"  "ibiem_computing.Rmd" "paths.md"           
    ## [10] "paths.nb.html"       "paths.Rmd"           "r_test_dir"

``` r
list.dirs()
```

    ## [1] "."               "./bash_test_dir" "./downloading"   "./r_test_dir"

Working Directory
-----------------

``` r
getwd()
```

    ## [1] "/home/guest/IBIEM_2018_2019/content/misc"

Environment Variables
---------------------

``` r
Sys.getenv()
```

    ## CLICOLOR_FORCE        1
    ## DISPLAY               :0
    ## EDITOR                vi
    ## GIT_ASKPASS           rpostback-askpass
    ## HOME                  /home/guest
    ## LANG                  en_US.UTF-8
    ## LD_LIBRARY_PATH       /usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/default-java/jre/lib/amd64/server:/usr/lib/R/lib::/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/default-java/jre/lib/amd64/server
    ## LN_S                  ln -s
    ## LOGNAME               guest
    ## MAKE                  make
    ## NOT_CRAN              true
    ## PAGER                 /usr/bin/pager
    ## PATH                  /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/conda/bin:/usr/lib/ChimeraSlayer
    ## PWD                   /home/guest/IBIEM_2018_2019/content/misc
    ## R_ARCH                
    ## R_BROWSER             xdg-open
    ## R_BZIPCMD             /bin/bzip2
    ## R_DOC_DIR             /usr/share/R/doc
    ## R_GZIPCMD             /bin/gzip -n
    ## R_HOME                /usr/lib/R
    ## R_INCLUDE_DIR         /usr/share/R/include
    ## R_LIBS                /home/guest/R/x86_64-pc-linux-gnu-library/3.4:/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library
    ## R_LIBS_SITE           /usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library
    ## R_LIBS_USER           ~/R/x86_64-pc-linux-gnu-library/3.4
    ## RMARKDOWN_MATHJAX_PATH
    ##                       /usr/lib/rstudio-server/resources/mathjax-26
    ## RMARKDOWN_PREVIEW_DIR
    ##                       /tmp/RtmpG4GbYV
    ## R_PAPERSIZE           letter
    ## R_PAPERSIZE_USER      letter
    ## R_PDFVIEWER           /usr/bin/xdg-open
    ## R_PLATFORM            x86_64-pc-linux-gnu
    ## R_PRINTCMD            /usr/bin/lpr
    ## R_RD4PDF              times,inconsolata,hyper
    ## R_SESSION_TMPDIR      /tmp/RtmpbqR7zF
    ## R_SHARE_DIR           /usr/share/R/share
    ## RS_RPOSTBACK_PATH     /usr/lib/rstudio-server/bin/rpostback
    ## RSTUDIO               1
    ## RSTUDIO_CONSOLE_COLOR
    ##                       256
    ## RSTUDIO_CONSOLE_WIDTH
    ##                       116
    ## RSTUDIO_HTTP_REFERER
    ##                       https://ipyn-az-09.oit.duke.edu:50000/
    ## RSTUDIO_PANDOC        /usr/lib/rstudio-server/bin/pandoc
    ## RSTUDIO_SESSION_STREAM
    ##                       guest-d
    ## RSTUDIO_USER_IDENTITY
    ##                       guest
    ## RSTUDIO_WINUTILS      bin/winutils
    ## R_SYSTEM_ABI          linux,gcc,gxx,gfortran,?
    ## R_TEXI2DVICMD         /usr/bin/texi2dvi
    ## R_UNZIPCMD            /usr/bin/unzip
    ## R_ZIPCMD              /usr/bin/zip
    ## SED                   /bin/sed
    ## SHELL                 /bin/bash
    ## SHLVL                 0
    ## SSH_AGENT_PID         1000
    ## SSH_ASKPASS           rpostback-askpass
    ## SSH_AUTH_SOCK         /tmp/ssh-0xSSpvIipdRR/agent.999
    ## TAR                   /bin/tar
    ## TERM                  xterm-256color
    ## USER                  guest

``` r
Sys.getenv("HOME")
```

    ## [1] "/home/guest"

``` r
Sys.getenv("PWD")
```

    ## [1] "/home/guest/IBIEM_2018_2019/content/misc"

Path Expansion and Normalization
--------------------------------

``` r
path.expand("~////")
```

    ## [1] "/home/guest////"

``` r
normalizePath("~////")
```

    ## [1] "/home/guest"

``` r
normalizePath("~/not_a_directory")
```

    ## Warning in normalizePath("~/not_a_directory"): path[1]="/home/guest/
    ## not_a_directory": No such file or directory

    ## [1] "/home/guest/not_a_directory"

``` r
path.expand("~/not_a_directory")
```

    ## [1] "/home/guest/not_a_directory"

``` r
path.expand("./paths.Rmd")
```

    ## [1] "./paths.Rmd"

``` r
normalizePath("./paths.Rmd")
```

    ## [1] "/home/guest/IBIEM_2018_2019/content/misc/paths.Rmd"

Path Maniputation
-----------------

``` r
file.path("some_dir", "some_subdir")
```

    ## [1] "some_dir/some_subdir"

``` r
file.path(Sys.getenv("HOME"), "IBIEM_2018_2019")
```

    ## [1] "/home/guest/IBIEM_2018_2019"

``` r
this.file = "/home/guest/IBIEM_2018_2019/content/misc/paths.Rmd"
basename(this.file)
```

    ## [1] "paths.Rmd"

``` r
dirname(this.file)
```

    ## [1] "/home/guest/IBIEM_2018_2019/content/misc"

Make a New Directory
--------------------

``` r
dir.create("r_test_dir")
```

    ## Warning in dir.create("r_test_dir"): 'r_test_dir' already exists

``` r
list.dirs()
```

    ## [1] "."               "./bash_test_dir" "./downloading"   "./r_test_dir"

RStudio
=======

-   Home
-   Project Directory (R cube)
-   ".."
-   "..."
-   Git Pane

References
==========

-   [Navigating Files and Directories](http://swcarpentry.github.io/shell-novice/02-filedir/index.html)
-   [Working With Files and Directories](http://swcarpentry.github.io/shell-novice/03-create/index.html)
