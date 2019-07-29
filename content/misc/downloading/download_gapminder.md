This will download the gapminder data into a directory named `r_lessons` in your home directory

``` r
# setup download directory
data.dir = path.expand("~/r_lessons/")
dir.create(data.dir)
```

    ## Warning in dir.create(data.dir): '/Users/josh/r_lessons' already exists

``` r
data.url = "https://raw.githubusercontent.com/resbaz/r-novice-gapminder-files/master/data/gapminder-FiveYearData.csv"

download.file(data.url,
              file.path(data.dir, basename(data.url)))
```
