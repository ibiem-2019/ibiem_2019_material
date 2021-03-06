Let's take a look at the voting history of countries in the United Nations General Assembly. We will be using data from the **unvotes** package. Additionally, we will make use of the **tidyverse** and **lubridate** packages for the analysis, and the **DT** package for interactive display of tabular output.

``` r
library(tidyverse)
library(unvotes)
library(lubridate)
library(DT)
```

Data
----

The unvotes package provides three datasets we can work with: `un_roll_calls`, `un_roll_call_issues`, and `un_votes`. Each of these datasets contains a variable called `rcid`, the roll call id, which can be used as a unique identifier to join them with each other.

-   The `un_votes` dataset provides information on the voting history of the United Nations General Assembly. It contains one row for each country-vote pair.

``` r
un_votes
```

-   The `un_roll_calls` dataset contains information on each roll call vote of the United Nations General Assembly.

``` r
un_roll_calls
```

-   The `un_roll_call_issues` dataset contains (topic) classifications of roll call votes of the United Nations General Assembly. Many votes had no topic, and some have more than one.

``` r
un_roll_call_issues
```

Analysis
--------

First, let's take a look at how often each country voted "Yes" on a resolution in each year. We'll visualize the results, so let's pick a few countries of interest first,

``` r
country_list <- c("Ghana",
                  "Japan",
                  "Turkey",
                  "United States of America")
```

and focus our analysis on them.

``` r
un_votes %>%
  filter(country %in% country_list) %>%
  inner_join(un_roll_calls, by = "rcid") %>%
  group_by(year = year(date), country) %>%
  summarize(
    votes = n(),
    percent_yes = mean(vote == "yes")
    ) %>%
  ggplot(mapping = aes(x = year, y = percent_yes, color = country)) +
    geom_line() +
    labs(
      title = "Percentage of 'Yes' votes in the UN General Assembly",
      subtitle = "1946 to 2015",
      y = "% Yes",
      x = "Year",
      color = "Country"
    )
```

![](020_unvotes_files/figure-markdown_github/plot-yearly-yes-1.png)

Next, let's create a visualization that displays how the voting record of the United States changed over time on a variety of issues, and compares it to another country. The other country we'll display is Turkey.

``` r
un_votes %>%
  filter(country %in% c("United States of America", "Turkey")) %>%
  inner_join(un_roll_calls, by = "rcid") %>%
  inner_join(un_roll_call_issues, by = "rcid") %>%
  group_by(country, year = year(date), issue) %>%
  summarize(
    votes = n(),
    percent_yes = mean(vote == "yes")
    ) %>%
  filter(votes > 5) %>%  # only use records where there are more than 5 votes
  ggplot(mapping = aes(x = year, y = percent_yes, color = country)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    facet_wrap(~ issue) +
    labs(
      title = "Percentage of 'Yes' votes in the UN General Assembly",
      subtitle = "1946 to 2015",
      y = "% Yes",
      x = "Year",
      color = "Country"
    )
```

![](020_unvotes_files/figure-markdown_github/plot-yearly-yes-issue-1.png)

Knitting
--------

1.  Knit to HTML
2.  Knit to PDF
3.  Knit to Word

Challenge: Your Turn
--------------------

You can easily change which countries are being plotted by changing which countries the code above `filter`s for.

Your job is to change one or more of the countries plotted in each of the plots above. Note that the country name should be spelled and capitalized exactly the same way as it appears in the data. See the [Appendix](#appendix) for a list of the countries in the data.

References
----------

1.  David Robinson (2017). unvotes: United Nations General Assembly Voting Data. R package version 0.2.0. <https://CRAN.R-project.org/package=unvotes>.
2.  Erik Voeten "Data and Analyses of Voting in the UN General Assembly" Routledge Handbook of International Organization, edited by Bob Reinalda (published May 27, 2013).
3.  Much of the analysis has been modeled on the examples presented in the [unvotes package vignette](https://cran.r-project.org/web/packages/unvotes/vignettes/unvotes.html).
4.  This notebook was shamelessly borrowed from [Mine Çetinkaya-Rundel](https://github.com/mine-cetinkaya-rundel/teach-data-sci-icots2018/blob/master/01-01-start/unvotes.Rmd)

Appendix
--------

Below is a list of countries in the dataset:

``` r
un_votes %>% 
  select(country) %>%
  arrange(country) %>% 
  distinct() %>%
  datatable()
```

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
