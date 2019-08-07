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

    ## # A tibble: 738,764 x 4
    ##     rcid country                  country_code vote 
    ##    <int> <chr>                    <chr>        <fct>
    ##  1     3 United States of America US           yes  
    ##  2     3 Canada                   CA           no   
    ##  3     3 Cuba                     CU           yes  
    ##  4     3 Haiti                    HT           yes  
    ##  5     3 Dominican Republic       DO           yes  
    ##  6     3 Mexico                   MX           yes  
    ##  7     3 Guatemala                GT           yes  
    ##  8     3 Honduras                 HN           yes  
    ##  9     3 El Salvador              SV           yes  
    ## 10     3 Nicaragua                NI           yes  
    ## # … with 738,754 more rows

-   The `un_roll_calls` dataset contains information on each roll call vote of the United Nations General Assembly.

``` r
un_roll_calls
```

    ## # A tibble: 5,429 x 9
    ##     rcid session importantvote date       unres  amend  para short  descr  
    ##    <int>   <dbl>         <dbl> <date>     <chr>  <dbl> <dbl> <chr>  <chr>  
    ##  1     3       1             0 1946-01-01 R/1/66     1     0 AMEND… TO ADO…
    ##  2     4       1             0 1946-01-02 R/1/79     0     0 SECUR… TO ADO…
    ##  3     5       1             0 1946-01-04 R/1/98     0     0 VOTIN… "TO AD…
    ##  4     6       1             0 1946-01-04 R/1/1…     0     0 DECLA… TO ADO…
    ##  5     7       1             0 1946-01-02 R/1/2…     1     0 GENER… "TO AD…
    ##  6     8       1             0 1946-01-05 R/1/2…     1     0 ECOSO… TO ADO…
    ##  7     9       1             0 1946-02-05 R/1/3…     0     0 POST-… TO OPE…
    ##  8    10       1             0 1946-02-05 R/1/3…     1     1 U.N. … "TO AD…
    ##  9    11       1             0 1946-02-05 R/1/3…     0     0 TRUST… TO ADO…
    ## 10    12       1             0 1946-02-06 R/1/3…     1     1 COUNC… TO ADO…
    ## # … with 5,419 more rows

-   The `un_roll_call_issues` dataset contains (topic) classifications of roll call votes of the United Nations General Assembly. Many votes had no topic, and some have more than one.

``` r
un_roll_call_issues
```

    ## # A tibble: 5,281 x 3
    ##     rcid short_name issue               
    ##    <int> <chr>      <chr>               
    ##  1  3372 me         Palestinian conflict
    ##  2  3658 me         Palestinian conflict
    ##  3  3692 me         Palestinian conflict
    ##  4  2901 me         Palestinian conflict
    ##  5  3020 me         Palestinian conflict
    ##  6  3217 me         Palestinian conflict
    ##  7  3298 me         Palestinian conflict
    ##  8  3429 me         Palestinian conflict
    ##  9  3558 me         Palestinian conflict
    ## 10  3625 me         Palestinian conflict
    ## # … with 5,271 more rows

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

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-e0603c3529247b91d4b8">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200"],["Afghanistan","Albania","Algeria","Andorra","Angola","Antigua and Barbuda","Argentina","Armenia","Australia","Austria","Azerbaijan","Bahamas","Bahrain","Bangladesh","Barbados","Belarus","Belgium","Belize","Benin","Bhutan","Bolivia (Plurinational State of)","Bosnia and Herzegovina","Botswana","Brazil","Brunei Darussalam","Bulgaria","Burkina Faso","Burundi","Cabo Verde","Cambodia","Cameroon","Canada","Central African Republic","Chad","Chile","China","Colombia","Comoros","Congo","Costa Rica","Côte D'Ivoire","Croatia","Cuba","Cyprus","Czech Republic","Czechoslovakia","Democratic People's Republic of Korea","Democratic Republic of the Congo","Denmark","Djibouti","Dominica","Dominican Republic","Ecuador","Egypt","El Salvador","Equatorial Guinea","Eritrea","Estonia","Ethiopia","Federal Republic of Germany","Fiji","Finland","France","Gabon","Gambia (Islamic Republic of the)","Georgia","German Democratic Republic","Germany","Ghana","Greece","Grenada","Guatemala","Guinea","Guinea Bissau","Guyana","Haiti","Honduras","Hungary","Iceland","India","Indonesia","Iran (Islamic Republic of)","Iraq","Ireland","Israel","Italy","Jamaica","Japan","Jordan","Kazakhstan","Kenya","Kiribati","Kuwait","Kyrgyzstan","Lao People's Democratic Republic","Latvia","Lebanon","Lesotho","Liberia","Libya","Liechtenstein","Lithuania","Luxembourg","Madagascar","Malawi","Malaysia","Maldives","Mali","Malta","Marshall Islands","Mauritania","Mauritius","Mexico","Micronesia (Federated States of)","Monaco","Mongolia","Montenegro","Morocco","Mozambique","Myanmar","Namibia","Nauru","Nepal","Netherlands","New Zealand","Nicaragua","Niger","Nigeria","Norway","Oman","Pakistan","Palau","Panama","Papua New Guinea","Paraguay","Peru","Philippines","Poland","Portugal","Qatar","Republic of Korea","Republic of Moldova","Romania","Russian Federation","Rwanda","Saint Kitts and Nevis","Saint Lucia","Saint Vincent and the Grenadines","Samoa","San Marino","Sao Tome and Principe","Saudi Arabia","Senegal","Seychelles","Sierra Leone","Singapore","Slovakia","Slovenia","Solomon Islands","Somalia","South Africa","South Sudan","Spain","Sri Lanka","Sudan","Suriname","Swaziland","Sweden","Switzerland","Syrian Arab Republic","Taiwan, Province of China","Tajikistan","Thailand","The former Yugoslav Republic of Macedonia","Timor-Leste","Togo","Tonga","Trinidad and Tobago","Tunisia","Turkey","Turkmenistan","Tuvalu","Uganda","Ukraine","United Arab Emirates","United Kingdom of Great Britain and Northern Ireland","United Republic of Tanzania","United States of America","Uruguay","Uzbekistan","Vanuatu","Venezuela, Bolivarian Republic of","Viet Nam","Yemen","Yemen Arab Republic","Yemen People's Republic","Yugoslavia","Zambia","Zanzibar","Zimbabwe"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>country<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"order":[],"autoWidth":false,"orderClasses":false,"columnDefs":[{"orderable":false,"targets":0}]}},"evals":[],"jsHooks":[]}</script>
<!--/html_preserve-->
