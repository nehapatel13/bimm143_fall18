class7rmd
================
Neha Patel
November 1, 2018

``` r
source("http://tinyurl.com/rescale-R")
```

    ## Warning in file(filename, "r", encoding = encoding): "internal" method
    ## cannot handle https redirection to: 'https://bioboot.github.io/bggn213_f17/
    ## class-material/rescale.R'

    ## Warning in file(filename, "r", encoding = encoding): "internal" method
    ## failed, so trying "libcurl"

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
#rescale(c(1, 5, "string"))

#rescale2(c(1, 5, "string"))


##  write a both_na() function
# Start with a simple setup where we know 
# what the answer should be 

# Lets define an example x and y 
x <- c( 1, 2, NA, 3, NA) 
y<-c(NA,3,NA,3, 4)

is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum(is.na(x))
```

    ## [1] 2

``` r
# this works
# basically counts the number of TRUES
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
# Our first function for NA finding 
both_na <- function(x, y) {
  sum(is.na(x) & is.na(y))
}
both_na(x, y)
```

    ## [1] 1

``` r
# eejit proofing
x <- c(NA, NA, NA)
y1 <- c(1, NA, NA)
y2 <- c(1, NA, NA, NA)

#both_na2(x, y2)
#res <- both_na3(x, y1)

# find common genes from 2 data sets 
# simplify further to single vectors
x <- df1$IDs
y <- df2$IDs

# this just tells you what genes are in both data sets
# not the actual values 
intersect(x, y)
```

    ## [1] "gene2" "gene3"

``` r
# gives which values of x are in y as a T/F 
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
inds <- x %in% y
x[inds]
```

    ## [1] "gene2" "gene3"

``` r
y[y %in% x]
```

    ## [1] "gene2" "gene3"

``` r
# only gives which IDs of Y are also in X
y %in% x
```

    ## [1]  TRUE FALSE  TRUE FALSE

``` r
# gives you IDs with the values of which ones are in both groups 
df1[df1$IDs %in% df2$IDs, ]
```

    ##     IDs exp
    ## 2 gene2   1
    ## 3 gene3   1

``` r
df2[df2$IDs %in% df1$IDs, ]
```

    ##     IDs exp
    ## 1 gene2  -2
    ## 3 gene3   1
