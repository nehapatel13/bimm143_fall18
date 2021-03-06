# Revisit function things
# below is sourcing our old 

source("http://tinyurl.com/rescale-R")

rescale(1:10)
rescale(c(1, 5, "string"))

rescale2(c(1, 5, "string"))


##  write a both_na() function
# Start with a simple setup where we know 
# what the answer should be 

# Lets define an example x and y 
x <- c( 1, 2, NA, 3, NA) 
y<-c(NA,3,NA,3, 4)

is.na(x)
is.na(y)

is.na(x) & is.na(y)

sum(is.na(x))

# this works
# basically counts the number of TRUES
sum(is.na(x) & is.na(y))

# Our first function for NA finding 
both_na <- function(x, y) {
  sum(is.na(x) & is.na(y))
}
both_na(x, y)

# eejit proofing
x <- c(NA, NA, NA)
y1 <- c(1, NA, NA)
y2 <- c(1, NA, NA, NA)

both_na2(x, y2)
res <- both_na3(x, y1)

# find common genes from 2 data sets 
# simplify further to single vectors
x <- df1$IDs
y <- df2$IDs

# this just tells you what genes are in both data sets
# not the actual values 
intersect(x, y)

# gives which values of x are in y as a T/F 
x %in% y

inds <- x %in% y
x[inds]

y[y %in% x]
# only gives which IDs of Y are also in X
y %in% x

# gives you IDs with the values of which ones are in both groups 
df1[df1$IDs %in% df2$IDs, ]
df2[df2$IDs %in% df1$IDs, ]
