# Baby weight data input 
weight <- read.table("bimm143_05_rstats/weight_chart.txt" , header = TRUE)

# Make a custom plot 
plot(weight , pch=15 , cex=1.5 , lwd=2 , ylim=c(2,10) , xlab="Age (months)" , ylab="Weight (kg)" , main= "Weight of Elephants" , type="o")

# Section 2
# 2B barplot
counts <- read.table("bimm143_05_rstats/feature_counts.txt", sep = "\t", header=TRUE)
barplot(counts$Count)

# bar plot with correct labels 
par(mar=c(3.1, 11.1, 4.1, 2) , mgp=c(10, 1, 0))
barplot(counts$Count , names.arg = counts$Feature , horiz = TRUE , ylab = "Features" , main = "Number of features in the mouse GRCm38 genome" , las=1 , xlim = c(0,80000))

# histogram plot 

par(mgp=c(3, 1, 0))
x <- c(rnorm(10000) , rnorm(10000) + 4)
hist(x, breaks = 80)

# Section 3
#Color Vectors
male_female <- read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE , sep = "\t")
barplot(male_female$Count, col = c(rainbow(10)))

#Only 2 colors
barplot(male_female$Count, names.arg = male_female$Sample, las=2, col=c("blue2", "red2"))


#Using NROW
barplot(male_female$Count, names.arg = male_female$Sample, las=2, col=rainbow(nrow(male_female)))

        
 genes <- read.table("bimm143_05_rstats/up_down_expression.txt", header = T, sep = "\t")       
nrow(genes)

# Dynamic use of color
meth <- read.table("bimm143_05_rstats/expression_methylation.txt", header = TRUE, sep = "\t")
nrow(meth)

# focus on the data we want to examine
inds <- meth$expression > 0 
mycols2 <- densCols(meth$gene.meth[inds], meth$expression[inds])

plot(meth$gene.meth[inds], meth$expression[inds], col= mycols2)


