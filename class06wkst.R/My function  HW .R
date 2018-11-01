# OUR FUNCTION

# Kinase with drug 
library(bio3d)
x <- "4AKE"
pdb <- read.pdb(x)
chainA <- trim.pdb(pdb, chain = "A", elety = "CA")
sb1 <- chainA$atom$b
plotb3(sb1, sse=chainA, typ="l", ylab="Bfactor")

# Kinase without Drug
library(bio3d)
x <- "1AKE"
pdb <- read.pdb(x)
chainA <- trim.pdb(pdb, chain = "A", elety = "CA")
sb2 <- chainA$atom$b
plotb3(sb2, sse=chainA, typ="l", ylab="Bfactor")

# Kinase with Drug
library(bio3d)
x <- "1E4Y"
pdb <- read.pdb(x)
chainA <- trim.pdb(pdb, chain = "A", elety = "CA")
sb3 <- chainA$atom$b
plotb3(sb3, sse=chainA, typ="l", ylab="Bfactor")

# Plot all together 

# FUNCTION FOR ASSIGNMENT
library(bio3d)
my_function <- function(x) {
  pdb <- read.pdb(x)
  chainA <- trim.pdb(pdb, chain = "A", elety = "CA")
  b <- chainA$atom$b
  plotb3(b, sse=chainA, typ="l", ylab="Bfactor")
}

