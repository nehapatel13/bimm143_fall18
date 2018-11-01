map.colors3 <- function (x, low.high = range(x), palette= cm.colors(100)) {
  ## Description: map the values of the input vector 'x'
  ## to the input colors vector 'palette'
  
  # Determine percent values of the 'high.low' range
  percent <- ((x-low.high[1])/(low.high[2]-low.high[1]))
  # Find corresponding index position in the color 'palette'
  # note , catch for - percent values to 1 
  index <- round (length(palette-1)*percent) + 1
  return (palette[index])
}

# New Function
add <- function(x, y=1) {
  x + y
}

# NEW FUNCTION 
# changing na.rm to TRUE allows you to neglect the NA and use the function normally
rescale <- function(x, na.rm = TRUE) {
    rng <- range(x, na.rm = na.rm)
    answer <- (x - rng[1]) / (rng[2] - rng[1])
    plot(answer, type = "o")
}
rescale(1:20)

# lecture function
rescale <- function(x, na.rm=TRUE, plot=FALSE) {
  if(na.rm) {
    rng <-range(x, na.rm=TRUE)
  } else {
    rng <-range(x)
  }
  print("Hello")
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  print("is it me you are looking for?")
  if(plot) {
    plot(answer, typ="b", lwd=4)
  }
  print("I can see it in ...")
  return(answer)
}
rescale(1:10, plot = TRUE)
