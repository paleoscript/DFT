## supplementary material for Virág & Karádi (in prep.)       ##
## (submitted to Palaeontology during September 2022)         ##

## written by Attila VIRAG (1, 2)                             ##
## 1: Dept. of Mineralogy and Geology, University of Debrecen ##
## 2: MTA-MTM-ELTE Research Group for Paleontology, Budapest  ##


#---------------------------------------------------------#
## FUNCTIONS ----

#---------------------------------------------------------#
## complex numbers ----

# class Complex
setClass("Complex", slots = list(re = "numeric", im = "numeric"))

# algebraic operations on complex numbers

# addition
addC <- function(c1, c2) {
  re <- c1@re + c2@re
  im <- c1@im + c2@im
  return(new("Complex", re = re, im = im))
}

# subtraction
subC <- function(c1, c2) {
  re <- c1@re - c2@re
  im <- c1@im - c2@im
  return(new("Complex", re = re, im = im))
}

# multiplication
mulC <- function(c1, c2) {
  re <- c1@re * c2@re - c1@im * c2@im
  im <- c1@re * c2@im + c1@im * c2@re
  return(new("Complex", re = re, im = im))
}

# division
divC <- function(c1, c2) {
  re <- (c1@re * c2@re + c1@im * c2@im) / (c1@im^2 + c2@im^2)
  im <- (c1@im * c2@re - c1@re * c2@im) / (c1@im^2 + c2@im^2)
}


#---------------------------------------------------------#
## shape functions ----

# shape is a 2-column matrix containing a set of coordinates
centerShape <- function(shape, resize = TRUE, size = 5000) {
  
  resShape <- matrix(NA, nrow = nrow(shape), ncol = 2)
  
  resShape[, 1] <- shape[, 1] - mean(shape[, 1])
  resShape[, 2] <- shape[, 2] - mean(shape[, 2])
  
  if (resize == TRUE) {
    cSize         <- sum(sqrt(resShape[, 1]^2 + resShape[, 2]^2))
    resShape[, 1] <- resShape[, 1] * (size/cSize)
    resShape[, 2] <- resShape[, 2] * (size/cSize)
  }
  
  return(resShape)
}


#---------------------------------------------------------#
## complex dft ----

# SMITH 2002, CHAPTER 31, EQUATION 31-6
# DOI: https://doi.org/10.1016/B978-0-7506-7444-7.50073-X

# forward complex discrete Fourier transform
# x is a 2-column matrix containing a set of coordinates
dft <- function(x) {
  
  N <- nrow(x)
  
  X <- matrix(NA, ncol = 3, nrow = N)
  colnames(X) <- c("re", "im", "k")
  
  # array indices start at 1 in R
  for (k in 1:N) {
    Xk <- new("Complex", re = 0, im = 0)
    
    for (n in 1:N) {
      xn <- new("Complex",
                re = x[n, 1],
                im = x[n, 2])
      
      expr <- 2*pi*(k-1)*n/N
      comp <- new("Complex",
                  re =  cos(expr), 
                  im = -sin(expr)) 
      
      Xk <- addC(Xk, mulC(xn, comp))
    }
    
    # scaling factor of 1/N is introduced 
    # to make the reconstructed signal be 
    # identical to the original signal
    X[k, 1] <- (1/N) * Xk@re
    X[k, 2] <- (1/N) * Xk@im
    X[k, 3] <- k - 1
  }
  
  return(X)
}


# SMITH 2002, CHAPTER 31, EQUATION 31-8
# simplified in accordance with the detailed explanation therein
# DOI: https://doi.org/10.1016/B978-0-7506-7444-7.50073-X

# inverse complex discrete Fourier transform
# X is a 3-column matrix containing the result of a dft
invDft <- function(X, epi = Inf) {
  
  N <- nrow(X)
  epi <- min(epi, N)
  
  x <- matrix(NA, ncol = 2, nrow = N)
  colnames(x) <- c("x", "y")
  
  # array indices start at 1 in R
  for (n in 1:N) {
    xn <- new("Complex", re = 0, im = 0)
    
    for (k in 1:epi) {
      Xk <- new("Complex", 
                re = X[k, 1], 
                im = X[k, 2])
      
      expr  <- 2*pi*(X[k, 3])*n/N
      comp <- new("Complex", 
                  re = cos(expr), 
                  im = -sin(expr))
      
      xn <- addC(xn, mulC(Xk, comp))
    }
    
    x[n, 1] <- xn@re
    x[n, 2] <- xn@im
  }
  return(x)
}

# X is a 3-column matrix containing the result of a dft
normDft <- function(X) {
  
  N <- nrow(X)
  
  # reordering the result of a dft by
  # pairing up the respective descriptors 
  # from the opposite ends of the frequency spectrum
  order <- numeric()
  for(i in 2:(N/2)) {
    order <- c(order, i, N + 2 - i)
  }
  order <- c(order, 1, N/2)
  
  X <- X[order, ]
  
  # rescaling the major axis of the first ellipse to 1
  # and the smaller ellipses accordingly
  unit <- sqrt(X[1, 1]^2 + X[1, 2]^2) +
    sqrt(X[2, 1]^2 + X[2, 2]^2)
  X[, 1] <- X[, 1] / unit
  X[, 2] <- X[, 2] / unit
  
  # using the first ellipse as a reference
  ref <- invDft(X, epi = 2) 
  
  # searching the first ellipse 
  # for its farthest point from the centre
  dist <- ref[, 1]^2 + ref[, 2]^2
  # narrowing the focus of the search
  # to the lower half of the coordinate system
  f <- ref[, 2] < 0
  i <- which(dist[f] == max(dist[f]))
  i <- which(f)[i]
  
  # rotational correction
  c1 <- atan2(ref[i, 2], ref[i, 1])
  # starting phase correction
  c2  <- 2*pi*(X[, 3])*(i-1)/N
  
  amp   <- sqrt(X[, 1]^2 + X[, 2]^2)
  phase <- atan2(X[, 2], X[, 1])
  phase <- phase - c1 - c2
  X[, 1] <- amp * cos(phase)
  X[, 2] <- amp * sin(phase)
  
  return(X)
}


#---------------------------------------------------------#
## EXAMPLE ----

#---------------------------------------------------------#
## working directory ----

# setting working directory (in RStudio)
setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path)))


#---------------------------------------------------------#
## packages ----

# equidistantCurve()
lib <- suppressWarnings(require("Morpho"))
if(!lib) {install.packages("Morpho"); library("Morpho")}

# load.image()
lib <- suppressWarnings(require("imager"))
if(!lib) {install.packages("imager"); library("imager")}


#---------------------------------------------------------#
## shape descriptors ----

# loading the raw image of the..
# A. uniformis depicted on Fig. 1
img <- load.image("uniformis.jpg")

# loading raw shape descriptors from .csv files

# element outline
elem <- read.csv("uniformis_elem.csv", 
                 header = TRUE, 
                 sep = ",", dec = ".")[, 6:7]

# keel outline
keel <- read.csv("uniformis_keel.csv", 
                 header = TRUE, 
                 sep = ",", dec = ".")[, 6:7]

# plotting the specimen using its shape descriptors
# SLOW!
{
  plot(img, axes = F)
  points(elem[-(1:2), ], pch = "×", cex = 0.8, 
         col = "yellow")
  points(keel[-(1:2), ], pch = "×", cex = 0.8, 
         col = "yellow")
  points(keel[1:2, ], pch = "+",
         col = "yellow")
}


#---------------------------------------------------------#
## equidistant points ----

# equidistant resampling of the shape descriptors
points <- 100   # a pre-defined number of points
shape  <- elem  # alternatively: sheep <- keel
{
  # pit and cusp
  LM1   <- shape[1, ]
  LM2   <- shape[2, ]
  
  # contour
  shape <- shape[-(1:2), ]
  shape <- equidistantCurve(as.matrix(shape), 
                            n = points + 1, 
                            open = FALSE)
  shape <- as.data.frame(shape[-(points + 1), ])
  
  colnames(shape) <- c("X", "Y")
  shape <- rbind(LM1, LM2, shape)
  
  # centering and scaling of the contour
  shape <- centerShape(shape)
}


#---------------------------------------------------------#
## dft() ----

# forward complex discrete Fourier transformation
dftRes  <- dft(shape[-(1:2), ])

# normalization
normRes <- normDft(dftRes)

# first ellipse
# (retreived by an inverse complex DFT..
# based on the first 2 epicycles)
firstE  <- invDft(normRes, epi = 2)

# smoothed contour
# (retreived by an inverse complex DFT..
# based on the first 20 epicycles)
invRes  <- invDft(normRes, epi = 20)

# plotting the selected contour.. 
# using its Fourier descriptors
{
  plot(invRes, asp = 1, type = "n")
  polygon(invRes, lwd = 2, col = "grey80")
  points(invRes[1, 1], invRes[1, 2], 
         pch = 16, col = "red")
  polygon(firstE, lwd = 2, lty = 2)
  abline(v = 0, lty = 3)
  abline(h = 0, lty = 3)
}

