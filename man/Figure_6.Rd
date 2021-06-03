\name{Figure_6}
\docType{data}
\alias{Figure_6}
    
\title{Figure_6}

\description{
  
  Speed (time in seconds) and memory (peak RAM in MB) performance of crqa() 
  and piecewiseRQA() for simulated data of increasing size 
  (from 3000 to 7000 data points), compared in  blocks of different 
  sizes (from 1000 to 6500 in increments of 500)

}

\usage{Figure_6}
\format{A data frame with 28 rows and 2 variables:
 \describe{
   \item{speed}{The time in seconds to perform a crqa() analysis}
   \item{memory}{The peak MB occupied in the RAM to perform a crqa() analysis}
   \item{typeRQA}{Whether the analysis contained all data points (full) or was
   computed piecewise (piece)}
  \item{datapoint}{The number of datapoints}
  \item{blocksize}{The dimension of the block}

  }
 }

\keyword{datasets}