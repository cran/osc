\name{osc.buffer}
\alias{osc.buffer}
\title{
Simple Buffer algorithm
}
\description{
Simple buffer based on euclidean distance are created around all cells equal to one.
}
\usage{
osc.buffer(input, width) 
}
\arguments{
  \item{input}{
    Matrix or Raster containing 1 indicating a cluster, no NA values are allowed
}
  \item{width}{
    Width of the buffer in cells
}
}

\value{
 Returns matrix or raster, depending on input, with 1 for cluster and -1 for surrounding buffer
}

\examples{
data(landcover)
landcover[landcover[]>1] <- 0
plot(osc.buffer(landcover, 4))
}
\keyword{utils}
