\name{shelltox}
\alias{shelltox}
\docType{data}
\title{Shell Toxicology data}
\description{
This is a classical developmental toxicology data set. Pregnant banded Dutch rabbits
were treated with one of four levels of a chemical. The actual doses are not known,
instead the groups are designated as Control, Low, Medium, and High. Before term the
animals were sacrificed, and the total number of fetuses, as well as the number affected
by the treatment was recorded.
}

\usage{data(shelltox)}
\format{
  A 'CBData' object, that is a data frame  with the following variables
  \tabular{ll}{
    Trt          \tab factor giving treatment group\cr
    ClusterSize  \tab the size of the litter\cr
    NResp        \tab the number of affected fetuses\cr
    Freq         \tab the number of litters with the given ClusterSize/NResp combination\cr
  }
}
\source{
Paul, S. R. (1982) Analysis of proportions of affected foetuses in teratological experiments. \emph{Biometrics}, 38, 361-370.  

This data set has been analyzed (and listed) in numerous papers, including

Rao, J. N. K. & Scott, A. J. (1992) A Simple Method for the Analysis of Clustered Data. \emph{Biometrics}, 48, 577-586.

George, E. O. & Kodell, R. L. (1996) Tests of Independence, Treatment Heterogeneity, and Dose-Related Trend With Exchangeable Binary Data.
\emph{Journal of the American Statistical Association}, 91, 1602-1610.

Lee, S. (2003) Analysis of the Binary Littermate Data in the One-Way Layout. \emph{Biometrical Journal}, 45, 195-206.
}
\examples{
  data(shelltox)
  stripchart(I(NResp/ClusterSize)~Trt, cex=sqrt(shelltox$Freq), data=shelltox, pch=1, 
             method="jitter", vertical=TRUE, ylab="Proportion affected")
}
\keyword{datasets}
