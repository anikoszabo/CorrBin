\name{egde}
\alias{egde}
\docType{data}
\title{EGDE data}
\description{
The data set is based on a developmental toxicity experiment on the effect of ethylene glycol 
diethyl ether (EGDE) on fetal development of New Zealand white rabbits.
In the study, four groups of pregnant does
were randomly assigned to dose levels $0, 25, 50$, and $100$
milligrams per kilogram body weight of EGDE. For each litter
and at each dose level, the adverse response used is the combined
number of fetal malformation and fetal death.
}

\usage{data(egde)}
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
Krewski, D., Zhu, Y., and Fung, K. (1995). Statistical analysis of overdispersed multinomial data from 
developmental toxicity studies. In \emph{Statistics in Toxicology}, Ed. B. Morgan, pp.\ 151--179.  New York: Oxford University Press.
}
\examples{
  data(egde)
  stripchart(I(NResp/ClusterSize)~Trt, cex=sqrt(egde$Freq), data=egde, pch=1, 
             method="jitter", vertical=TRUE, ylab="Proportion affected")
}
\keyword{datasets}
