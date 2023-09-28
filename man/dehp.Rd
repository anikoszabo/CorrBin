\name{dehp}
\docType{data}
\alias{dehp}
\title{Developmental toxicology study of DEHP in mice}
\description{
This data set is based on a National Toxicology Program study on diethylhexyl phthalate, DEHP.
Pregnant CD-1 mice were randomly assigned to receive 0, 250, 500, 1000, or 1500 ppm of DEHP in their feed
during gestational days 0-17. The uterine contents of the mice were examined for toxicity endpoints prior to
normal delivery. The possible outcomes are 1) malformation, 2) death or resorption, 3) no adverse event.
}
\usage{dehp}
\format{  A 'CMData' object, that is a data frame  with the following variables
  \tabular{ll}{
    Trt          \tab factor giving treatment group\cr
    ClusterSize  \tab the size of the litter\cr
    NResp.1        \tab the number of fetuses with a type 1 outcome (malformation) \cr
    NResp.2        \tab the number of fetuses with a type 2 outcome (death or resorption) \cr
    NResp.3        \tab the number of fetuses with a type 3 outcome (normal) \cr
    Freq         \tab the number of litters with the given ClusterSize/NResp.1-NResp.3 combination\cr
  }}
\source{National Toxicology Program, NTP Study TER84064.}
\references{
Tyl, R. W., Price, C. J., Marr, M. C., and Kimmel, C. A. (1988). Developmental toxicity evaluation of dietary
di(2-ethylhexy)phthalate in Fischer 344 rats and CD-1 mice. \emph{Fundamental and Applied Toxicology}
 10, 395-412.
}
\examples{
data(dehp)
library(lattice)
pl <- xyplot(NResp.1/ClusterSize + NResp.2/ClusterSize + NResp.3/ClusterSize ~ Trt, 
        data=dehp, outer=TRUE, type=c("p","a"), jitter.x=TRUE)
pl$condlevels[[1]] <- c("Malformed", "Dead", "Normal")
print(pl)
}
\keyword{datasets}

