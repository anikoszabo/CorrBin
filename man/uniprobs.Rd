% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExchMultinomial.R
\name{uniprobs}
\alias{uniprobs}
\title{Extract univariate marginal probabilities from joint probability arrays}
\usage{
uniprobs(jp, type = attr(jp, "type"))
}
\arguments{
\item{jp}{the output of \code{\link{jointprobs}} - a list of joint probability arrays by treatment}

\item{type}{one of c("averaged","cluster","mc") - the type of joint probability. By default,
the \code{type} attribute of \code{jp} is used.}
}
\value{
a list of estimated probability of each outcome by treatment group. The elements are either
matrices or vectors depending on whether cluster-size specific estimates were requested 
(\code{type="cluster"}) or not.
}
\description{
Calculates the marginal probability of each event type for exchangeable correlated multinomial
data based on joint probability estimates calculated by the \code{\link{jointprobs}} function.
}
\examples{
data(dehp)
tau <- jointprobs(dehp, type="averaged")
uniprobs(tau)

#separately for each cluster size
tau2 <- jointprobs(dehp, type="cluster")
uniprobs(tau2)
}
\seealso{
\code{\link{jointprobs}} for calculating the joint probability arrays
}
