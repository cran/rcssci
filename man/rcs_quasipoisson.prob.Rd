% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcs_quasipoisson.prob.R
\name{rcs_quasipoisson.prob}
\alias{rcs_quasipoisson.prob}
\title{rcs_quasipoisson.prob}
\arguments{
\item{data}{data.frame.Rdata}

\item{knot}{knot=3-7 or automatic calculate by AIC min}

\item{y}{outcome=0,1}

\item{covs}{covariables, univariate analysis  without "covs" command, multivariable analysis  with "covs" command}

\item{prob}{position parameter,range from 0-1}

\item{x}{main exposure and X-axis when plotting}

\item{filepath}{path of plots output.}
}
\value{
message.print PH assumption and other message
}
\description{
restricted cubic splines (RCS) published in SCI.
}
\details{
quasipoisson models with RCS splines were performed to explore the shape linear or nonlinear(U, inverted U,J,S,L,log,-log,temporary plateau shape)
}
\examples{
library(rcssci)
rcs_quasipoisson.prob(data=sbpdata, y = "status",x = "sbp",
prob=0.1,filepath=tempdir())
# library(rcssci)
# rcs_quasipoisson.prob(knot=4,data=sbpdata, y = "status",x = "sbp",
# covs=c("age","gender"),prob=0.1,filepath="D:/temp")

}
\author{
Zhiqiang Nie, \email{niezhiqiang@gdph.org.cn}
}
