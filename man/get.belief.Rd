\name{get.belief}
\alias{get.belief}
\title{Get beliefs in a gRaven domain}
\description{Get beliefs (marginal probabilities) in a gRaven domain}
\usage{get.belief(dom, n)}
\arguments{
  \item{dom}{name of gRaven domain}
  \item{n}{name of node}
}
\details{Emulates function of the same name in the \code{RHugin} package by calls to \code{gRain} functions}
\value{Vector of marginal probabilities}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo('chest',package="gRaven",echo=FALSE)
compile(chest)
chest
set.finding(chest,"asia","yes")
set.finding(chest,"dysp","no")
propagate(chest)
get.belief(chest,'asia')
get.belief(chest,'tub')
}
