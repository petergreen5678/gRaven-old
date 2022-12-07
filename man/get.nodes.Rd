\name{get.nodes}
\alias{get.nodes}
\title{Get the nodes of a gRaven domain}
\description{Get the nodes of a gRaven domain}
\usage{get.nodes(dom)}
\arguments{
  \item{dom}{name of gRaven domain}
}
\details{Emulates function of the same name in the \code{RHugin} package by calls to \code{gRain} functions}
\value{Character vector of node names.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo('chest',package="gRaven",echo=FALSE)
get.nodes(chest)
chest$states
sapply(get.nodes(chest),function(x) get.parents(chest,x))
}
