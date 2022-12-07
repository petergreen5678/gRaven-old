\name{add.edge}
\alias{add.edge}
\title{Add Edge(s) to a gRaven domain}
\description{Add a directed edge from one or more parents.}
\usage{add.edge(dom, child, parent)}
\arguments{
  \item{dom}{name of a gRaven domain}
  \item{child}{a character string containing the name of the child node.}
  \item{parent}{a character string specifying by name the parent nodes}
}
\details{Emulates function of the same name in the \code{RHugin} package by calls to \code{gRain} functions}
\value{a NULL value is invisibly returned.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo(chest,package="gRaven")
chest
}
