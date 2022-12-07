\name{add.node}
\alias{add.node}
\title{Add Node to a gRaven domain}
\description{Add a node to a gRaven domain.}
\usage{add.node(dom, n, states, subtype)}
\arguments{
  \item{dom}{name of a gRaven domain}
  \item{n}{a character string containing the name of the node to be added}
  \item{states}{a vector listing the states for this node}
  \item{subtype}{a character string, 'labeled', 'numbered' or'boolean'}
}
\details{Emulates function of the same name in the \code{RHugin} package by calls to \code{gRain} functions}
\value{a NULL value is invisibly returned.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo(chest,package="gRaven")
chest
}
