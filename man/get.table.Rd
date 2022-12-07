\name{get.table}
\alias{get.table}
\title{Get a CPT in a gRaven domain}
\description{Get a conditional probability table for a node in a gRaven domain}
\usage{
get.table(dom, n)
}
\arguments{
  \item{dom}{name of gRaven domain}
  \item{n}{name of node}
}
\details{Emulates function of the same name in the \code{RHugin} package by calls to \code{gRain} functions}
\value{The CPT as a data frame.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
chest<-hugin.domain()
add.node(chest,"asia",c("yes","no"))
add.node(chest,"tub",c("yes","no"))
add.edge(chest,"tub","asia")
compile(chest)
chest
get.table(chest,'asia')
}
