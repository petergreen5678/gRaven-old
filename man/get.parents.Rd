\name{get.parents}
\alias{get.parents}
\title{Get the parents of a node in a gRaven domain}
\description{Get the parents of a node in a gRaven domain}
\usage{
get.parents(dom, n)
}
\arguments{
  \item{dom}{name of gRaven domain}
  \item{n}{name of node}
}
\value{A named list of character vectors of names of parents.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}


\examples{
demo('chest',package="gRaven",echo=FALSE)
get.nodes(chest)
chest$states
sapply(get.nodes(chest),function(x) get.parents(chest,x))
}
