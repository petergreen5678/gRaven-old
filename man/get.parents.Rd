\name{get.parents}
\alias{get.parents}
\title{Get the parents of a node in a gRaven domain}
\description{Get the parents of a node in a gRaven domain}
\usage{
get.parents(dname, n)
}
\arguments{
  \item{dname}{name of gRaven domain}
  \item{n}{name of node}
}
\value{A named list of character vectors of names of parents.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}


\examples{
demo('chest',echo=FALSE)
get.nodes(chest)
get(chest)$states
sapply(get.nodes(chest),function(x) get.parents(chest,x))
}
