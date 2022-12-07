\name{delete.node}
\alias{delete.node}
\title{Delete node from a gRaven domain}
\description{Delete node from a gRaven domain}
\usage{
delete.node(dom, n)
}
\arguments{
  \item{dom}{domain name}
  \item{n}{node of to be deleted}
}
\details{Emulates function of the same name in the \code{RHugin} package by calls to \code{gRain} functions}
\value{a NULL value is invisibly returned.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo('chest',package="gRaven",echo=FALSE)
get.nodes(chest)
chest$states
sapply(get.nodes(chest),function(x) get.parents(chest,x))

delete.node(chest,'tub')

get.nodes(chest)
chest$states
sapply(get.nodes(chest),function(x) get.parents(chest,x))
}
