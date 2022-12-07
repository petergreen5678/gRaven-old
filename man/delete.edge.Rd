\name{delete.edge}
\alias{delete.edge}
\title{Delete edge from a gRaven domain}
\description{
Delete edge from a gRaven domain
}
\usage{delete.edge(dom, n, p)}
\arguments{
  \item{dom}{domain name}
  \item{n}{child node of edge to be deleted}
  \item{p}{parent node of edge to be deleted}
}
\value{a NULL value is invisibly returned.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo('chest',package="gRaven",echo=FALSE)
get.nodes(chest)
chest$states
sapply(get.nodes(chest),function(x) get.parents(chest,x))

delete.edge(chest,"xray","either")

get.nodes(chest)
chest$states
sapply(get.nodes(chest),function(x) get.parents(chest,x))
}
