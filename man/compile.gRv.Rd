\name{compile.gRv}
\alias{compile.gRv}
\title{Compile a gRaven domain}
\description{Compile a gRaven domain}
\usage{\method{compile}{gRv}(object, ...)}
\arguments{
  \item{object}{name of gRv domain}
  \item{...}{additional arguments to \code{\link{compile}}}
}
\value{a NULL value is invisibly returned.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
chest<-hugin.domain()
add.node(chest,"asia",c("yes","no"))
add.node(chest,"tub",c("yes","no"))
add.edge(chest,"tub","asia")
chest
compile(chest)
chest
set.finding(chest,"asia","yes")
set.finding(chest,"dysp","no")
propagate(chest)
chest$nodes
chest$states
chest$parents
}
