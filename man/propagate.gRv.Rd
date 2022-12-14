\name{propagate.gRv}
\alias{propagate.gRv}
\title{Propagate method for a gRaven domain}
\description{Propagate method for a gRaven domain}
\usage{\method{propagate}{gRv}(object, ...)}
\arguments{
  \item{object}{character string, name of gRaven domain}
  \item{...}{additional arguments to \code{\link{print}}}
}
\value{a NULL value is invisibly returned.}
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
