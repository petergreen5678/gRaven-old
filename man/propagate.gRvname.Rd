\name{propagate.gRvname}
\alias{propagate.gRvname}
\title{Propagate method for a gRaven domain}
\description{Propagate method for a gRaven domain}
\usage{\method{propagate}{gRvname}(object, ...)}
\arguments{
  \item{object}{character string, name of gRaven domain}
  \item{...}{additional arguments to \code{\link{print}}}
}
\value{a NULL value is invisibly returned.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo('chest',echo=FALSE)
compile(chest)
chest
set.finding(chest,"asia","yes")
set.finding(chest,"dysp","no")
propagate(chest)
get.belief(chest,'asia')
get.belief(chest,'tub')
}
