\name{clone.domain}
\alias{clone.domain}
\title{Clone a gRaven domain Object}
\description{Create a gRaven domain.}
\usage{clone.domain(dom1)}
\arguments{
  \item{dom1}{name of gRaven domain}
}
\details{Makes a copy of an existing gRaven domain}
\value{character string naming a gRaven domain}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo(chest,package="gRaven")
chest
chest2<-clone.domain(chest)
chest2

set.finding(chest,"asia","yes")
set.finding(chest,"dysp","no")
propagate(chest)
get.belief(chest,'asia')
get.belief(chest,'tub')

propagate(chest2)
get.belief(chest2,'asia')
get.belief(chest2,'tub')

list.domains()
}
