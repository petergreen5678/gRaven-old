\name{list.domains}
\alias{list.domains}
\title{List gRaven domain Objects}
\description{List gRaven domain Objects, and optionally output or tidy this list.}
\usage{list.domains(print=TRUE,clean=FALSE)}
\arguments{
  \item{clean}{logical, should domains considered to be "orphans" be deleted, see details.}
  \item{print}{logical, should list of domains be output to console}
}
\details{gRaven domains are R list objects, with named slots typically including nodes, parents, cptables and net, where the net component is a \code{gRain} \code{grain} object. These domains are by convention in the .GlobalEnv environment, and have names beginning ".gRvd". The domains that this function finds and lists are all the objects in .GlobalEnv of class "gRv" (whatever their names). The function also lists "owners" of such objects. These owners are R list objects in .GlobalEnv that have a slot named \code{domains}, which is itself a list of gRaven domains.}
\value{The function invisibly returns a character vector of the domain names found.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}
\examples{
demo(chest)
chest
get(chest)
chest2<-clone.domain(chest)
get(chest2)

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