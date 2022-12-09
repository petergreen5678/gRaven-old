\name{list.domains}
\alias{list.domains}
\title{List gRaven domain Objects}
\description{List gRaven domain Objects}
\usage{list.domains(print=TRUE)}
\arguments{
  \item{print}{logical, should list of domains be output to console}
}
\details{gRaven domains are R environment objects, with named slots typically including nodes, parents, cptables and net, where the net component is a \code{gRain} \code{grain} object. The domains that this function finds and lists are all the objects in .GlobalEnv of class "gRv", and all those objects in this class that are components of a component named 'domains' of a list object in .GlobalEnv.}
\value{The function invisibly returns a character vector of the domain names found.}
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
