\name{set.finding}
\alias{set.finding}
\alias{retract}
\alias{get.finding}
\title{Set, retract and insspect findings}
\description{Set, retract and insspect findings}
\usage{
set.finding(dom, node, finding)
retract(dom, node=dom$nodes)
get.finding(domain, nodes=domain$nodes, type = c("entered", "propagated"), namestates=TRUE)
}
\arguments{
  \item{dom}{name of gRaven domain}
  \item{node}{name of node}
  \item{nodes}{character vector naming one or mode nodes}
  \item{finding}{either a state of node \code{node} or a probability distribution over the states}
  \item{type}{character vector of one or more options to select between evidence entered before or after most recent propagation}
  \item(namestates}{logical, should output have the states labelled}
}
\details{Unlike with the \code{RHugin} version, evidence can be set at only one node at a time with \code{set.finding}. For \code{type}, partial matching is used.}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}

\value{a NULL value is invisibly returned.}

\examples{
demo(chest,package="gRaven",echo=FALSE)
chest
set.finding(chest,"asia","yes")
set.finding(chest,"dysp","no")
propagate(chest)

get.belief(chest,"asia")
get.belief(chest,"tub")
get.finding(chest)

retract(chest,"dysp")
get.belief(chest,"tub")
get.finding(chest)
}
