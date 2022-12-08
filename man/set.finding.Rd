\name{set.finding}
\alias{set.finding}
\alias{retract}
\title{Set and retract finding}
\description{Set findings}
\usage{
set.finding(dom, node, finding)
retract(dom, node=dom$nodes)
}
\arguments{
  \item{dom}{name of gRaven domain}
  \item{node}{name of node}
  \item{nodes}{character vector naming one or mode nodes}
  \item{finding}{either a state of node \code{node} or a probability distribution over the states}
}
\details{Unlike with the \code{RHugin} version, evidence can be set at only one node at a time with \code{set.finding}}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}

\value{a NULL value is invisibly returned.}

\examples{
demo('chest',package="gRaven",echo=FALSE)
compile(chest)
chest
set.finding(chest,"asia","yes")
set.finding(chest,"dysp","no")
propagate(chest)

get.belief(chest,"asia")
get.belief(chest,"tub")

retract(chest,"dysp")
get.belief(chest,"tub")
}
