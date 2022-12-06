\name{set.finding}
\alias{set.finding}
\title{Set finding}
\description{Set finding}
\usage{
set.finding(dname, node, finding)
}
\arguments{
  \item{dname}{name of gRaven domain}
  \item{node}{name of node}
  \item{finding}{either a state of node \code{node} or a probability distribution over the states}
}
\author{Peter J. Green, \email{P.J.Green@bristol.ac.uk}}

\value{a NULL value is invisibly returned.}

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