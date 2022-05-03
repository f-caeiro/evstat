\name{adapt.shape2}
\alias{adapt.shape2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adaptative estimation of the shape second order (tail) parameter
}
\description{
This function computes the adaptative estimate of the shape second order parameter
}
\usage{
adapt.shape2(x,
  log.x=FALSE,
  shape2=c("fagh","cg"), 
  threshold=floor(length(x[x>0])^0.999), 
  range=floor(length(x[x>0])^(c(0.995,0.999))), 
  param=c(0,1), 
  adjust="neg.abs"
  )
}
\arguments{
  \item{x}{a (non-empty) numeric data vector.}
  \item{log.x}{logical; if TRUE, values of \code{x} are given as \code{log(x)}.}
  \item{shape2}{a character string with the estimator of the shape second order parameter.}
  \item{threshold}{the threshold used to estimate the second order parameter.}
  \item{range}{range of thresholds used to selected.}
  \item{param}{a vector of tuning parameters.}
  \item{adjust}{a character string with the type of correction applied in order to have consistent estimates.}
}
\details{
The available \code{shape2} estimators are \code{cg} and \code{fagh} (default), the estimators introduced in Caeiro and Gomes (2014) and Fraga Alves et al. (2003), respectively.
}
\value{
.
}
\references{
Gomes, M.I. and Pestana, D.  (2007). A sturdy reduced bias extreme quantile (VaR) estimator, \emph{J. American Statistical Association}, \bold{102}:477,  280--292. 
}
\author{
Frederico Caeiro
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}},
\code{\link{shape2}}
}
%\examples{
%%---- Should be DIRECTLY executable !! ----
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Algorithm}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
