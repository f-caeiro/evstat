\name{shape2}
\alias{shape2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of the shape second order (tail) parameter
}
\description{
This function computes the shape second order parameter
}
\usage{
shape2(x, log.x=FALSE, method="fagh", threshold=NA, tau=0, adjust="neg.abs")
}
\arguments{
  \item{x}{a (non-empty) numeric data vector.}
  \item{log.x}{logical; if TRUE, values of \code{x} are given as \code{log(x)}.}
  \item{method}{a character string with the estimator name; see `Details'.}
  \item{threshold}{the threshold(s) used to estimate the second order parameter}
  \item{tau}{a parameter.}
  \item{adjust}{a character string with the type of correction applied in order to have consistent estimates.}
}
\details{
The available methods \code{cg}, \code{fagh} (default) and \code{hw} provide the estimators in Caeiro and Gomes (2014), Fraga Alves \emph{et al.} (2003) and Hall and Welsh (1985), respectively.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
}
\references{
Caeiro, F. and Gomes, M.I. (2014) A semi-parametric estimator of a shape second order parameter. In Pacheco, A., Santos, R., Rosario Oliveira, M. and Paulino, C.D. (Eds.), \emph{New Advances in Statistical Modeling and Applications}, 137--144, Springer.

Fraga Alves, M.I., Gomes, M.I. and de Haan, L. (2003). A new class of
semi-parametric estimators of the second order parameter, \emph{Portugaliae Mathematica} \bold{60}, 193--213.

Hall, P.  and Welsh, A.H.  (1985). Adaptive estimates of parameters of regular variation, \emph{The Annals of Statistics} \bold{13}, 331--341. doi:10.1214/aos/1176346596. \url{http://projecteuclid.org/euclid.aos/1176346596}.
}
\author{
Frederico Caeiro
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}},
\code{\link{adapt.shape2}}
}
%\examples{
%%---- Should be DIRECTLY executable !! ----
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{estimator}
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
