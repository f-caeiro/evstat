\name{evindex}
\alias{evindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimators of the extreme value index 
}
\description{
This function computes the extreme value index with one of the available semi-parametric estimator. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
evindex(x, log.x=FALSE, method="hill", threshold=NA, 
shape2 = NA, scale2 = NA, param=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a (non-empty) numeric data vector.}
  \item{log.x}{logical; if TRUE, values of \code{x} are given as \code{log(x)}.}
  \item{method}{a character string with the estimator name; see `Details'.}
  \item{threshold}{a vector with the threshold(s) used to estimate the Extreme Value Index.}
  \item{shape2}{NA.} 
  \item{scale2}{NA.}
  \item{param}{a vector with aditional parameters for the Extreme Value Index estimators.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
We shall use the notation 
\eqn{X_{1:n}\le X_{2:n}\le \dots \le X_{n:n}}{X(1:n)\le X(2:n)\le \dots \le X(n:n)} 
for the ascending order statistics of a sample of size \eqn{n}.

Method "\code{hill}" is an implementation of the Hill estimator (Hill, 1975), the average of the log excesses over the high threshold \eqn{X_{n-k:n}}{X(n-k:n)},
\deqn{\hat\xi_{k,n}^{H}=\frac{1}{k}\sum_{i=1}^k (\ln X_{n-i+1:n}-\ln X_{n-k:n}), \qquad k=1,2,\ldots,n-1.}{H(k)=(1/k)\sum ((\ln X(n-i+1:n)-\ln X(n-k:n)), k=1,...,n-1.} 

Method "\code{chill}" give us
the simple class of second-order minimum-variance reduced-bias (MVRB) estimator  in Caeiro et al. (2005), defined by
\deqn{
\hat \xi^{CH}_{k,n}=\hat\xi_{k,n}^{H}\left(1-\frac{\hat\beta(n/k)^{\hat\rho}}{1-\hat\rho}\right),
\qquad k=1,2,\ldots,n-1,}
with  \eqn{\hat\beta, \hat\rho} adequate estimators of the 
`scale' and `shape' second-order parameters  \eqn{(\beta, \rho)}. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
 ...
}
\references{
Brilhante, F., Gomes, M.I. and Pestana, D. (2013). A Simple Generalization of the Hill Estimator. \emph{Computational Statistics and Data Analysis} \bold{57}:1, 518--535.

Caeiro, F., Gomes, M.I. and Pestana, D.D. (2005).  Direct reduction of bias of the classical Hill estimator.  \emph{Revstat} \bold{3}:2, 111--136.

Hill, B. (1975). A simple general approach to inference about the tail of a
distribution, \emph{Ann. Statist.} \bold{3}, 1163--1174.
}
\author{
Frederico Caeiro
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}},
.
}
%\examples{
%%---- Should be DIRECTLY executable !! ----
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ estimator }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
