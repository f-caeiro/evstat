\name{threshold.hill}
\alias{threshold.hill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimators of the extreme value index 
}
\description{
This function computes the threshold and the extreme value index with the Hill estimator. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
threshold.hill(x, method="gh", shape2=NA, scale2=NA, 
ccrit=1.25, plot=FALSE,  legend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a (non-empty) numeric data vector.}
  \item{method}{a character string with the estimator name; see `Details'.}
  \item{shape2}{..}
  \item{scale2}{..}
  \item{ccrit}{critical value used in \code{gh} method; default critical value value is 1.25.}  
  \item{plot}{logical; if TRUE, a graph with the values of the statistic used in \code{gh} method will be plotted.}
  \item{legend}{logical; if TRUE (default) a legend is added to the Plot.}
  \item{\dots}{further arguments to be passed to the Plot.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
method minmse

\deqn{
\hat k_0^H:=  \bigg\lfloor\bigg(\frac{(1-\hat\rho)^2n^{-2\hat\rho}}{-2\hat\rho\hat\beta^2}\bigg)^{1/(1-2\hat\rho)}\bigg\rfloor.}

For a more objective choice of the nuisance parametere \code{ccrit}, see Gomes and Pestana (2011). 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A list with the following components:
\itemize{
    \item \code{threshold}, the threshould estimate.
    \item \code{hill}, the estimate of the Extreme Value Index.
}
}
\references{
Gomes, M.I. and Pestana D. (2011). A Note on the Adaptive Choice of the Optimal Threshold in Extreme Value Analysis. In International Statistical Institute: Proceedings 58th World Statistics Congress, Dublin 4161--4166 (Session CPS014). \url{http://2011.isiproceedings.org/papers/950425.pdf}

Guillou, A. and Hall, P. (2001). A diagnostic for selecting the threshold in extreme-value analysis, \emph{Journal of the Royal Statistical Society  B} \bold{63}:2, 293--305.

Hall, P. (1982). On some simple estimates of an exponent of regular variation, \emph{Journal of the Royal Statistical Society B} \bold{44}, 37--42.
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
\code{\link{evindex}}, \code{\link{shape2}}, \code{\link{scale2}}.
}
\examples{
%%---- Should be DIRECTLY executable !! ----
set.seed(2)
rho <- -0.7
gamma <- 0.3
x<-((runif(500)^rho)-1)^(-gamma/rho)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Algorithm}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
