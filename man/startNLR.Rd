% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/startNLR.R
\name{startNLR}
\alias{startNLR}
\title{Calculates starting values for non-linear regression DIF models.}
\usage{
startNLR(Data, group, model, match = "zscore", parameterization = "alternative",
         simplify = FALSE)
}
\arguments{
\item{Data}{Data data.frame or matrix: dataset which rows represent
scored examinee answers (\code{"1"} correct, \code{"0"}
incorrect) and columns correspond to the items.}

\item{group}{numeric: binary vector of group membership. \code{"0"}
for reference group, \code{"1"} for focal group.}

\item{model}{character: generalized logistic regression model for
which starting values should be estimated. See \strong{Details}.}

\item{match}{character or numeric: matching criterion to be used as
estimate of trait. Can be either \code{"zscore"} (default,
standardized total score), \code{"score"} (total test score), or
numeric vector of the same length as number of observations in
\code{Data}.}

\item{parameterization}{character: parameterization of regression
coefficients. Possible options are \code{"classic"} (IRT
parameterization), \code{"alternative"} (default) and
\code{"logistic"} (logistic regression). See \strong{Details}.}

\item{simplify}{logical: should initial values be simplified into
the matrix? This is only applicable when parameterization is the
same for all items.}
}
\value{
A list containing elements representing items. Each element is a
named numeric vector of length 8 with initial values for
generalized logistic regression model.
}
\description{
Calculates starting values for \code{difNLR()} function based
on linear approximation.
}
\details{
The unconstrained form of 4PL generalized logistic regression model
for probability of correct answer (i.e., \eqn{y = 1}) is
\deqn{P(y = 1) = (c + cDif * g) + (d + dDif * g - c - cDif * g) / (1 + exp(-(a + aDif * g) * (x - b - bDif * g))), }
where \eqn{x} is by default standardized total score (also called
Z-score) and \eqn{g} is a group membership. Parameters \eqn{a},
\eqn{b}, \eqn{c}, and \eqn{d} are discrimination, difficulty,
guessing, and inattention. Terms \eqn{aDif}, \eqn{bDif},
\eqn{cDif}, and \eqn{dDif} then represent differences between two
groups (reference and focal) in relevant parameters.

The \code{model} argument offers several predefined models. The options are as follows:
\code{Rasch} for 1PL model with discrimination parameter fixed on value 1 for both groups,
\code{1PL} for 1PL model with discrimination parameter fixed for both groups,
\code{2PL} for logistic regression model,
\code{3PLcg} for 3PL model with fixed guessing for both groups,
\code{3PLdg} for 3PL model with fixed inattention for both groups,
\code{3PLc} (alternatively also \code{3PL}) for 3PL regression model with guessing parameter,
\code{3PLd} for 3PL model with inattention parameter,
\code{4PLcgdg} for 4PL model with fixed guessing and inattention parameter for both groups,
\code{4PLcgd} (alternatively also \code{4PLd}) for 4PL model with fixed guessing for both groups,
\code{4PLcdg} (alternatively also \code{4PLc}) for 4PL model with fixed inattention for both groups,
or \code{4PL} for 4PL model.

Three possible parameterization can be specified in
\code{"parameterization"} argument: \code{"classic"} returns IRT
parameters of reference group and differences in these parameters
between reference and focal group. \code{"alternative"} returns IRT
parameters of reference group, the differences in parameters
\code{"a"} and \code{"b"} between two groups and parameters
\code{"c"} and \code{"d"} for focal group. \code{"logistic"}
returns parameters in logistic regression parameterization.
}
\examples{
# loading data
data(GMAT)
Data <- GMAT[, 1:20] # items
group <- GMAT[, "group"] # group membership variable

# starting values for 3PL model
startNLR(Data, group, model = "3PL")

# starting values for 3PL model
# simplified into single table
startNLR(Data, group, model = "3PL", simplify = TRUE)

# starting values for 3PL model
# with score as matching criterion
startNLR(Data, group, model = "3PL", match = "score")

# starting values for model specified for each item
startNLR(Data, group,
  model = c(
    rep("1PL", 5), rep("2PL", 5),
    rep("3PL", 5), rep("4PL", 5)
  )
)
}
\references{
Drabinova, A. & Martinkova, P. (2017). Detection of differential item functioning with nonlinear regression:
A non-IRT approach accounting for guessing. Journal of Educational Measurement, 54(4), 498--517,
\doi{10.1111/jedm.12158}.

Hladka, A. (2021). Statistical models for detection of differential item functioning. Dissertation thesis.
Faculty of Mathematics and Physics, Charles University.
}
\seealso{
\code{\link[difNLR]{difNLR}}
}
\author{
Adela Hladka (nee Drabinova) \cr
Institute of Computer Science of the Czech Academy of Sciences \cr
Faculty of Mathematics and Physics, Charles University \cr
\email{hladka@cs.cas.cz} \cr

Patricia Martinkova \cr
Institute of Computer Science of the Czech Academy of Sciences \cr
\email{martinkova@cs.cas.cz} \cr
}
