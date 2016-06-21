\name{startNLR}
\alias{startNLR}
\title{Calculates Starting Values}
\usage{
startNLR(data, group, parameterization = "IRT")
}
\arguments{
  \item{data}{a matrix or data frame of binary data.}
  \item{group}{a binary vector of group membership. "0" represents reference group, "1" represents focal group.}
  \item{parameterization}{character: parameterization of regression coefficients. Possible options are "IRT" (Item response theory, default option) and "LR" (Logistic regression)}
}
\description{
Calculates starting values for function \code{difNLR} function based on linear approximation.
}
\author{Adela Drabinova \cr
Institute of Computer Science, The Czech Academy of Sciences \cr
Faculty of Mathematics and Physics, Charles University in Prague \cr
adela.drabinova@gmail.com\cr\cr
Patricia Martinkova \cr
Institute of Computer Science, The Czech Academy of Sciences \cr
martinkova@cs.cas.cz\cr\cr
Karel Zvara \cr
Faculty of Mathematics and Physics, Charles University in Prague\cr
}
\references{
Drabinova, A. and Martinkova P. (2016). Detection of Differenctial Item Functioning Based on Non-Linear Regression, Technical Report, V-1229, \url{http://hdl.handle.net/11104/0259498}
}
\examples{
# loading data based on GMAT
data(GMAT)

data  <- GMAT[, colnames(GMAT) != "group"]
group <- GMAT[, "group"]

# starting values
startNLR(data, group)
}
