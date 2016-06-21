\name{GMAT}
\alias{GMAT}
\docType{data}
\title{Dichotomous Data Set Based on Graduate Management Admission Test}
\description{
The \code{GMAT} data set is generated data set based on parameters from Graduate Management Admission Test (GMAT) data set (Kingston et al., 1985). First two items were considered to function differently in uniform and non-uniform way respectively. The data set represents responses of 1,000 subjects to multiple-choice test of 20 items. A correct answer is coded as 1 and incorrect answer as 0. The column \code{group} represents group membership, where 0 represents reference group and 1 represents focal group.

}
\usage{data(GMAT)}
\format{
  A \code{GMAT} data frame consists of 1,000 observations on the following 21 variables.

  The first 20 columns represent dichotomously scored items of the test. The 21st column is vector of group membership; values 0 and 1 refer to reference and focal group.
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
Kingston, N., Leary, L., and Wightman, L. (1985). An Exploratory Study of the Applicability of Item Response Theory Methods to the Graduate Management Admission Test. ETS Research Report Series, 1985(2) : 1-64.
}
\examples{
data(GMAT)
}
\keyword{datasets}
