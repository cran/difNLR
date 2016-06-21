\name{GMATtest}
\alias{GMATtest}
\docType{data}
\title{Data Set Based on Graduate Management Admission Test}
\description{
The \code{GMATtest} data set is generated data set based on parameters from Graduate Management Admission Test (GMAT) data set (Kingston et al., 1985). First two items were considered to function differently in uniform and non-uniform way respectively. The data set represents responses of 1,000 subjects to multiple-choice test of 20 items. Aditionally, 4 possible answers on all items were generated, coded A, B, C and D. The column \code{group} represents group membership, where 0 represents reference group and 1 represent focal group.
}
\usage{data(GMATtest)}
\format{
  A \code{GMAT} data frame consists of 1,000 observations on the following 21 variables.

  The first 20 columns represents answers of subject to an items of the test. The 21st column is vector of group membership; values 0 and 1 refer to reference and focal group.

  Correct answers were D, C, C, C, B, B, B, A, A, D, A, A, A, C, C, B, C, B, D, A respectively. For more details see data \code{GMATkey}

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
\seealso{
\code{\link{GMATkey}}, \code{\link{GMAT}}
}
\examples{
data(GMATtest)
}
\keyword{datasets}
