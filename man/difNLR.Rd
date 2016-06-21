\name{difNLR}
\alias{difNLR}
\alias{print.difNLR}
\alias{plot.difNLR}
\alias{fitted.difNLR}
\alias{predict.difNLR}
\title{Non-Linear Regression DIF Method}
\description{Performs DIF detection using Non-Linear Regression method.}
\usage{
difNLR(data, group,  type = "both", p.adjust.method = "BH", start,...)

\method{print}{difNLR}(x, ...)
\method{plot}{difNLR}(x, plot.type = "cc", item = "all", col = c("tomato", "turquoise"),
alpha = .5, shape = 21, size = .8, linetype = c(1, 2), title, ...)
\method{fitted}{difNLR}(object, item = "all", ...)
\method{predict}{difNLR}(object, item = "all", score, group,...)
}
\arguments{
  \item{data}{numeric: binary data matrix. See \bold{Details}.}
  \item{group}{numeric: binary vector of group membership. "0" for reference group, "1" for focal group.}
  \item{type}{character: type of DIF to be tested (either "both" (default), "udif" or "nudif"). See \bold{Details}.}
  \item{p.adjust.method}{character: method for multiple comparison correction. See \bold{Details}.}
  \item{start}{numeric: matrix with n rows (where n is the number of items) and at most 5 columns containing initial item parameters estimates. See \bold{Details}.}
  \item{x, object}{an object of \code{"difNLR"} class}
  \item{plot.type}{character: type of plot to be plotted (either \code{"cc"} for characteristic curve, or "stat" for F-test statistics). See \bold{Details}.}
  \item{item}{either character (\code{"all"}), or numeric vector, or single number corresponding to column indicators. See \bold{Details}.}
  \item{col}{character: single value, or vector of two values representing colors for plot.}
  \item{alpha}{numeric: transparency parameter for plot.}
  \item{shape}{integer: shape parameter for plot.}
  \item{size}{numeric: line width in plot.}
  \item{linetype}{line type for reference and focal group in plot.}
  \item{title}{string corresponding to title of plot.}
  \item{score}{numeric: standardized total score of subject.}
  \item{\dots}{ not used. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
DIF detection procedure based on Non-Linear Regression is the extension of Logistic Regression procedure (Swaminathan and Rogers, 1990).

The \code{data} is a matrix whose rows represents examinee answers (\code{"1"} correct, \code{"0"} incorrect) and columns correspond to the items. The \code{group} must be a vector of the same length as \code{nrow(data)}.

The \code{type} corresponds to type of DIF to be tested. Possible values are \code{"both"} to detect any DIF (uniform and/or non-uniform), \code{"udif"} to detect only uniform DIF or \code{"nudif"} to detect only non-uniform DIF.

The \code{start} is a matrix with a number of rows equal to number of items. The number of columns correspond to number of parameters in model in alternative hypothesis (5 for values \code{"both"} and \code{"nudif"} in \code{type}, 4 for \code{"udif"} in \code{type}). If \code{start} missing, initial values are calculated by \code{startNLR()} function.

The \code{p.adjust.method} is a character for \code{p.adjust} function from the \code{stats} package. Possible values are \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.


The output of the \code{difNLR} is displayed by the \code{print.difNLR} function.


Two types of plots are available.  The first one is obtained by setting \code{plot.type = "cc"}(default). The characteristic curve for item specified in \code{item} option is plotted. For default option \code{"all"} of \code{item}, characteristic curves of all converged items are plotted. The drawn curves represent best model.

The second plot is obtained by setting \code{plot.type = "stat"}. The F-statistics are displayed on the Y axis, for each coverged item. The detection threshold is displayed by a horizontal line and items detected as DIF are printed with the red color. Only parameters \code{size} and \code{title} are used.

Fitted values are extracted by the \code{fitted.difNLR} function for item(s) specified in \code{item} argument.

Predicted values are produced by the \code{predict.difNLR} function for item(s) specified in \code{item} argument. \code{score} represents standardized total score of new subject and \code{group} argument represents group membership of new subject.

Missing values are not allowed either for responses or for group membership.
}
\value{
A list of class \code{"difNLR"} with the following arguments:
\item{DIF}{either the column indicators of the items which were detected as DIF, or \code{"NONE"}.}
\item{Fval}{the values of F-test statistics.}
\item{Pval }{the p-values by F-test.}
\item{df }{the degrees of freedom of F-test.}
\item{coef }{the matrix of estimated item parameters.}
\item{group }{the vector of group membership.}
\item{data }{the binary data matrix.}
\item{type }{character: type of DIF that was tested.}
\item{conv_fail }{numeric: number of convergence issues.}
\item{conv_fail_which }{the column indicators of the items which did not converge.}
\item{p.adjust.method }{character: method for multiple comparison correction which was applied.}
}
\references{
Drabinova, A. and Martinkova P. (2016). Detection of Differenctial Item Functioning Based on Non-Linear Regression, Technical Report, V-1229, \url{http://hdl.handle.net/11104/0259498}

Swaminathan, H. and Rogers, H. J. (1990). Detecting Differential Item Functioning Using Logistic Regression Procedures. Journal of Educational Measurement, 27, 361-370.
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
\seealso{
\code{\link{p.adjust}}
}
\examples{
\dontrun{
# loading data based on GMAT
data(GMAT)

data  <- GMAT[, colnames(GMAT) != "group"]
group <- GMAT[, "group"]

# Testing both DIF effects simultaneously
x <- difNLR(data, group)

# Testing both DIF effects with none multiple comparison correction
difNLR(data, group, type = "both", p.adjust.method = "none")

# Testing uniform DIF effects
difNLR(data, group, type = "udif")

# Testing non-uniform DIF effects
difNLR(data, group, type = "nudif")

# Graphical devices
plot(x)
plot(x, item = x$DIF)
plot(x, plot.type = "stat")

# Fitted values
fitted(x)
fitted(x, item = 1)

# Predicted values
predict(x)
predict(x, item = 1)

# Predicted values for new subjects
predict(x, item = 1, score = 0, group = 1)
predict(x, item = 1, score = 0, group = 0)
}}
\keyword{DIF}

