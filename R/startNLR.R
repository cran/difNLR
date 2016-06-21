###############################################################################
### FUNCTION startNLR() #######################################################
###    13.04.2016       #######################################################
###############################################################################
### calculation of starting values for function difNLR()                    ###
### 3 parametric logistic regression model                                  ###
### INPUT: 'data' matrix or data frame of binary vectors which represents   ###
###               answers of students                                       ###
###        'group' binary vector of group membership,                       ###
###                "1" represents reference group                           ###
###                "0" represents focal group                               ###
###        'parametrization' a character string specifying for which        ###
###                          parametrization starting values should be      ###
###                          calculated. Possible values "IRT" (default)    ###
###                          and "logistic"                                 ###
### OUTPUT: matrix consists of 5 columns                                    ###
###                1st column - discrimination parameter "a"                ###
###                2nd column - difficulty parameter "b"                    ###
###                3rd column - guessing parameter "c"                      ###
###                4th column - difference in discrimination parameter      ###
###                             for reference and focal group               ###
###                5th column - difference in difficulty parameter          ###
###                             for reference and focal group               ###
###############################################################################
###############################################################################
###############################################################################
startNLR <- function(data, group, parameterization = "IRT"){

  startNLS_line <- function(DATA){
    stand_total_score <- c(scale(apply(DATA, 1, sum)))
    ### standardized total score
    Q3 <- cut(stand_total_score, quantile(stand_total_score, (0:3)/3),
              c("I", "II", "III"),
              include.lowest = TRUE)
    ### points
    x <- cbind(mean(stand_total_score[Q3 == "I"]), apply(DATA[Q3 == "I", ], 2, mean))
    y <- cbind(mean(stand_total_score[Q3 == "III"]), apply(DATA[Q3 == "III", ], 2, mean))
    u1 <- y[, 1] - x[, 1]
    u2 <- y[, 2] - x[, 2]
    ### intercept of line
    c <- -(-u1*y[, 2]+ u2*y[, 1])/u1
    ### slope of line
    t <- u2/u1
    ### our line is p(x) = t*x + c, t is slope, c is intercept, p is probability
    results <- as.data.frame(cbind(t, c))
    return(results)
  }
  line <- startNLS_line(data)
  ### guessing parameter influences discrimination and difficulty parameters
  ### we have to apply some correction!

  ### guessing parameter "g"
  ###### as we have new approach in estimating discrimination parameter and
  ###### parameters are connected, we estimate guessing parameters as
  ###### limit p(-infty), but because we have linear function as p(x) we should
  ###### take something finite and negative, e.g. seems reasonable -4 (because we work
  ###### with standardized total score)
  ###### hence: guess = p(-4) = t*(-4) + c
  ###### we consider only non-negative values for guessing parameter
  g <- apply(cbind(0, line$t*(-4) + line$c), 1, max)

  ### now we split data set to two dataset by group and calculate other parameters
  ### for these groups separately
  group <- as.factor(group)
  levels(group) <- c("R", "F")

  data_R <- data[group == "R", ] ### reference group
  data_F <- data[group == "F", ] ### foal group
  line_R <- startNLS_line(data_R)
  line_F <- startNLS_line(data_F)

  ### difficulty parameter "b"
  ###### when no guessing parameter, difficulty parameter "b" is defined as
  ###### p(b) = 1/2
  ###### when considering guessing, we have to apply this correction
  ###### p(b) = (1 + g)/2
  ###### as we have p(x) = t*x + c, then (1 + g)/2 = p(b) = t*b + c
  ###### t*b = (1 + g)/2 - c
  ###### b = ((1 + g)/2 - c)/t

  b_R <- ((1 + g)/2 - line_R$c)/line_R$t
  b_F <- ((1 + g)/2 - line_F$c)/line_F$t
  ### discrimination parameter "a"
  ###### when considering guessing, we have to apply this correction
  ###### p'(b) = a*(1 - g)/4
  ###### as we have p(x) = t*x + c, then p'(x) = t (for each x) and we have to solve
  ###### a*(1 - g)/4 = p'(b) = t
  ###### a = 4*t/(1 - g)
  a_R <- 4*line_R$t/(1 - g)
  a_F <- 4*line_F$t/(1 - g)
  ### difference in difficulties
  results <- switch(parameterization,
                    IRT = data.frame(a_R, b_R, g, a_R - a_F, b_R - b_F),
                    logistic = data.frame(a_R, -a_R*b_R, g, a_R - a_F, -a_R*b_R + a_F*b_F))
  colnames(results) <- c("a", "b", "c", "aDif", "bDif")
  return(results)
}

