
###############################################################################
### FUNCTION difNLR() #########################################################
###   02.05.2016      #########################################################
###############################################################################
### detection of DIF by function NLS - non-linear regression model          ###
### 3 parametric logistic regression model                                  ###
### INPUT: 'data' matrix or data frame of binary vectors which represents   ###
###               answers of students                                       ###
###        'group' binary vector of group membership,                       ###
###                "0" represents reference group                           ###
###                "1" represents focal group                               ###
###        'start' matrix or data frame of starting values                  ###
###                must containt of 5 columns (for type = "both"            ###
###                or "nudif"), each column stands for one parameter        ###
###                1st column - discrimination parameter "a"                ###
###                2nd column - difficulty parameter "b"                    ###
###                3rd column - guessing parameter "c"                      ###
###                4th column - difference in discrimination parameter      ###
###                             for reference and focal group               ###
###                5th column - difference in difficulty parameter          ###
###                             for reference and focal group               ###
###                for type = "udif" 4 columns are required. 4th column     ###
###                stands for difference in difficulty parameter for        ###
###                reference and focal group. Rows represents items.        ###
###                If no starting values are obtained, they are calculated  ###
###                function startNLS()                                      ###
###         'type' a character string specifying which DIF effects          ###
###                must be tested. Possible values are "both" (default),    ###
###                "udif" and "nudif"                                       ###
###         'p.adjust.method' correction method for adjusting p-values      ###
###                           Possible values are "holm", "hochberg"        ###
###                           (default), "hommel", "bonferroni", "BH",      ###
###                           "BY", "fdr", "none"                           ###
### OUTPUT: list of 5 elements                                              ###
###         'DIF' numeric vector of items which are detected as DIF         ###
###               If there is no DIF, value of 'DIF' is character string    ###
###               "NONE"                                                    ###
###         'Fval' F-value of ANOVA test of model and submodel              ###
###         'Pval' P-value of ANOVA test after multiple comparison          ###
###                correction (Hochberg method)                             ###
###         'coef' list of coeficients for model or submodel,               ###
###                one that fits data the best                              ###
###         'conv_fail' counts of items with convergence failure            ###
###         'group' as in input                                             ###
###         'data'  as in input                                             ###
###         'type'  as in input                                             ###
###############################################################################
###############################################################################
###############################################################################




difNLR <- function(data, group, type = "both", p.adjust.method = "BH", start, ...)
{
  if (length(levels(as.factor(group))) != 2)
    stop("'group' must be binary vector",
         call. = FALSE)
  if (length(levels(as.factor(group))) != 2)
    stop("'group' must be binary vector",
         call. = FALSE)
  if(is.matrix(data) | is.data.frame(data)){
    if (!all(apply(data, 2, function(i) {length(levels(as.factor(i))) == 2})))
      stop("'data' must be data frame or matrix of binary vectors",
           call. = FALSE)
    if(nrow(data) != length(group))
      stop("'data' must have the same number of rows as is length of vector 'group'",
           call. = FALSE)
  } else {
    stop("'data' must be data frame or matrix of binary vectors",
         call. = FALSE)
  }
  if (!type %in% c("udif", "nudif", "both") | !is.character(type))
    stop("'type' must be either 'udif', 'nudif' or 'both'",
         call. = FALSE)
  ### STARTING VALUES
  if (missing(start)){
    start <- startNLR(data, group)
    start_m0 <- switch(type, both = start, nudif = start, udif = start[, -4])
    start_m1 <- switch(type, both = start[, -c(4,5)],
                       nudif = start[, -4], udif = start[, -c(4,5)])
  } else {
    if (ncol(start) != 5 & type != "udif")
      stop("'start' must be data frame or matrix with 5 columns",
           call. = FALSE)
    if (ncol(start) != 4 & type == "udif")
      stop("'start' must be data frame or matrix with 4 columns for detecting uniform DIF",
           call. = FALSE)
    if (nrow(start) != ncol(data))
      stop("'start' must be data frame or matrix with starting values for each item",
           call. = FALSE)
    start_m0 <- start
    start_m1 <- switch(type, both = start[, -c(4, 5)],
                       nudif = start[, -4], udif = start[, -4])
  }
  ### NA values omitted
  DATA <- data.frame(data, group)
  DATA <- na.omit(DATA)
  data  <- DATA[, !(names(DATA) %in% "group")]
  group <- DATA$group

  ### MODELS
  ### reg function for non-uniform DIF
  regFce_nonUDIF <- deriv3( ~ c + (1 - c)/(1 + exp(-(a + aDif*group)*(x - (b + bDif*group)))),
                            namevec = c("a", "b", "c", "aDif", "bDif"),
                            function.arg = function(x, group, a, b, c, aDif, bDif){})
  ### reg function for uniform DIF
  regFce_UDIF <- deriv3( ~ c + (1 - c)/(1 + exp(-a*(x - (b + bDif*group)))),
                         namevec = c("a", "b", "c", "bDif"),
                         function.arg = function(x, group, a, b, c, bDif){})
  ### reg function for NO DIF
  regFce_noDIF <- deriv3( ~ c + (1 - c)/(1 + exp(-a*(x - b))),
                          namevec = c("a", "b", "c"),
                          function.arg = function(x, group, a, b, c){})
  ### STANDARDIZED TOTAL SCORE
  stand_total_score <- c(scale(apply(data, 1, sum)))
  m                 <- ncol(data)
  ### number of items
  n                 <- nrow(data)
  ### number of observations

  ### ESTIMATES
  ### model m0
  conv_fail <- 0
  k         <- Inf
  estim_m0  <- lapply(1:m, function(i) NA)
  hv <- which(!(sapply(1:m, function(i) is(try(switch(type, both =  nls(data[, i] ~ regFce_nonUDIF(stand_total_score, group, a, b, c, aDif, bDif),
                                                                        algorithm = "port", start = start_m0[i, ],
                                                                        lower = c(-k, -k, 0, -k, -k), upper = c(k, k, 1, k, k)),
                                                      nudif = nls(data[, i] ~ regFce_nonUDIF(stand_total_score, group, a, b, c, aDif, bDif),
                                                                  algorithm = "port", start = start_m0[i, ],
                                                                  lower = c(-k, -k, 0, -k, -k), upper = c(k, k, 1, k, k)),
                                                      udif = nls(data[, i] ~ regFce_UDIF(stand_total_score, group, a, b, c, bDif),
                                                                 algorithm = "port", start = start_m0[i, ],
                                                                 lower = c(-k, -k, 0, -k), upper = c(k, k, 1, k))), silent = T), "try-error"))))
  estim_m0[hv] <- lapply(hv, function(i) switch(type, both =  nls(data[, i] ~ regFce_nonUDIF(stand_total_score, group, a, b, c, aDif, bDif),
                                                                  algorithm = "port", start = start_m0[i, ],
                                                                  lower = c(-k, -k, 0, -k, -k), upper = c(k, k, 1, k, k)),
                                                nudif = nls(data[, i] ~ regFce_nonUDIF(stand_total_score, group, a, b, c, aDif, bDif),
                                                            algorithm = "port", start = start_m0[i, ],
                                                            lower = c(-k, -k, 0, -k, -k), upper = c(k, k, 1, k, k)),
                                                udif = nls(data[, i] ~ regFce_UDIF(stand_total_score, group, a, b, c, bDif),
                                                           algorithm = "port", start = start_m0[i, ],
                                                           lower = c(-k, -k, 0, -k), upper = c(k, k, 1, k))))
  ### model m1
  estim_m1 <- lapply(1:m, function(i) NA)
  hv <- which(!(sapply(1:m, function(i) is(try(switch(type, both =  nls(data[, i] ~ regFce_noDIF(stand_total_score, group, a, b, c),
                                                                        algorithm = "port", start = start_m1[i, ],
                                                                        lower = c(-k, -k, 0), upper = c(k, k, 1)),
                                                      nudif = nls(data[, i] ~ regFce_UDIF(stand_total_score, group, a, b, c, bDif),
                                                                  algorithm = "port", start = start_m1[i, ],
                                                                  lower = c(-k, -k, 0, -k), upper = c(k, k, 1, k)),
                                                      udif = nls(data[, i] ~ regFce_noDIF(stand_total_score, group, a, b, c),
                                                                 algorithm = "port", start = start_m1[i, ],
                                                                 lower = c(-k, -k, 0), upper = c(k, k, 1))), silent = T), "try-error"))))
  estim_m1[hv] <- lapply(hv, function(i) switch(type, both =  nls(data[, i] ~ regFce_noDIF(stand_total_score, group, a, b, c),
                                                                  algorithm = "port", start = start_m1[i, ],
                                                                  lower = c(-k, -k, 0), upper = c(k, k, 1)),
                                                nudif = nls(data[, i] ~ regFce_UDIF(stand_total_score, group, a, b, c, bDif),
                                                            algorithm = "port", start = start_m1[i, ],
                                                            lower = c(-k, -k, 0, -k), upper = c(k, k, 1, k)),
                                                udif = nls(data[, i] ~ regFce_noDIF(stand_total_score, group, a, b, c),
                                                           algorithm = "port", start = start_m1[i, ],
                                                           lower = c(-k, -k, 0), upper = c(k, k, 1))))

  conv_fail <- conv_fail + sum(is.na(estim_m1) | is.na(estim_m0))
  conv_fail_which <- which(is.na(estim_m1) | is.na(estim_m0))
  if (conv_fail > 0){
    warning("Convergence failure")
  }
  pval <- Fval <- rep(NA, m)

  # pval[which(!((is.na(estim_m1))|(is.na(estim_m0))))] <- sapply(which(!((is.na(estim_m1))|(is.na(estim_m0)))),
  #                                                               function(l) anova(estim_m1[[l]], estim_m0[[l]])$Pr[2])
  # Fval[which(!((is.na(estim_m1))|(is.na(estim_m0))))] <- sapply(which(!((is.na(estim_m1))|(is.na(estim_m0)))),
  #                                                               function(l) anova(estim_m1[[l]], estim_m0[[l]])$F[2])

  df <- switch(type, both = c(2, n - 5), udif = c(1, n - 4), nudif = c(1, n - 5))
  Fval[which(!((is.na(estim_m1))|(is.na(estim_m0))))] <- sapply(which(!((is.na(estim_m1))|(is.na(estim_m0)))),
                                                                function(l) ((estim_m1[[l]]$m$deviance() -
                                                                                estim_m0[[l]]$m$deviance()) / df[1]) / (estim_m0[[l]]$m$deviance()/df[2]))
  pval[which(!((is.na(estim_m1))|(is.na(estim_m0))))] <- sapply(which(!((is.na(estim_m1))|(is.na(estim_m0)))),
                                                                function(l) (1 - pf(Fval[l], df[1], df[2])))
  # critical value
  # Fval_critical <- qf(0.95, df[1], df[2])

  pval_adj <- p.adjust(pval, method = p.adjust.method)
  significant <- which(pval_adj < 0.05)


  if ( length(significant) > 0 ) {
    DIF <- significant
  } else {
    DIF <- "NONE"
  }
  coef <- ifelse(pval_adj < 0.05, lapply(estim_m0[which(!is.na(estim_m0))], coef), lapply(estim_m1[which(!is.na(estim_m1))], coef))
  mat <- t(sapply(coef, "[", i = 1:5))
  mat[is.na(mat)] <- 0
  coef <- switch(type, both = mat, nudif = mat[, c(1:3, 5, 4)], udif = mat[, c(1:3, 5, 4)])
  colnames(coef) <- c(letters[1:3], "aDif", "bDif")


  ### results
  results <- list(DIF = DIF, Fval = Fval, Pval = pval_adj, df = df, coef = coef, group = group,
                  data = data, conv_fail = conv_fail, type = type, conv_fail_which = conv_fail_which,
                  p.adjust.method = p.adjust.method)
  class(results) <- "difNLR"
  results
}


###############################################################################
### FUNCTION print.difNLR() ###################################################
###       18.04.2016        ###################################################
###############################################################################
### prints information about DIF detection of object of class difNLS        ###
### INPUT: 'x' object of class difNLS                                       ###
###############################################################################
print.difNLR <- function(x, ...){
  ### printing output
  title <- switch(x$type, both = "Detection of both types of Differential Item Functioning using Non-Linear Regression method",
                  udif = "Detection of uniform Differential Item Functioning using Non-Linear Regression method",
                  nudif = "Detection of non-uniform Differential Item Functioning using Non-Linear Regression method")
  cat(paste(strwrap(title, width = 60), collapse = "\n"))
  cat("\n\nNon-linear regression DIF statistics\n")
  if (x$p.adjust.method == "none"){
    cat("p-values with none multiple comparison correction\n\n")
  } else {
    cat(paste("p-values adjusted using", switch(x$p.adjust.method,
                                                holm = "Holm", hochberg = "Hochberg",
                                                hommel = "Hommel", bonferroni = "Bonferroni",
                                                BH = "Benjamini-Hochberg", BY = "Benjamini-Yekutieli",
                                                fdr = "FDR"), "method\n\n"))
  }

  tab  <- format(round(cbind(x$Fval, x$Pval), 4))
  sign <- ifelse(is.na(x$Pval), " ", ifelse(x$Pval < 0.001, "***",
                                            ifelse(x$Pval < 0.01, "**",
                                                   ifelse(x$Pval < 0.05, "*",
                                                          ifelse(x$Pval < 0.1, ".", " ")))))
  tab <- matrix(cbind(tab, sign), ncol = 3)
  rownames(tab) <- paste("Item", 1:length(x$Pval))
  colnames(tab) <- c("F-value", "P-value", "")
  print(tab, quote = F, digits = 4, zero.print = F)
  cat("\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  n <- nrow(x$data)
  df <- switch(x$type, both = c(2, n - 5), udif = c(1, n - 4), nudif = c(1, n - 5))
  Fval_critical <- qf(0.95, df[1], df[2])
  cat(paste("\nDetection thresholds:", round(Fval_critical, 3)))
  if (is.character(x$DIF)){
    switch(x$type, both = cat("\nNone of items is detected as DIF"),
           udif = cat("\nNone of items is detected as uniform DIF"),
           nudif = cat("\nNone of items is detected as non-uniform DIF"))
  } else {
    switch(x$type, both = cat("\nItems detected as DIF items:\n"),
           udif = cat("\nItems detected as uniform DIF items:\n"),
           nudif = cat("\nItems detected as non-uniform DIF items:\n"))
    cat("\n", paste("Item ", x$DIF, "\n", sep = ""))
  }
}


###############################################################################
### FUNCTION plot.difNLR() ####################################################
###      18.04.2016        ####################################################
###############################################################################
### plots of characteristic curves from difNLR() usign ggplot2              ###
### INPUT: 'x'     an object of 'difNLS'                                    ###
###        'type'  type of plot, possible values 'cc' for characteristic    ###
###                curve, 'stat' for F-statistics                           ###
###  Further input options only for 'cc' type, ignored if type 'stat'       ###
###        'item'  character string or an object of class numeric or        ###
###                integer specifying for which items characteristic curves ###
###                should be plotted                                        ###
###                Possible values "all" (default), numeric vector of items ###
###                or integer representing item                             ###
###        'col'   vector of two values or single integer of colors         ###
###        'alpha' single value of transparency parameter for ggplot        ###
###        'shape' single value of shape of points for ggplot               ###
###        'size'  single value of thickness of lines for ggplot            ###
###        'linetype' vector of two values or single integer/string for     ###
###                   line type for ggplot                                  ###
###############################################################################
###############################################################################
###############################################################################
### packages
#   list.of.packages <- c("ggplot2", "grid")
#   new.packages     <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
#   if(length(new.packages)) install.packages(new.packages)
#   lapply(list.of.packages, require, character.only = TRUE)



plot.difNLR <- function(x, plot.type = "cc", item = "all",
                        col = c("tomato", "turquoise"), alpha = .5, shape = 21, size = .8,
                        linetype = c(1, 2), title,
                        ...){
  plotstat <- function(x, size = size, title = title){
    if (x$conv_fail != 0){
      if (length(x$conv_fail) == length(x$Fval)){
        stop("None of items does converge. F-statistic values not plotted",
             call. = FALSE)
      } else {
        warning(paste("Item ", x$conv_fail_which, " does not converge. F-statistic value not plotted",
                      sep = "", collapse = "\n"),
                call. = FALSE)
      }
    }

    if(missing(title)){
      title <- "Non-Linear Regression DIF Detection \n with None Multiple Comparison Correction"
    }

    n <- nrow(x$data)
    df <- switch(x$type, both = c(2, n - 5), udif = c(1, n - 4), nudif = c(1, n - 5))
    Fval_critical <- qf(0.95, df[1], df[2])
    g <- as.factor(ifelse(x$Fval > Fval_critical, 1, 0))
    items <- setdiff(1:length(x$Fval), x$conv_fail_which)
    hv <- na.omit(as.data.frame(cbind(1:length(x$Fval), x$Fval, g)))
    plot_stat <- ggplot(hv, aes_string(x = "V1", y = "V2",
                                       label = as.character("V1"), col = as.factor(g))) +
      ### points
      geom_text() +
      scale_color_manual(values = c("black", "red")) +
      ### critical value
      geom_hline(yintercept = Fval_critical, size = size) +
      ### theme
      ggtitle(title) +
      labs(x = "Item", y = "F-statistic") +
      theme_bw() +
      theme(text = element_text(size = 10),
            plot.title = element_text(size = 10, face = "bold", vjust = 1.5),
            axis.line  = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.position = "none")

    print(plot_stat)
  }

  plotCC <- function(x, item = item,
                     col = col, alpha = alpha, shape = shape, size = size,
                     linetype = linetype, title = title){
    m <- nrow(x$coef)
    if (class(item) == "character"){
      if (item != "all")
        stop("'item' must be either numeric vector or character string 'all' ",
             call. = FALSE)
    } else {
      if (class(item) != "integer" & class(item) != "numeric")
        stop("'item' must be either numeric vector or character string 'all' ",
             call. = FALSE)
    }
    if (class(item) == "numeric" & !all(item %in% 1:m))
      stop("invalid number of 'item'",
           call. = FALSE)
    if (class(item) == "integer" & !all(item %in% 1:m))
      stop("'item' must be either numeric vector or character string 'all' ",
           call. = FALSE)
    if (length(col) == 1){
      col <- rep(col, 2)
    } else {
      if (length(col) > 2){
        col <- col[1:2]
        warning("Length of 'col' is greater than 2. Only first two values are used",
                call. = FALSE)
      }
    }
    if (class(item) == "character"){
      items <- 1:m
    } else {
      items <- item
    }
    if (any(x$conv_fail_which %in% items)){
      if (length(setdiff(items, x$conv_fail_which)) == 0){
        stop(paste("Item ", intersect(x$conv_fail_which, items), " does not converge. Characteristic curve not plotted",
                   sep = "", collapse = "\n"),
             call. = FALSE)
      } else {
        warning(paste("Item ", intersect(x$conv_fail_which, items), " does not converge. Characteristic curve not plotted",
                      sep = "", collapse = "\n"),
                call. = FALSE)
        items <- setdiff(items, x$conv_fail_which)
      }
    }
    if (length(linetype) != 2){
      if (length(linetype) == 1){
        linetype <- rep(linetype, 2)
      } else {
        linetype <- linetype[1:2]
        warning("Length of 'linetype' is greater than 2. Only first two values are used",
                call. = FALSE)
      }
    }
    if (!missing(title)){
      TITLE <- title
    }

    ### functions
    regFce_nonUDIF <- function(STS, group, a, b, c, aDif, bDif){
      return(c + (1 - c)/(1 + exp(-(a + aDif*group)*(STS - (b + bDif*group)))))
    }

    ### data
    stand_total_score_R <- c(scale(apply(x$data[x$group == 0, ], 1, sum)))
    stand_total_score_F <- c(scale(apply(x$data[x$group == 1, ], 1, sum)))
    max_sts <- max(as.numeric(levels(as.factor(stand_total_score_R))),
                   as.numeric(levels(as.factor(stand_total_score_F))))
    min_sts <- min(as.numeric(levels(as.factor(stand_total_score_R))),
                   as.numeric(levels(as.factor(stand_total_score_F))))

    for (i in items){
      hv_R <- data.frame(cbind(as.numeric(levels(as.factor(stand_total_score_R))),
                               tapply(x$data[x$group == 0, i], as.factor(stand_total_score_R), mean)))
      hv_F <- data.frame(cbind(as.numeric(levels(as.factor(stand_total_score_F))),
                               tapply(x$data[x$group == 1, i], as.factor(stand_total_score_F), mean)))
      hv   <- data.frame(rbind(cbind(hv_R, Group = "Reference"), cbind(hv_F, Group = "Focal")))
      rownames(hv) <- 1:dim(hv)[1]
      hv$size <- c(table(stand_total_score_R), table(stand_total_score_F))
      if (missing(title)){
        TITLE <- paste("Item", i)
      }
      plot_CC <- ggplot(hv, aes_string("X1", "X2")) +
        ### points
        geom_point(aes_string(colour = "Group", fill = "Group",
                       size = "size"),
                   alpha = alpha, shape = shape) +
        ### lines
        stat_function(aes(colour = "Reference"),
                      fun = function(xx) regFce_nonUDIF(xx, 0, x$coef[i, 1], x$coef[i, 2], x$coef[i, 3], x$coef[i, 4],
                                                        x$coef[i, 5]), size = size, linetype = linetype[1], geom = "line") +
        stat_function(aes(colour = "Focal"),
                      fun = function(xx) regFce_nonUDIF(xx, 1, x$coef[i, 1], x$coef[i, 2], x$coef[i, 3], x$coef[i, 4],
                                                        x$coef[i, 5]), size = size, linetype = linetype[2], geom = "line") +
        ### style
        scale_size_continuous(name = "Counts")  +
        scale_colour_manual(breaks = hv$Group, values = col) +
        scale_fill_manual(values = col) +
        scale_linetype_manual(values = linetype,
                              guide = FALSE) +
        ### theme
        ggtitle(TITLE) +
        labs(x = "Standardized total score", y = "Probability of correct answer") +
        ylim(c(0, 1)) + theme_bw() +
        theme(text = element_text(size = 18),
              plot.title = element_text(size = 18, face = "bold", vjust = 1.5),
              axis.line  = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
        ### legend
        theme(legend.box.just = "left",
              legend.justification = c(1, 0),
              legend.position = c(1, 0),
              legend.margin = unit(0, "lines"),
              legend.box = "vertical",
              legend.key.size = unit(1, "lines"),
              legend.text.align = 0,
              legend.title.align = 0)
      print(plot_CC)
    }
  }
  ### checking input
  if (!(plot.type %in% c("cc", "stat"))){
    stop("Possible values of 'plot.type' is 'cc' or 'stat' ",
         call. = FALSE)
  } else {
    if (plot.type == "cc"){
      plotCC(x, item = item,
             col = col, alpha = alpha, shape = shape, size = size,
             linetype = linetype, title = title)
    } else {
      plotstat(x, size = size, title = title)
    }
  }
}

###############################################################################
### FUNCTION fitted.difNLR() ##################################################
###      13.04.2016          ##################################################
###############################################################################
### generic function which extracts fitted values from object difNLR        ###
### INPUT: 'x'     an object of 'difNLR'                                    ###
###        'item'  character string or an object of class numeric or        ###
###                integer specifying for which items characteristic curves ###
###                should be plotted                                        ###
###                Possible values "all" (default), numeric vector of items ###
###                or integer representing item                             ###
### OUTPUT: list of fitted values for all or specified items                ###
###############################################################################
###############################################################################
###############################################################################

fitted.difNLR <- function(object, item = "all", ...){

  ### checking input
  m <- nrow(object$coef)
  if (class(item) == "character"){
    if (item != "all")
      stop("'item' must be either numeric vector or character string 'all' ",
           call. = FALSE)
  } else {
    if (class(item) != "integer" & class(item) != "numeric")
      stop("'item' must be either numeric vector or character string 'all' ",
           call. = FALSE)
  }
  if (class(item) == "numeric" & !all(item %in% 1:m))
    stop("invalid number of 'item'",
         call. = FALSE)
  if (class(item) == "integer" & !all(item %in% 1:m))
    stop("'item' must be either numeric vector or character string 'all' ",
         call. = FALSE)
  if (class(item) == "character"){
    items <- 1:m
  } else {
    items <- item
  }
  if (any(object$conv_fail_which %in% items)){
    if (length(setdiff(items, object$conv_fail_which)) == 0){
      stop(paste("Item ", intersect(object$conv_fail_which, items), " does not converge. No values displayed",
                 sep = "", collapse = "\n"),
           call. = FALSE)
    } else {
      warning(paste("Item ", intersect(object$conv_fail_which, items), " does not converge. No values displayed",
                    sep = "", collapse = "\n"),
              call. = FALSE)
      items <- setdiff(items, object$conv_fail_which)
    }
  }


  ### functions
  NLR <- function(STS, group, a, b, c, aDif, bDif){
    return(c + (1 - c)/(1 + exp(-(a + aDif*group)*(STS - (b + bDif*group)))))
  }

  ### data
  STS <- c(scale(apply(object$data, 1, sum)))

  ### fitted values
  FV <- lapply(items, function(i) NLR(STS, object$group, object$coef[i, 1], object$coef[i, 2],
                                      object$coef[i, 3], object$coef[i, 4], object$coef[i, 5]))
  FV <- lapply(FV, setNames, NULL)
  names(FV) <- paste("Item", items)
  return(FV)

}


###############################################################################
### FUNCTION predict.difNLR() #################################################
###      18.04.2016           #################################################
###############################################################################
###  generic function for predictions from the results of object difNLR     ###
### INPUT: 'x'     an object of 'difNLR'                                    ###
###        'item'  character string or an object of class numeric or        ###
###                integer specifying for which items characteristic curves ###
###                should be plotted                                        ###
###                Possible values "all" (default), numeric vector of items ###
###                or integer representing item                             ###
### OUTPUT: list of predicted values for all or specified items             ###
###         generally for all fitted data or for examinee with specified    ###
###         standard total score and group                                  ###
###############################################################################
###############################################################################
###############################################################################


predict.difNLR <- function(object, item = "all",
                           score, group, ...){

  ### checking input
  m <- nrow(object$coef)
  if (class(item) == "character"){
    if (item != "all")
      stop("'item' must be either numeric vector or character string 'all' ",
           call. = FALSE)
  } else {
    if (class(item) != "integer" & class(item) != "numeric")
      stop("'item' must be either numeric vector or character string 'all' ",
           call. = FALSE)
  }
  if (class(item) == "numeric" & !all(item %in% 1:m))
    stop("invalid number of 'item'",
         call. = FALSE)
  if (class(item) == "integer" & !all(item %in% 1:m))
    stop("'item' must be either numeric vector or character string 'all' ",
         call. = FALSE)
  if (missing(score)){
    score <- c(scale(apply(object$data, 1, sum)))
  }
  if (missing(group)){
    group <- object$group
  }
  if(length(score) != length(group)){
    stop("'score' and 'group' must be the same length",
         call. = FALSE)
  }
  if (class(item) == "character"){
    items <- 1:m
  } else {
    items <- item
  }
  if (any(object$conv_fail_which %in% items)){
    if (length(setdiff(items, object$conv_fail_which)) == 0){
      stop(paste("Item ", intersect(object$conv_fail_which, items), " does not converge. No values displayed",
                 sep = "", collapse = "\n"),
           call. = FALSE)
    } else {
      warning(paste("Item ", intersect(object$conv_fail_which, items), " does not converge. No values displayed",
                    sep = "", collapse = "\n"),
              call. = FALSE)
      items <- setdiff(items, object$conv_fail_which)
    }
  }
  ### functions
  NLR <- function(STS, group, a, b, c, aDif, bDif){
    return(c + (1 - c)/(1 + exp(-(a + aDif*group)*(STS - (b + bDif*group)))))
  }

  ### data
  STS <- score


  ### predicted values
  PV <- lapply(items, function(i) NLR(STS, group, object$coef[i, 1], object$coef[i, 2],
                                      object$coef[i, 3], object$coef[i, 4], object$coef[i, 5]))
  PV <- lapply(PV, setNames, NULL)
  names(PV) <- paste("Item", items)
  return(predict = PV)
}


