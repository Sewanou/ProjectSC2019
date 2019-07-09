#' @name gen.data
#' @title Function for data generation
#'
#' @description generates data from a one-way random effects ANOVA model.
#'
#' @usage     gen.data (ng, ns, moy, varg, vars)
#'
#' @param ng the number of groups.
#' @param ns the number of subject in each group. All groups have the same number of subjects.
#' @param moy the overall mean.
#' @param varg the variance of group means.
#' @param vars the variance of subjects within groups.
#'#'
#' @return Returns a data frame with ng * ns rows, and two columns called "group" and "value", giving the group id for each subject (1 to ng), and the value for that subject.
#'
#' @export
#'
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>} and Romain Glèlè Kakaï \url{<glele.romain@gmail.com>} / LABEF_07_2019

gen.data <- function(ng, ns, moy, varg, vars)
{
  if(ng < 1 || ns < 1 || varg < 0 || vars < 0)
    stop("Bad argument for gen.data")

  g <- c()
  v <- c()

  for(i in 1 : ng){

    g <- c(g, rep(i, ns))

    v <- c(v, moy + rnorm(1, 0, sqrt(varg)) + rnorm(ns, 0, sqrt(vars)))

  }

  k <- data.frame(group = g, value = v)

  return(k)

}

#' @name def.t
#' @title Standard t test
#'
#' @description Performs the Standard t-test of mean = 0 for a one-way random effects model.
#'
#' @usage     def.t (X)
#'
#' @param X a \code{gen.data} object.
#'
#' @return Returns a list including:
#' \item{p.value}{the probability for the Standard two-sided t-test of the overall mean being zero.}
#' \item{moy.group}{the vector of means of values within each group.}
#'
#' @export
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>} and Romain Glèlè Kakaï \url{<glele.romain@gmail.com>}/ LABEF_07_2019

def.t <- function(X)
{
  g.vec <- unique(X$group)
  ng <- length(g.vec)
  ns <- length(X$group) / ng

  if(ng < 2) stop("Cannot do a t-test with less than two groups")

  k1 <- numeric(ng)

  for(i in g.vec){

    if(sum(X$group == i) != ns) stop("The number of subjects is not the same for all groups")

    k1[i] <- mean(X$value[X$group == i])

  }

  t.statistic <- mean(k1) / sqrt(var(k1) / ng)
  p <- 2 * pt(-abs(t.statistic), ng - 1)

  return(  list(p.value = p, moy.group = k1))

}

#' @name mod.t
#' @title Modified t test
#'
#' @description Performs the Modified t-test of mean = 0 for a one-way random effects model.
#'
#' @usage     mod.t (X)
#'
#' @param X a \code{gen.data} object.
#'
#' @return Returns a list including:
#' \item{p.value}{the probability for the Modified two-sided t-test of the overall mean being zero.}
#' \item{moy.group}{the vector of means of values within each group.}
#' \item{vars.group}{the vector of variances of values within each group}
#'
#' @export
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>} and Romain Glèlè Kakaï \url{<glele.romain@gmail.com>}/ LABEF_07_2019

mod.t <- function(X)
{
  g.vec <- unique(X$group)
  ng <- length(g.vec)
  ns <- length(X$group) / ng

  if(ng < 2) stop("Cannot perform a t-test with less than two groups")

  if(ns < 2) stop("Cannot perform the modified t-test with less than two subjects per group")

  k2 <- numeric(ng)
  k3 <- numeric(ns)

  for(i in g.vec){

    if(sum(X$group == i) != ns) stop("The number of subjects is not the same for all groups")

    k2[i] <- mean(X$value[X$group == i])
    k3[i] <- var(X$value[X$group == i])

  }

  nu <- max(var(k2), mean(k3) / ns)

  t.statistic <- mean(k2) / sqrt(nu / ng)
  p <- 2 * pt(-abs(t.statistic), ng - 1)

  return(list(p.value = p, moy.group = k2, vars.group = k3))

}


#' @name sim.test
#' @title A test on simulated data
#'
#' @description applies a test to simulated one-way ANOVA data sets.
#'
#' @usage     sim.test (r, ng, ns, moy, varg, vars, method = "def.t")
#'
#' @param r the number of simulated data sets.
#' @param ng the number of groups in a data set.
#' @param ns the number of subjects per group.
#' @param moy the overall mean.
#' @param varg the variance of group means.
#' @param vars the variance within a group.
#' @param method a character value specifying which t-test to apply. "def.t"  (by default) indicates the standard t-test and "mod.t", the modified t-test.
#'
#' @return Returns a vector of p-values from the selected test.
#'
#' @export
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>} and Romain Glèlè Kakaï \url{<glele.romain@gmail.com>}/ LABEF_07_2019

sim.test <- function(r, ng, ns, moy, varg, vars, method = "def.t")
{
  if(r < 1 || ng < 2 || ns < 1 || varg < 0 || vars < 0) stop("Bad argument for sim. test")

  p.values <- numeric(r)

  for(i in 1 : r){

    q <- gen.data(ng, ns, moy, varg, vars)

    if(method == "def.t")
       p.values[i] <- def.t(q)$p.value

    if(method == "mod.t")
      p.values[i] <- mod.t(q)$p.value

  }

  return(p.values)

}

#' @name comp.test
#' @title Rejection probability test
#'
#' @description estimates the rejection probability for a test at a given level.
#'
#' @details The test is performed for various values of the overall mean, and the rejection probability estimated for each such value.
#'
#' @usage     comp.test (r, level, ng, ns, moys, varg, vars, method = "def.t")
#'
#' @param r the number of simulated data sets to use in finding each estimate.
#' @param level the significance level.
#' @param ng the number of groups in a data set.
#' @param ns the number of subjects per group.
#' @param moys the vector of means for the data sets.
#' @param varg the variance of group means.
#' @param vars the variance within a group.
#' @param method a character value specifying which t-test to apply. "def.t" (by default) indicates the standard t-test and "mod.t", the modified t-test.
#'
#' @return Returns a vector of estimated rejection probabilities for each mean.
#'
#' @export
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>} and Romain Glèlè Kakaï \url{<glele.romain@gmail.com>}/ LABEF_07_2019

comp.test <- function(r, level, ng, ns, moys, varg, vars, method)
{
  if(r < 1 || level <= 0 || level > 1 || ng < 2 || ns < 1 || varg < 0 || vars < 0) stop("Bad argument for comp. test")

  prob <- numeric(length(moys))

  for(i in 1 : length(moys)){

    p <- sim.test (r, ng, ns, moys[i], varg, vars, method = "def.t")
    prob[i] <- mean(p <= level)

  }

  return(prob)

}

