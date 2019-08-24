#' @name iwls.bnreg
#' @title Iterative Weighted Least Squares for logistic models
#'
#' @description Performs logisitc regressions for bionomial data.
#'
#' @usage     iwls.bnreg (y, X, n = NULL, intercept = FALSE, tol = 1e-07, max.iter = 1000)
#'
#' @param y response vector of observed proportions.
#' @param X model matrix of \code{p} covariates.
#' @param n the vector of binomial denominators. When NULL n is equal to a 1s vector.
#' @param intercept a logic value to indicate whether the intercept will be added directly by the function to the model matrix or not. The default value is FALSE i.e. add a vector of 1s to X.
#' @param tol the convergence tolerance criterion.
#' @param max.iter Maximum of iteration to be the limit if convergence is not attained. The default valueis 1000 iterations.
#'
#' @details  The function performs the IWLS algorithm applied to binomial logistic regression (Fox, 2002)
#'
#' @return Returns a list including the following:
#' \item{estimates}{the maximum likelihood estimates of the coefficients.}
#' \item{var}{the covariance matrix of coeeficients.}
#' \item{n.iter}{the number of iterations at convergence.}
#'
#' @export
#'
#' @references Fox John (2002). An R and S-Plus companion to applied regression. Sage Publications.
#'
#' @author Sewanou Honfo \url{<honfosewanou@gmail.com>} and Romain Glèlè Kakaï \url{<glele.romain@gmail.com>}/ LABEF_07_2019
#'
#' @examples
#'
#' # Example on Mroz data from the package car
#' ## Load packages and data
#' library(car)
#' data(Mroz)
#' str(Mroz)
#'
#' ## Convert non numeric variables to numeric
#' Mroz$lfp <- recode(Mroz$lfp, "'yes' = 1; 'no'= 0", as.factor = F)
#' Mroz$wc <- recode(Mroz$wc, "'yes' = 1; 'no'= 0", as.factor = F)
#' Mroz$hc <- recode(Mroz$hc, "'yes' = 1; 'no'= 0", as.factor = F)
#'
#' ## Run the binary logistic regression using the iterative weighted least squares estimation methods
#' attach(Mroz)
#' m.log <- iwls.bnreg(y = lfp, X = cbind(k5, k618, age, wc, hc, lwg, inc), n = NULL, intercept = FALSE, tol = 1e-07, max.iter = 1000)
#'
#' ## View the result
#' m.log$estimates
#' m.log$var
#' m.log$n.iter
#' update

iwls.bnreg <- function (y, X, n = NULL, intercept = FALSE, tol = 1e-07, max.iter = 1000)
{

  X <- as.matrix(X)
  y <- as.matrix(y)
  xnames <- colnames(X)
  if (is.null(xnames)) {
    if (ncol(x) == 1L)
      xnames <- "X"
    else xnames <- paste0("X", 1L:ncol(x))
  }


  if(intercept != FALSE) X <- X

  else {
    X <- cbind(rep(1, length(y)), X)
    xnames <- c("Intercept", xnames)
  }

  yname <- colnames(y)
  if (is.null(xnames)) yname <- "Y"

  if(is.null(n)) n <- rep(1, length(y))

  if(is.vector(n) == 1) n <- n

  if(class(n) != "numeric") stop("The weight vector is not numeric.")

  if(!is.vector(n)== 1) stop("The weight vector is not a vector or is not numeric.")


  complet <- complete.cases(X, y, n)
  if (any(!complet)) {
    warning(sprintf(ngettext(sum(!complet), "%d missing value deleted",
                             "%d missing values deleted"), sum(!complet)),
            domain = NA)
    X <- as.matrix(X)[complet, , drop = FALSE]
    y <- as.matrix(y)[complet, , drop = FALSE]
    n <- n[complet]
  }
  nrx <- dim(X)[1]
  ncx <- dim(X)[2]
  nry <- dim(y)[1]
  ncy <- dim(y)[2]
  ns <- length(n)
  if (nry != nrx)
    stop(sprintf(paste0(ngettext(nrx, "'X' matrix has %d case (row)",
                                 "'X' matrix has %d cases (rows)"), ", ",
                        ngettext(nry, "'Y' has %d case (row)", "'Y' has %d cases (rows)")),
                 nrx, nry), domain = NA)
  if (nry < ncx)
    stop(sprintf(paste0(ngettext(nry, "only %d case", "only %d cases"), ", ", ngettext(ncx,
                                                                                       "but %d variable", "but %d variables")), nry, ncx), domain = NA)
  if (!is.null(n)) {
    if (any(n <= 0))
      stop("Elements of the denominator vector must be positve and not equal to zero")
    if (ns != nry)
      stop(gettextf("number of denominator values = %d should equal %d (number of observations)",
                    ns, nry), domain = NA)
  }


  beta_0 <- beta_t1 <- rep(0, dim(X)[2])

  iter <- 1

  while(iter < max.iter) {

    eta <- X %*% beta_0

    mu <- apply(eta, 1, FUN = function(x){solve(1 + exp(-x))})

    nu <- mu * (1 - mu)

    z <- eta + ((y - mu) / nu)

    w <- n * nu

    mod <- lsfit(x = X, y = z, wt = c(w), intercept = FALSE, tolerance = 1e-07,
                 yname = NULL)

    beta_t1 <- mod$coefficients

    if(max(beta_t1 - beta_0) < tol)
      break

    beta_0 <- beta_t1

    iter <- iter + 1

  }


  V_beta <- solve(t(X) %*% diag(w) %*% X)

  return(list(estimates = beta_t1, var = V_beta, n.iter = iter))


}
