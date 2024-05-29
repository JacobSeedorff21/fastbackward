#' Performs Fast Backward Elimination by AIC
#' @description Performs backward elimination by AIC, backward elimination is
#' performed with a bounding algorithm to make it faster.
#' @param object an object representing a model of an appropriate class. This is
#'  used as the initial model in the stepwise search.
#' @param scope defines the range of models examined in the stepwise search. This
#' should be missing or a single formula. If a formula is included, all of the
#' components on the right-hand-side of the formula are always included in the model.
#' If missing, then only the intercept (if included) is always included in the model.
#' @param scale used in the definition of the AIC statistic for selecting the models,
#' currently only for [lm] and [aov] models (see [extractAIC] for details).
#' @param trace if positive, information is printed during the running of `fastbackward`.
#' Larger values may give more detailed information. If `trace` is greater than 1, then
#' information about which variables at each step are not considered for removal due
#' to the bounding algorithm are printed.
#' @param keep a filter function whose input is a fitted model object and the associated `AIC` statistic,
#' and whose output is arbitrary. Typically `keep` will select a subset of the components
#' of the object and return them. The default is not to keep anything.
#' @param steps the maximum number of steps to be considered. The default is 1000
#' (essentially as many as required). It is typically used to stop the process early.
#' @param k the multiple of the number of degrees of freedom used for the penalty.
#' Only `k = 2` gives the genuine AIC: `k = log(n)` is sometimes referred to as BIC or SBC.
#' @param ... any additional arguments to [extractAIC].
#' @seealso [MASS::stepAIC], [MASS::dropterm], and [extractAIC]
#' @return The stepwise-selected model is returned, with up to two additional components.
#' There is an "`anova`" component corresponding to the steps taken in the search,
#' as well as a "`keep`" component if the `keep=` argument was supplied in the call.
#' The "`Resid. Dev`" column of the analysis of deviance table refers to a constant
#' minus twice the maximized log likelihood: it will be a deviance only in cases
#' where a saturated model is well-defined (thus excluding `lm`, `aov` and `survreg` fits,
#' for example)
#' @details
#'
#' The bounding algorithm allows us to avoid fitting models that cannot possibly
#' provide an improvement in AIC. At a high-level, the algorithm basically
#' works by identifying important predictors whose removal from the current model
#' cannot possibly improve upon the current AIC.
#'
#' Test statistics, p-values, and confidence intervals from the final selected
#' model are not reliable due to the selection process. Thus, it is not recommended
#' to use these quantities.
#'
#' See more details at [MASS::stepAIC].
#'
#' @examples
#' # Loading fastbackward
#'library(fastbackward)
#'
#'# Using examples provided in MASS::stepAIC, but with fastbackward instead
#'## aov with quine dataset
#'quine.hi <- aov(log(Days + 2.5) ~ .^4, MASS::quine)
#'quine.nxt <- update(quine.hi, . ~ . - Eth:Sex:Age:Lrn)
#'quine.stp <- fastbackward(quine.nxt, trace = FALSE)
#'quine.stp$anova
#'
#'## lm with cpus dataset
#'cpus1 <- MASS::cpus
#'for(v in names(MASS::cpus)[2:7])
#'  cpus1[[v]] <- cut(MASS::cpus[[v]], unique(quantile(MASS::cpus[[v]])),
#'                    include.lowest = TRUE)
#'cpus0 <- cpus1[, 2:8]  # excludes names, authors' predictions
#'cpus.samp <- sample(1:209, 100)
#'cpus.lm <- lm(log10(perf) ~ ., data = cpus1[cpus.samp,2:8])
#'cpus.lm2 <- fastbackward(cpus.lm, trace = FALSE)
#'cpus.lm2$anova
#'
#'## glm with bwt dataset
#'example(birthwt, package = "MASS")
#'birthwt.glm <- glm(low ~ ., family = binomial, data = bwt)
#'birthwt.step <- fastbackward(birthwt.glm, trace = FALSE)
#'birthwt.step$anova
#'
#' @export

fastbackward <- function (object, scope, scale = 0, trace = 1, keep = NULL,
                  steps = 1000, k = 2, ...) {
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev))
      dev
    else extractAIC(x, k = 0)[2L]
  }
  cut.string <- function(string) {
    if (length(string) > 1L)
      string[-1L] <- paste("\n", string[-1L], sep = "")
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr,
                                                           namc))
  }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, abs(diff(rdf)))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
                 "\nInitial Model:", deparse(formula(object)), "\nFinal Model:",
                 deparse(formula(fit)), "\n")
    aod <- if (usingCp)
      data.frame(Step = change, Df = ddf, Deviance = dd,
                 `Resid. Df` = rdf, `Resid. Dev` = rd, Cp = AIC,
                 check.names = FALSE)
    else data.frame(Step = change, Df = ddf, Deviance = dd,
                    `Resid. Df` = rdf, `Resid. Dev` = rd, AIC = AIC,
                    check.names = FALSE)
    attr(aod, "heading") <- heading
    class(aod) <- c("Anova", "data.frame")
    fit$anova <- aod
    fit
  }
  Terms <- terms(object)
  object$formula <- Terms
  if (inherits(object, "lme"))
    object$call$fixed <- Terms
  else if (inherits(object, "gls"))
    object$call$model <- Terms
  else object$call$formula <- Terms
  if (missing(scope)) {
    fdrop <- numeric()
  }else {
    fdrop <- attr(terms(update.formula(object, scope)), "factors")
  }
  models <- vector("list", steps)
  if (!is.null(keep))
    keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  bAIC <- extractAIC(fit, scale, k = k, ...)
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]
  if (is.na(bAIC))
    stop("AIC is not defined for this model, so 'fastbackward' cannot proceed")
  if (bAIC == -Inf)
    stop("AIC is -infinity for this model, so 'fastbackward' cannot proceed")
  nm <- 1
  if (trace) {
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(formula(fit))),
        "\n\n", sep = "")
    flush.console()
  }
  models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n -
                         edf, change = "", AIC = bAIC)
  if (!is.null(keep))
    keep.list[[nm]] <- keep(fit, bAIC)
  usingCp <- FALSE
  AICs <- NULL
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    if (!is.null(sp <- attr(Terms, "specials")) && !is.null(st <- sp$strata))
      ffac <- ffac[-st, ]
    scope <- factor.scope(ffac, list(drop = fdrop))
    if(!is.null(AICs)){
      df <- aod$Df[which.min(aod[, nc])]
      AICs <- AICs - df * k
      if(trace > 1){
        if(any(AICs[scope$drop] > AIC + 1e-6)){
          skipind <- scope$drop[AICs[scope$drop] > AIC + 1e-6]
          for(x in skipind){
            cat(paste0("Not considering ", x, " for removal because LB(", format(round(AICs[x], 2)),
                       ") > Best(", format(round(AIC, 2)), ")\n"))
            flush.console()
          }
        }
      }
      scope$add <- c(scope$add, names(which(AICs[scope$drop] > AIC + 1e-6)))
      scope$drop <- names(which(AICs[scope$drop] <= AIC + 1e-6))
    }else{
      AICs <- rep(-Inf, ncol(attr(Terms, "factors")))
      names(AICs) <- colnames(attr(Terms, "factors"))
    }
    aod <- NULL
    change <- NULL
    if(length(scope$drop)) {
      aod <- fastdropterm(object = fit, scope = scope$drop, scale = scale, trace = max(0, trace - 1),
                   k = k, LBs = AICs, AIC = AIC, ...)
      rn <- row.names(aod)
      AICsind <- match(names(AICs), rn)
      AICsind <- AICsind[!is.na(AICsind)]
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L]))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        nc <- match(c("Cp", "AIC"), names(aod))
        nc <- nc[!is.na(nc)][1L]
        ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
        if (any(is.finite(ch) & ch)) {
          warning("0 df terms are changing AIC")
          zdf <- zdf[!ch]
        }
        if (length(zdf) > 0L)
          change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df))
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0)
        break
      nc <- match(c("Cp", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      AICs[rn[AICsind]] <- aod[AICsind, nc]
      o <- order(aod[, nc])
      if (trace)
        print(aod[o, ])
        flush.console()
      if (o[1L] == 1)
        break
      change <- rownames(aod)[o[1L]]
    }
    usingCp <- match("Cp", names(aod), 0) > 0
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)
    nnew <- nobs(fit, use.fallback = TRUE)
    if (all(is.finite(c(n, nnew))) && nnew != n)
      stop("number of rows in use has changed: remove missing values?")
    Terms <- terms(fit)
    bAIC <- extractAIC(fit, scale, k = k, ...)
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    if (trace) {
      cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n",
          cut.string(deparse(formula(fit))), "\n\n", sep = "")
      flush.console()
    }
    if (bAIC >= AIC + 1e-07)
      break
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n -
                           edf, change = change, AIC = bAIC)
    if (!is.null(keep))
      keep.list[[nm]] <- keep(fit, bAIC)
  }
  if (!is.null(keep))
    fit$keep <- re.arrange(keep.list[seq(nm)])
  step.results(models = models[seq(nm)], fit, object, usingCp)
}

#' Helper Function for fastbackward Function
#' @description Performs a fast backward elimination step, this is mainly used
#' by the fastbackward function.
#' @param object a fitted model object.
#' @param scope a character vector of names giving the terms to be considered for dropping.
#' @param LBs lower bounds of the AIC for the removal of each variable.
#' @param AIC the best AIC value observed thus far, used in combination with LBs
#' to avoid fitting unnecessary models.
#' @param trace if trace is greater than 0, then information about which variables
#' are not considered for removal due to the bounding algorithm are printed.
#' @param ... any additional arguments to [MASS::dropterm].
#' @seealso [fastbackward], [MASS::dropterm], and [extractAIC]
#' @return An object of class "anova" summarizing the differences in fit between
#' the necessary models.
#' @noRd

fastdropterm <- function(object, scope, LBs, AIC, trace = 0, ...){
  fitorder <- scope[order(LBs[scope])]
  aod <- NULL
  for(i in fitorder){
    if(LBs[i] <= AIC + 1e-6){
      tempscope <- i
      aod1 <- MASS::dropterm(object, tempscope, trace = trace, ...)
      nc <- match(c("Cp", "AIC"), names(aod1))
      nc <- nc[!is.na(nc)][1L]
      tempAIC <- min(aod1[, nc], na.rm = TRUE)
      if(any(!is.na(aod1[, nc])) && tempAIC < AIC){
        AIC <- tempAIC
        if(trace > 0){
          cat(paste0("Updating the best observed AIC to be ", format(round(AIC, 2)), "\n"))
          flush.console()
        }
      }
      if(is.null(aod)){
        aod <- aod1
      }else{
        aod <- rbind(aod, aod1[-1, ])
      }
    }else if(trace > 0){
      cat(paste0("Not considering ", i, " for removal because LB(", format(round(LBs[i], 2)),
                 ") > Best(", format(round(AIC, 2)), ")\n"))
      flush.console()
    }
  }
  return(aod)
}
