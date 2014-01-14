WeibullReg <- function (formula, data = parent.frame(), conf.level = 0.95){

    alpha <- 1 - conf.level
    m <- survival::survreg(formula, data, dist = "weibull")
    mle <- ConvertWeibull(m, alpha)
    return(list(formula = formula, coef = mle$vars, HR = mle$HR, ETR = mle$ETR, summary = summary(m)))
}
