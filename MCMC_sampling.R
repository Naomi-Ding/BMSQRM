###### Step 3: Estimation via MCMC given {h^*_{k,i}(t), i=1,...,n, k=1,2} ######
require(GIGrvg)
# require(nimble)
require(MfUSampler)
source("utilities.R")

##### Likelihood #####
require(Rcpp)
sourceCpp("logLikelihood.cpp")

#### Priors ######
logprior.b <- function(b, a.DP) {
    sum(log(1 - b)) * (a.DP - 1)
}

logprior.mu <- function(mu, mean.mu = 0, sd.mu = 1) {
    -sum((mu - mean.mu)^2 / (2 * sd.mu^2))
} # {-sum(mu^2) / 2}

logprior.xi <- function(xi, shape.xi = 2, scale.xi = 2) {
    -sum((shape.xi + 1) * log(xi) + scale.xi / xi)
} # {-sum(3 * log(xi) + 1 / xi)}


##### Posteriors ######
logpost.theta <- function(theta, beta.t, sigma.t, e.t, link, fz.idx = 1,
                          z_mu = 0, z_sigma = 0.7,
                          pi.t, mu.t, xi.t, y, preds.all, k1, k2, r = 1) {
    logLikelihood_cpp(
        alpha.transform(theta, r), beta.t, sigma.t, e.t,
        link, fz.idx, z_mu, z_sigma, pi.t, mu.t, xi.t,
        y, preds.all, k1, k2
    )
}

logpost.b <- function(b, alpha.t, beta.t, sigma.t, e.t, mu.t, xi.t, y, preds.all, k1, k2, a.DP) {
    res <- logLikelihood_cpp(
        alpha.t, beta.t, sigma.t, e.t, FALSE, 2, 0, 0,
        pi.transform(b), mu.t, xi.t, y, preds.all, k1, k2
    ) + logprior.b(b, a.DP)
    return(res)
}

logpost.mu <- function(mu, alpha.t, beta.t, sigma.t, e.t, pi.t, xi.t, y, preds.all,
                       k1, k2, mean.mu = 0, sd.mu = 1) {
    res <- logLikelihood_cpp(
        alpha.t, beta.t, sigma.t, e.t, FALSE, 2, 0, 0,
        pi.t, mu, xi.t, y, preds.all, k1, k2
    ) + logprior.mu(mu, mean.mu, sd.mu)
    return(res)
}

logpost.xi <- function(xi, alpha.t, beta.t, sigma.t, e.t, pi.t, mu.t, y, preds.all,
                       k1, k2, shape.xi = 2, scale.xi = 2) {
    res <- logLikelihood_cpp(
        alpha.t, beta.t, sigma.t, e.t, FALSE, 2, 0, 0,
        pi.t, mu.t, xi, y, preds.all, k1, k2
    ) + logprior.xi(xi, shape.xi, scale.xi)
    return(res)
}

logpost.beta.kl <- function(beta.kl, idx, beta.t, alpha.t, sigma.t, e.t, link,
                            fz.idx = 2, z_mu = 0, z_sigma = 0.7,
                            pi.t, mu.t, xi.t, y, preds.all, k1, k2, eta.t.kl) {
    beta.t[idx] <- beta.kl
    logp.beta.kl <- logLikelihood_cpp(
        alpha.t, beta.t, sigma.t, e.t, link, fz.idx, z_mu, z_sigma,
        pi.t, mu.t, xi.t, y, preds.all, k1, k2
    ) - beta.kl^2 / (2 * eta.t.kl)
    return(logp.beta.kl)
}


##### MCMC ######
sample.theta <- function(theta.t, beta.t, sigma.t, e.t, link = TRUE, fz.idx = 2,
                         z_mu = 0, z_sigma = 0.7,
                         pi.t, mu.t, xi.t, y, preds.all, k1, k2,
                         lb.theta = NULL, ub.theta = NULL,
                         slice.w = 0.5, slice.m, r = 1) {
    p <- length(theta.t) + 1
    if (is.null(lb.theta)) {
        lb.theta <- c(rep(-pi / 2, p - 2), 0)
    }
    if (is.null(ub.theta)) {
        ub.theta <- c(rep(pi / 2, p - 2), 2 * pi)
    }
    theta.t <- MfU.Sample(theta.t, logpost.theta,
        uni.sampler = "slice",
        beta.t = beta.t, sigma.t = sigma.t, e.t = e.t, link = link, fz.idx = fz.idx,
        z_mu = z_mu, z_sigma = z_sigma, pi.t = pi.t,
        mu.t = mu.t, xi.t = xi.t, y = y, preds.all = preds.all,
        k1 = k1, k2 = k2, r = r,
        control = MfU.Control( # theta_j \in (-pi/2, pi/2), j=1,...,p-2
            n = p - 1,
            slice.lower = lb.theta, # theta_{p-1} \in (0, 2*pi)
            slice.upper = ub.theta,
            slice.w = slice.w, slice.m = slice.m
        )
    )
    return(theta.t)
}

sample.b <- function(b.t, alpha.t, beta.t, sigma.t, e.t, mu.t, xi.t,
                     y, preds.all, k1, k2, a.DP,
                     lb.b = NULL, ub.b = NULL, slice.w = 0.1, slice.m) {
    K <- length(b.t) + 1
    if (is.null(lb.b)) {
        lb.b <- 0
    }
    if (is.null(ub.b)) {
        ub.b <- 1
    }
    b.t <- MfU.Sample(b.t, logpost.b,
        uni.sampler = "slice",
        alpha.t = alpha.t, beta.t = beta.t, sigma.t = sigma.t, e.t = e.t,
        mu.t = mu.t, xi.t = xi.t, y = y, preds.all = preds.all, k1 = k1, k2 = k2,
        a.DP = a.DP,
        control = MfU.Control(
            n = K - 1, slice.lower = lb.b, slice.upper = ub.b,
            slice.w = slice.w, slice.m = slice.m
        )
    )
    return(b.t)
}

sample.mu <- function(mu.t, alpha.t, beta.t, sigma.t, e.t, pi.t,
                      xi.t, y, preds.all, k1, k2, mean.mu = 0, sd.mu = 1,
                      lb.mu = -Inf, ub.mu = Inf, slice.m) {
    if (is.null(lb.mu)) {
        lb.mu <- -Inf
    }
    if (is.null(ub.mu)) {
        ub.mu <- Inf
    }

    mu.t <- MfU.Sample(mu.t, logpost.mu,
        uni.sampler = "slice",
        alpha.t = alpha.t, beta.t = beta.t, sigma.t = sigma.t, e.t = e.t,
        pi.t = pi.t, xi.t = xi.t, y = y, preds.all = preds.all,
        k1 = k1, k2 = k2, mean.mu = mean.mu, sd.mu = sd.mu,
        control = MfU.Control(n = length(mu.t), slice.m = slice.m)
    )
    return(mu.t)
}

sample.xi <- function(xi.t, alpha.t, beta.t, sigma.t, e.t, pi.t,
                      mu.t, y, preds.all, k1, k2, shape.xi = 2, scale.xi = 2,
                      lb.xi = NULL, ub.xi = NULL, slice.m) {
    K <- length(xi.t)
    if (is.null(lb.xi)) {
        lb.xi <- 0
    }
    if (is.null(ub.xi)) {
        ub.xi <- Inf
    }
    xi.t <- MfU.Sample(xi.t, logpost.xi,
        uni.sampler = "slice",
        alpha.t = alpha.t, beta.t = beta.t, sigma.t = sigma.t, e.t = e.t,
        pi.t = pi.t, mu.t = mu.t, y = y, preds.all = preds.all,
        k1 = k1, k2 = k2, shape.xi = shape.xi, scale.xi = scale.xi,
        control = MfU.Control(
            n = K, slice.lower = lb.xi, slice.upper = ub.xi,
            slice.m = slice.m
        )
    )
    return(xi.t)
}

sample.beta.kl <- function(alpha.t, beta.t, idx, sigma.t, e.t, link = TRUE, fz.idx = 2,
                           z_mu = 0, z_sigma = 0.7, fz = NULL, fz.prime = NULL,
                           pi.t, mu.t, xi.t, eta.t.kl, y, preds.all, k1, k2, p.kl = NULL, # n.f,
                           slice.m) {
    ### Input:
    ### Data:
    ###     y: response, vector of length n
    ###     pred.all: all the predictors, (p + L0, n)
    ###             the first p rows represent scalar predictors, (p, n)
    ###             while the rest rows represent decomposed functional covaraites (L0, n)
    ###     d: number of functional predictors
    ###     L: number of basis functions
    ###     x: scalar covariates, matrix of size (p,n)
    ###     w: reconstructed covariates by stacking L functional covariates,
    ###        matrix of size (d*L, n)
    ### Parameters:
    ###     alpha.t: vector of length p
    ###     beta.t: vector of length 2*L
    ###     e.t: vector of length n, samples of {e_1, ..., e_n}
    ###     pi.t: vector of length K, (pi_1, ..., pi_K)
    ###     mu.t: vector of length K
    ###     xi.t: vector of length K
    ###     eta.t: vector of length 2*L, variance related to componentwise beta.t
    ###     k1: scalar
    ###     k2: scalar
    ###     p.kl: scalar in (0,1), prior for delta_0(beta_{kl}), 1 - 1 / (2^(L-1))
    ###     n.f: # of functional predictors
    ### Output:
    ###     beta.t: updated beta_{k,l}, a vector of length 2 * L
    p <- length(alpha.t)
    x <- preds.all[1:p, ]
    w <- preds.all[-(1:p), ]
    n <- length(y)
    L <- length(beta.t)
    if (is.null(p.kl)) {
        p.kl <- 1 - 1 / (2^(L - 1)) # default value
    }
    if (p.kl == 1) {
        beta.kl <- 0
    } else if (p.kl == 0) {
        # sample from the proportional posterior distribution
        beta.kl <- MfU.Sample(beta.t[idx], logpost.beta.kl,
            uni.sampler = "slice", idx = idx, beta.t = beta.t,
            alpha.t = alpha.t, sigma.t = sigma.t, e.t = e.t,
            link = link, fz.idx = fz.idx, z_mu = z_mu, z_sigma = z_sigma,
            pi.t = pi.t, mu.t = mu.t,
            xi.t = xi.t, y = y, preds.all = preds.all,
            k1 = k1, k2 = k2, eta.t.kl = eta.t.kl,
            control = MfU.Control(n = 1, slice.m = slice.m)
        )
        # }
    } else { # p.kl in (0,1)
        ax <- alpha.t %*% x
        # eta.inv <- 1 / eta.t # vector of length 2 * L
        sigma.w.log <- log(k2) + log(sigma.t) + log(e.t) # (n)
        sigma.w <- exp(sigma.w.log)
        # const.w2.vec.log <- log(1 - p.kl) - log(2 * pi) * n / 2 - sum(sigma.w.log) / 2 + log(eta.inv) / 2

        # for (idx in 1:L) {
        z <- ax + beta.t %*% w # vector of length n
        z.kl <- z - beta.t[idx] * w[idx, ] # vector of length n
        if (link) {
            f.z <- fz(z)
            f.z.kl <- fz(z.kl)
            f.prime.z <- fz.prime(z)
        } else {
            f.z <- numeric(n) # {f(z_i)}, vector of length n
            f.prime.z <- numeric(n)
            f.z.kl <- numeric(n)
            for (i in 1:n) {
                f.z[i] <- sum(pi.t * pnorm(z[i], mu.t, sqrt(xi.t))) # f(z_i)
                f.prime.z[i] <- sum(pi.t * dnorm(z[i], mu.t, sqrt(xi.t)))
                f.z.kl[i] <- sum(pi.t * pnorm(z.kl[i], mu.t, sqrt(xi.t)))
            }
        }
        ## w_1
        w1.log <- log(p.kl) + sum(dnorm(y, f.z.kl + k1 * e.t, sqrt(sigma.w), log = TRUE))
        ## w_2
        x.star <- w[idx, ] * f.prime.z # vector of length n
        y.star <- y - f.z + beta.t[idx] * x.star - k1 * e.t # vector of length n
        sigma2.beta.inv <- 1 / eta.t.kl + sum(x.star^2 / sigma.w)
        log.w2 <- -sum(y.star^2 / sigma.w) / 2 + sum(x.star * y.star / sigma.w)^2 / (2 * sigma2.beta.inv)
        const.w2.vec.log.kl <- log(1 - p.kl) - log(2 * pi) * n / 2 - sum(sigma.w.log) / 2 - log(eta.t.kl) / 2
        const.w2.log <- const.w2.vec.log.kl - log(sigma2.beta.inv) / 2
        w2.log <- const.w2.log + log.w2
        ## \tilde{p_{k,l}} # post.p.kl <- 1 / (1 + w2 / w1)
        w21.log <- w2.log - w1.log
        if (is.infinite(exp(w21.log))) {
            p.kl.tilde <- 0
        } else {
            log.post.p.kl <- -log(1 + exp(w21.log))
            p.kl.tilde <- exp(log.post.p.kl)
        }
        u <- rbinom(1, 1, p.kl.tilde)
        if (u == 1) {
            beta.kl <- 0
        } else { # sample from the proportional posterior distribution
            beta.kl <- MfU.Sample(beta.t[idx], logpost.beta.kl,
                uni.sampler = "slice", idx = idx, beta.t = beta.t,
                alpha.t = alpha.t, sigma.t = sigma.t, e.t = e.t,
                link = link, fz.idx = fz.idx, z_mu = z_mu, z_sigma = z_sigma,
                pi.t = pi.t, mu.t = mu.t,
                xi.t = xi.t, y = y, preds.all = preds.all,
                k1 = k1, k2 = k2, eta.t.kl = eta.t.kl,
                control = MfU.Control(n = 1, slice.m = slice.m)
            )
        }
    }
    return(beta.kl)
}



sample.beta <- function(alpha.t, beta.t, sigma.t, e.t, link = TRUE, fz.idx = 2,
                        z_mu = 0, z_sigma = 0.7, fz = NULL, fz.prime = NULL,
                        pi.t, mu.t, xi.t, eta.t, y, preds.all, k1, k2, p.kl = NULL, # n.f,
                        slice.m) {
    ### Input:
    ### Data:
    ###     y: response, vector of length n
    ###     pred.all: all the predictors, (p + L0, n)
    ###             the first p rows represent scalar predictors, (p, n)
    ###             while the rest rows represent decomposed functional covaraites (L0, n)
    ###     d: number of functional predictors
    ###     L: number of basis functions
    ###     x: scalar covariates, matrix of size (p,n)
    ###     w: reconstructed covariates by stacking L functional covariates,
    ###        matrix of size (d*L, n)
    ### Parameters:
    ###     alpha.t: vector of length p
    ###     beta.t: vector of length 2*L
    ###     e.t: vector of length n, samples of {e_1, ..., e_n}
    ###     pi.t: vector of length K, (pi_1, ..., pi_K)
    ###     mu.t: vector of length K
    ###     xi.t: vector of length K
    ###     eta.t: vector of length 2*L, variance related to componentwise beta.t
    ###     k1: scalar
    ###     k2: scalar
    ###     p.kl: scalar or vector in (0,1), prior for delta_0(beta_{kl}), 1 - 1 / (2^(L-1))
    ###     n.f: # of functional predictors
    ### Output:
    ###     beta.t: updated beta_{k,l}, a vector of length 2 * L
    p <- length(alpha.t)
    if (p < dim(preds.all)[1]) {
        x <- preds.all[1:p, ]
        w <- preds.all[-(1:p), ]
        n <- length(y)
        # L <- round(length(beta.t) / n.f)
        L <- length(beta.t)
        if (is.null(p.kl)) {
            p.kl <- 1 - 1 / (2^(L - 1)) # default value
        }

        if (length(p.kl) == 1) { # uniform p.kl for all beta.kl
            if (p.kl == 1) {
                beta.t <- rep(0, L)
            } else if (p.kl == 0) {
                # sample from the proportional posterior distribution
                # for (k in 1:n.f) {
                #     for (l in 1:L) {
                #         idx <- L * (k - 1) + l # index for beta_{k,l}
                for (idx in 1:L) {
                    beta.kl <- MfU.Sample(beta.t[idx], logpost.beta.kl,
                        uni.sampler = "slice", idx = idx, beta.t = beta.t,
                        alpha.t = alpha.t, sigma.t = sigma.t, e.t = e.t,
                        link = link, fz.idx = fz.idx, z_mu = z_mu, z_sigma = z_sigma,
                        pi.t = pi.t, mu.t = mu.t,
                        xi.t = xi.t, y = y, preds.all = preds.all,
                        k1 = k1, k2 = k2, eta.t.kl = eta.t[idx],
                        control = MfU.Control(n = 1, slice.m = slice.m)
                    )
                    beta.t[idx] <- beta.kl
                }
                # }
            } else {
                ax <- alpha.t %*% x
                eta.inv <- 1 / eta.t # vector of length 2 * L
                sigma.w.log <- log(k2) + log(sigma.t) + log(e.t) # (n)
                sigma.w <- exp(sigma.w.log)
                const.w2.vec.log <- log(1 - p.kl) - log(2 * pi) * n / 2 - sum(sigma.w.log) / 2 + log(eta.inv) / 2

                # for (k in 1:n.f) {
                #     for (l in 1:L) {
                # idx <- L * (k - 1) + l # index for beta_{k,l}
                for (idx in 1:L) {
                    z <- ax + beta.t %*% w # vector of length n
                    z.kl <- z - beta.t[idx] * w[idx, ] # vector of length n
                    if (link) {
                        f.z <- fz(z)
                        f.z.kl <- fz(z.kl)
                        f.prime.z <- fz.prime(z)
                    } else {
                        f.z <- numeric(n) # {f(z_i)}, vector of length n
                        f.prime.z <- numeric(n)
                        f.z.kl <- numeric(n)
                        for (i in 1:n) {
                            f.z[i] <- sum(pi.t * pnorm(z[i], mu.t, sqrt(xi.t))) # f(z_i)
                            f.prime.z[i] <- sum(pi.t * dnorm(z[i], mu.t, sqrt(xi.t)))
                            f.z.kl[i] <- sum(pi.t * pnorm(z.kl[i], mu.t, sqrt(xi.t)))
                        }
                    }
                    ## w_1
                    w1.log <- log(p.kl) + sum(dnorm(y, f.z.kl + k1 * e.t, sqrt(sigma.w), log = TRUE))
                    ## w_2
                    x.star <- w[idx, ] * f.prime.z # vector of length n
                    y.star <- y - f.z + beta.t[idx] * x.star - k1 * e.t # vector of length n
                    sigma2.beta.inv <- eta.inv[idx] + sum(x.star^2 / sigma.w)
                    log.w2 <- -sum(y.star^2 / sigma.w) / 2 + sum(x.star * y.star / sigma.w)^2 / (2 * sigma2.beta.inv)
                    const.w2.log <- const.w2.vec.log[idx] - log(sigma2.beta.inv) / 2
                    w2.log <- const.w2.log + log.w2
                    ## \tilde{p_{k,l}} # post.p.kl <- 1 / (1 + w2 / w1)
                    w21.log <- w2.log - w1.log
                    if (is.infinite(exp(w21.log))) {
                        p.kl.tilde <- 0
                    } else {
                        log.post.p.kl <- -log(1 + exp(w21.log))
                        p.kl.tilde <- exp(log.post.p.kl)
                    }
                    u <- rbinom(1, 1, p.kl.tilde)
                    if (u == 1) {
                        beta.kl <- 0
                    } else { # sample from the proportional posterior distribution
                        beta.kl <- MfU.Sample(beta.t[idx], logpost.beta.kl,
                            uni.sampler = "slice", idx = idx, beta.t = beta.t,
                            alpha.t = alpha.t, sigma.t = sigma.t, e.t = e.t,
                            link = link, fz.idx = fz.idx, z_mu = z_mu, z_sigma = z_sigma,
                            pi.t = pi.t, mu.t = mu.t,
                            xi.t = xi.t, y = y, preds.all = preds.all,
                            k1 = k1, k2 = k2, eta.t.kl = eta.t[idx],
                            control = MfU.Control(n = 1, slice.m = slice.m)
                        )
                    }
                    beta.t[idx] <- beta.kl
                    # }
                }
            }
        } else if (length(p.kl) == L) { # unique p.kl for each beta.kl
            for (idx in 1:L) {
                beta.kl <- sample.beta.kl(
                    alpha.t, beta.t, idx, sigma.t, e.t, link, fz.idx,
                    z_mu, z_sigma, fz, fz.prime, pi.t, mu.t, xi.t,
                    eta.t[idx], y, preds.all, k1, k2, p.kl[idx], slice.m
                )
                beta.t[idx] <- beta.kl
            }
        }
        return(beta.t)
    }
}

sample.paras.single <- function(theta.t, alpha.t, beta.t, sigma.t, link = TRUE, fz.idx = 1, z_mu = 0, z_sigma = 0.7,
                                fz = NULL, fz.prime = NULL, b.t, pi.t, mu.t, xi.t, y, preds.all, k1, k2,
                                a.sigma = 1, b.sigma = 1, a.DP = 2, mean.mu = 0, sd.mu = 1,
                                shape.xi = 2, scale.xi = 2, lam = 1 / 2, p.kl = NULL, n.f,
                                lb.list, ub.list, slice.w, slice.m, r = 1, eta.scale = 2) {
    ### Input:
    ### Data:
    ###     y: response, vector of length n
    ###     pred.all: all the predictors, (p + L0, n)
    ###             the first p rows represent scalar predictors, (p, n)
    ###             while the rest rows represent decomposed functional covaraites (L0, n)
    ###     d: number of functional predictors
    ###     L: number of basis functions
    ###     x: scalar covariates, matrix of size (p,n)
    ###     w: reconstructed covariates by stacking L functional covariates,
    ###        matrix of size (d*L, n)
    ### Parameters:
    ###     alpha.t: vector of length p
    ###     beta.t: vector of length 2*L
    ###     sigma.t: scalar
    ###     link: known or unknown f(z)
    ###     fz: link function, f(z), if link=='known'
    ###     pi.t: vector of length K, (pi_1, ..., pi_K)
    ###     mu.t: vector of length K,
    ###     xi.t: vector of length K
    ###     k1: scalar
    ###     k2: scalar
    ###     a.sigma & b.sigma: parameters in the prior of sigma
    ###     a.DP: parameter in the prior of b
    ###     lam: the order in the posterior of e_i, default=1/2
    ###     p.kl: scalar, default=1 - 1/(2^(L-1))
    ### Output:
    ###     a list containing a new sample of parameters of interest

    n <- length(y)
    # alpha.t <- alpha.transform(theta.t, r = r)
    p <- length(alpha.t) # number of scalar covariates
    # print(alpha.t)
    if (p < dim(preds.all)[1]) {
        # L0 <- dim(preds.all)[1] - p
        x <- preds.all[1:p, ]
        w <- preds.all[-(1:p), ]
        z <- alpha.t %*% x + beta.t %*% w # vector of length n
        # L <- round(length(beta.t) / n.f)
        if (is.null(p.kl)) {
            # p.kl <- 1 - 1 / (2^(L - 1)) # default value
            p.kl <- 1 - 1 / (2^eta.scale + 1)
        }
    } else {
        z <- alpha.t %*% preds.all
        beta.t <- 0
        p.kl <- 0
        eta.t <- 0
    }

    if (link) {
        f.z <- fz(z)
    } else {
        f.z <- numeric(n) # {f(z_i)}, vector of length n
        for (i in 1:n) {
            f.z[i] <- sum(pi.t * pnorm(z[i], mu.t, sqrt(xi.t))) # f(z_i)
        }
    }
    y.diff <- y - f.z # vector of length n

    ## (1) e_i, i=1,...,n
    psi.e <- k1^2 / (k2 * sigma.t) + 2 / sigma.t # scalar, psi_ei, same for each ei
    chi.e <- y.diff^2 / (k2 * sigma.t) # vector of length n
    e.t <- numeric(n)
    for (i in 1:n) {
        # e.t[i] <- rgig(1, param = c(chi.e[i], psi.e, lam))
        e.t[i] <- GIGrvg::rgig(1, lambda = lam, chi = chi.e[i], psi = psi.e)
    }

    ## (2) sigma
    sigma.t <- 1 / (rgamma(1,
        shape = 3 * n / 2 + a.sigma + 1,
        rate = b.sigma + sum(e.t + (y.diff - k1 * e.t)^2 / (2 * k2 * e.t))
    ))

    ### (4) theta & alpha
    lb.theta <- lb.list$theta
    ub.theta <- ub.list$theta
    if (is.list(slice.w)) {
        w.theta <- slice.w$theta
    } else {
        w.theta <- slice.w
    }
    theta.t <- sample.theta(theta.t, beta.t, sigma.t, e.t,
        link, fz.idx, z_mu, z_sigma,
        pi.t, mu.t, xi.t, y, preds.all, k1, k2,
        lb.theta, ub.theta,
        slice.w = w.theta, slice.m = slice.m, r = r
    )
    alpha.t <- alpha.transform(theta.t, r)

    if (!link) {
        ### (5) b
        lb.b <- lb.list$b
        ub.b <- ub.list$b
        if (is.list(slice.w)) {
            w.b <- slice.w$b
        } else {
            w.b <- slice.w
        }
        b.t <- sample.b(b.t, alpha.t, beta.t, sigma.t, e.t,
            mu.t, xi.t, y, preds.all, k1, k2, a.DP,
            lb.b, ub.b,
            slice.w = w.b, slice.m = slice.m
        )
        pi.t <- pi.transform(b.t)

        ### (6) mu
        lb.mu <- lb.list$mu
        ub.mu <- ub.list$mu
        mu.t <- sample.mu(mu.t, alpha.t, beta.t, sigma.t, e.t,
            pi.t, xi.t, y, preds.all, k1, k2, mean.mu, sd.mu,
            lb.mu, ub.mu,
            slice.m = slice.m
        )

        ### (7) xi
        lb.xi <- lb.list$xi
        ub.xi <- ub.list$xi
        xi.t <- sample.xi(xi.t, alpha.t, beta.t, sigma.t, e.t,
            pi.t, mu.t, y, preds.all, k1, k2, shape.xi, scale.xi,
            lb.xi, ub.xi,
            slice.m = slice.m
        )
    }

    if (p < dim(preds.all)[1]) {
        ### (3) eta_{k,l}, l=1,...,L, k=1,2
        eta.t <- numeric(length(beta.t))
        eta.a <- 5 / 2 + 1
        # eta.b <- 2^(L - 2) + beta.t^2 / 2
        eta.b <- eta.scale + beta.t^2 / 2
        for (l in 1:(length(beta.t))) {
            eta.t[l] <- 1 / (rgamma(1, shape = eta.a, rate = eta.b[l]))
        }

        ### (8) beta_{k,l}, l=1,...,L, k=1,...,n.f
        beta.t <- sample.beta(
            alpha.t, beta.t, sigma.t, e.t, link, fz.idx,
            z_mu, z_sigma, fz, fz.prime, pi.t, mu.t, xi.t,
            eta.t, y, preds.all, k1, k2, p.kl, # n.f,
            slice.m = slice.m
        )
    }
    new.samples <- list(
        e.t = e.t, sigma.t = sigma.t, eta.t = eta.t,
        theta.t = theta.t, alpha.t = alpha.t, b.t = b.t,
        pi.t = pi.t, mu.t = mu.t, xi.t = xi.t, beta.t = beta.t
    )
    return(new.samples)
}

sample.paras.chains <- function(jchain, Nsamples, paras.init, link = TRUE, fz.idx = 1,
                                z_mu = 0, z_sigma = 0.7, fz = NULL, fz.prime = NULL,
                                y, preds.all, k1, k2, a.sigma = 1, b.sigma = 1, a.DP = 2,
                                mean.mu = 0, sd.mu = 1, shape.xi = 2, scale.xi = 2,
                                lam = 1 / 2, p.kl = NULL, n.f = 2,
                                lb.list, ub.list, slice.w, slice.m = 1e4, r = 1,
                                eta.scale = 2) {
    # jchain: chain index
    # Nsamples: number of iterations in each chain
    # paras.init: a list of the initials for the parameters
    ### Parameters of Interest:
    ###     theta.t: vector of length p-1
    ###     alpha.t: vector of length p
    ###     beta.t: vector of length 2*L
    ###     sigma.t: scalar
    ###     b.t: vector of length K-1
    ###     pi.t: vector of length K, (pi_1, ..., pi_K)
    ###     mu.t: vector of length K,
    ###     xi.t: vector of length K

    set.seed(jchain + 2023)
    sigma.t <- paras.init$sigma.t
    theta.t <- paras.init$theta.t
    alpha.t <- paras.init$alpha.t
    b.t <- paras.init$b.t
    pi.t <- paras.init$pi.t
    mu.t <- paras.init$mu.t
    xi.t <- paras.init$xi.t

    n <- length(y)
    K <- length(pi.t)
    p <- length(alpha.t)
    if (p < dim(preds.all)[1]) {
        beta.t <- paras.init$beta.t
        # L <- round(length(beta.t) / n.f)
        L <- length(beta.t)
        eta.samples <- beta.samples <- matrix(nrow = Nsamples, ncol = L)
        if (is.null(p.kl)) {
            p.kl <- 1 - 1 / (2^(L - 1)) # default value
        }
    } else {
        beta.t <- 0
        eta.samples <- beta.samples <- 0
        p.kl <- 0
    }

    e.samples <- matrix(nrow = Nsamples, ncol = n)
    sigma.samples <- numeric(Nsamples)
    theta.samples <- matrix(nrow = Nsamples, ncol = p - 1)
    alpha.samples <- matrix(nrow = Nsamples, ncol = p)
    if (link) {
        b.samples <- pi.samples <- mu.samples <- xi.samples <- NULL
    } else {
        b.samples <- matrix(nrow = Nsamples, ncol = K - 1)
        pi.samples <- mu.samples <- xi.samples <- matrix(nrow = Nsamples, ncol = K)
    }

    start <- Sys.time()
    for (i in 1:Nsamples) {
        if (i %% 1000 == 0) {
            print(paste("i =", i, ", elapsed time:", Sys.time() - start))
        }
        new.samples <- sample.paras.single(theta.t, alpha.t,
            beta.t, sigma.t, link, fz.idx,
            z_mu, z_sigma, fz, fz.prime,
            b.t, pi.t, mu.t, xi.t, y, preds.all, k1, k2, a.sigma, b.sigma, a.DP,
            mean.mu, sd.mu, shape.xi, scale.xi,
            lam = 1 / 2, p.kl, n.f, lb.list, ub.list, slice.w, slice.m, r,
            eta.scale
        )
        e.samples[i, ] <- new.samples$e.t
        sigma.samples[i] <- sigma.t <- new.samples$sigma.t
        theta.samples[i, ] <- theta.t <- new.samples$theta.t
        alpha.samples[i, ] <- alpha.t <- new.samples$alpha.t
        if (!link) {
            b.samples[i, ] <- b.t <- new.samples$b.t
            pi.samples[i, ] <- pi.t <- new.samples$pi.t
            mu.samples[i, ] <- mu.t <- new.samples$mu.t
            xi.samples[i, ] <- xi.t <- new.samples$xi.t
        }
        if (p < dim(preds.all)[1]) {
            eta.samples[i, ] <- new.samples$eta.t
            beta.samples[i, ] <- beta.t <- new.samples$beta.t
        }
    }
    all.samples <- list(
        e.samples = e.samples, sigma.samples = sigma.samples, eta.samples = eta.samples,
        theta.samples = theta.samples, alpha.samples = alpha.samples, b.samples = b.samples,
        pi.samples = pi.samples, mu.samples = mu.samples, xi.samples = xi.samples,
        beta.samples = beta.samples
    )
    return(all.samples)
}



MCMC_estimators_spike_slab <- function(x, y, w, nbasis = 15, n.f = 2, tau, r = 1,
                                       K = 10, seed = 2023, p.kl = 0.5, lam = 0.5,
                                       n.chains = 2, Nsamples = 1e5, slice.w = 1, slice.m = 1,
                                       lb.list = NULL, ub.list = NULL, link = FALSE, fz.idx = 2,
                                       fz = NULL, fz.prime = NULL,
                                       a.sigma = 1, b.sigma = 1, a.DP = 2, mean.mu = 0, sd.mu = 10,
                                       shape.xi = 2, scale.xi = 3, eta.scale = 2,
                                       burnin = 2e4, thinning = 100,
                                       psrf = FALSE, autocorr = FALSE, trace = FALSE, density = FALSE,
                                       effsize = TRUE, effsize.fz = FALSE, save.samples = FALSE,
                                       paras.init = NULL) {
    k1 <- (1 - 2 * tau) / (tau * (1 - tau))
    k2 <- 2 / (tau * (1 - tau))
    preds.all <- rbind(x, w)

    samples <- list()
    for (j.chain in 1:n.chains) {
        if (is.null(paras.init)) {
            print("initialize...")
            set.seed(seed + j.chain)
            paras.init <- list()
            paras.init$theta.t <- theta.t <- c(runif(p - 2, -pi / 2, pi / 2), runif(1, 0, 2 * pi))
            paras.init$alpha.t <- alpha.transform(theta.t, r = r)
            paras.init$beta.t <- rnorm(nbasis * n.f, 0, 1)
            paras.init$sigma.t <- 1 / (rgamma(1, shape = a.sigma + 1, rate = b.sigma))
            paras.init$b.t <- b.t <- rbeta(K - 1, 1, a.DP)
            paras.init$pi.t <- pi.transform(b.t)
            paras.init$mu.t <- rnorm(K, mean.mu, sd.mu) # rep(mu0, K) # rnorm(K, 0, 1)
            paras.init$xi.t <- 1 / (rgamma(K, shape = shape.xi, rate = scale.xi)) # rep(xi0, K) # rinvgamma(K, 1 + 1, 1)
        }
        print(paras.init)

        print("start sampling...")
        samples[[j.chain]] <- sample.paras.chains(
            j.chain,
            Nsamples = Nsamples, paras.init,
            link = link, fz.idx = fz.idx, z_mu = 0, z_sigma = 1,
            fz = fz, fz.prime = fz.prime, y, preds.all, k1, k2,
            a.sigma, b.sigma, a.DP, mean.mu, sd.mu,
            shape.xi, scale.xi, lam, p.kl = p.kl,
            n.f = n.f, lb.list = lb.list, ub.list = ub.list,
            slice.w = slice.w, slice.m = slice.m, r = r,
            eta.scale = eta.scale
        )
    }

    ### diagnosis on the samples ####
    paras <- c("sigma", "theta", "alpha")
    if (!is.null(w)) {
        paras <- c(paras, "eta", "beta")
    }
    if (!link) {
        paras <- c(paras, "b", "mu", "xi")
    }
    res0 <- mcmc.diagnose(samples, paras, Nsamples,
        burnin = burnin, thinning = thinning,
        psrf = psrf, autocorr = autocorr, trace = trace, effsize = effsize,
        est = "median"
    )

    #### (1) estimators for all chains ####
    sigma.est <- res0$all.estimators$sigma # (n.chains)
    p <- dim(x)[1]
    theta.est <- res0$all.estimators$theta # (n.chains, p-1)
    theta.est0 <- colMeans(theta.est)
    alpha.est <- res0$all.estimators$alpha # (n.chains, p)
    alpha.est0 <- colMeans(alpha.est)
    if (!is.null(w)) {
        beta.est <- res0$all.estimators$beta # (n.chains, nbasis * nf)
        beta.est0 <- colMeans(beta.est)
        z.est <- alpha.est0 %*% x + beta.est0 %*% w # (n.chains, n)
    } else {
        z.est <- alpha.est0 %*% x # (n.chains, n)
        beta.est <- NULL
    }

    #### (2) f(z) ####
    if (!link) {
        b.est <- res0$all.estimators$b # (n.chains, K-1)
        mu.est <- res0$all.estimators$mu # (n.chains, K)
        xi.est <- res0$all.estimators$xi # (n.chains, K)
        b.est0 <- colMeans(b.est)
        mu.est0 <- colMeans(mu.est)
        xi.est0 <- colMeans(xi.est)

        res.fz <- mcmc.linkf(samples, x, w, Nsamples, K, burnin, thinning,
            psrf = psrf, autoburnin = FALSE, multivariate = TRUE,
            autocorr = autocorr, trace = trace,
            density = density, effsize = effsize.fz
        )
        fz.est <- t(res.fz$fz.sample.median) # (n.chains, n)
        fz.est0 <- colMeans(fz.est)
    } else {
        fz.est0 <- b.est0 <- mu.est0 <- xi.est0 <- NULL
    }

    ##### results ####
    res.list <- list(
        res0 = res0, sigma.est = sigma.est, theta.est = theta.est0,
        alpha.est = alpha.est0, beta.est = beta.est0,
        z.est = z.est, fz.est = fz.est0, paras.init = paras.init,
        b.est = b.est0, mu.est = mu.est0, xi.est = xi.est0
    )
    if (save.samples) {
        res.list$samples <- samples
    }
    return(res.list)
}
