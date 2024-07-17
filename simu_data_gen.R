##### NOTE: FUNCTIONAL PREDICTORS ARE GENERATED USING gamma(t) AS THE MEAN CURVE.

#### Separate alpha & beta, constrain norm(alpha) = r1 & norm(beta) = r2 ####
## r = 1 is the mose common case, while in high dimensional case,
## we can release this constraint and set r to be some large constant

#### libraries ####
source("utilities.R")
require(pracma)
require(ald)

##### Data generation #####
simu.data.gen <- function(
    n = 200, N = 100, p = 4, r1 = 1, r2 = 10, func.p = TRUE, tau = 0.5,
    fz.idx = 1, n.f = 2, norder = 4, nbasis = 7, seed = 2023,
    umin = 0.5, umax = 1.5, z_mu = NULL, z_sigma = NULL, z0 = 0,
    srvf.lambda = 0,
    showplot.srvf = FALSE, curves.visualize = TRUE,
    response.plot = FALSE, theta2.case = 1, bss.sd1 = 0.5, bss.sd2 = 0.5,
    beta.constrain = TRUE, beta.sparse = TRUE, curve_center = FALSE) {
    #### (1) scalar predictors & coefficients ####
    ##### (a) scalar predictors ####
    set.seed(seed)
    x <- matrix(nrow = p, ncol = n)
    for (j in 1:p) {
        x[j, ] <- runif(n, -1, 1)
    }

    ##### (b) scalar coefficients ####
    if (p == 4) {
        theta1 <- c(pi / 6, 0.6154797, 5 * pi / 4)
    } else if (p == 2) {
        theta1 <- runif(1, 0, 2 * pi)
    } else if (p > 2) {
        theta1 <- c(runif(p - 2, -pi / 2, pi / 2), runif(1, 0, 2 * pi))
    }
    alpha.true <- alpha.transform(theta1, r = r1)
    ax <- crossprod(x, alpha.true)

    #### (2) Functional predictors & coefficients ####
    if (func.p) { # exists functional predictors
        require(fda)
        require(fdasrvf)
        ### (a) B-Spline coefficients ####
        if (beta.constrain) { # constrain norm(beta)=r2
            if (n.f == 1) {
                if (theta2.case == 1) {
                    theta2 <- c(
                        0.15 * pi, -pi / 6, -pi / 5, -pi * 0.1, 0.1 * pi, pi / 3
                    ) # 6
                } else if (theta2.case == 2) {
                    theta2 <- c(
                        0.1 * pi, 0, 0, 0, -0.3 * pi, pi / 3
                    )
                } else if (theta2.case == 3) {
                    # theta2 <- c(0.2, -0.3, -0.6, -0.2, 0.5, 1)
                    # theta2 <- c(0, -0.3, -0.6, -0.2, 0.5, 1)
                    theta2 <- c(0, -0.4, -0.7, -0.5, 0.5, 1.5)
                }
            } else if (n.f == 2) {
                if (theta2.case == 1) {
                    theta2 <- c(
                        0.05 * pi, -pi / 12, -pi / 8, -pi * 0.05, 0.07 * pi, pi / 10, pi / 12, # for p1
                        -0.5, -0.5, -0.3, 0.6, 1, 2.5
                    )
                } else if (theta2.case == 2) {
                    theta2 <- c(
                        0.05, -0.2, -0.4, -0.2, 0.25, 0.33, 0.21,
                        -0.35, -0.56, -0.33, 0.7, 1.2, 3
                    )
                } else if (theta2.case == 3) {
                    theta2 <- c(
                        0, -0.2, -0.35, -0.2, 0.25, 0.3, 0.15,
                        -0.35, -0.55, -0.33, 0.9, 1.2, 2.5
                    )
                }
            }
            beta.true <- alpha.transform(theta2, r = r2)
        } else {
            theta2 <- NULL
            # beta.true <- numeric(length = n.f * nbasis)
            if (beta.sparse) {
                beta.true <- c(
                    rep(0, nbasis %/% 3), c(3), rep(0, nbasis - nbasis %/% 3 - 1),
                    rep(0, nbasis %/% 3 * 2), c(-4), rep(0, nbasis - nbasis %/% 3 * 2 - 1)
                )
                beta.true <- beta.true[1:(nbasis * n.f)]
                # beta.true[3:4] = c(3,4)
                # beta.true[nbasis + (6:7)] = c(-4, -3)
            } else{
                beta.true <- numeric(nbasis * n.f)
                beta.true[3:5] <- c(-4,-5,4)
                beta.true[(3:5)+nbasis] <- c(-5,-3,4)
            }
        }

        ### (b) functional coefficients ####
        t <- seq(0, 1, length.out = N)
        psi0 <- create.bspline.basis(norder = norder, nbasis = nbasis)
        eval.psi0 <- eval.basis(t, psi0) # (T, L)
        gamma.coef <- matrix(beta.true, nrow = nbasis, ncol = n.f) # (n.f, L)
        gamma.true <- eval.psi0 %*% gamma.coef # (T, n.f)

        ### (c) functional predictors ####
        #### (i) initial curves ####
        bss <- rbind(sin(2 * pi * t), cos(2 * pi * t)) # * 2 # (2, T)
        curves <- array(dim = c(N, n, n.f))
        for (k in 1:n.f) {
            gamma.true[, k] <- tcrossprod(gamma.coef[, k], eval.psi0)
            if (curve_center) {
                curves.mu <- matrix(0, nrow = 1, ncol = N)
            } else {
                curves.mu <- matrix(gamma.true[, k], nrow = 1, ncol = N)
            }
            curves.eta <- mvrnorm(n, rep(0, nbasis), diag(rep(bss.sd1^2, nbasis))) %*% t(eval.psi0)
            # curves.eta <- cbind(rnorm(n, 0, bss.sd1), rnorm(n, 0, bss.sd2)) %*% bss
            # curves.err <- mvrnorm(N, rep(0, n), diag(0.001, n, n))
            curves.err <- 0
            curves[, , k] <- t(pracma::repmat(curves.mu, n, 1) + curves.eta) + curves.err
        }
        #### (ii) true predictors in the model ####
        start <- Sys.time()
        curves.true <- array(dim = c(N, n, n.f))
        for (k in 1:n.f) {
            srvf.res <- fdasrvf::multiple_align_functions(curves[, , k], t, gamma.true[, k],
                lambda = srvf.lambda, showplot = showplot.srvf
            )
            curves.true[, , k] <- srvf.res$fn
        }
        print(Sys.time() - start)
        bw <- numeric(n)
        for (k in 1:n.f) {
            bw <- bw + crossprod(curves.true[, , k], gamma.true[, k]) / N
        }

        #### (iii) random distorted curves as the observed curves ####
        curves.random <- distort.curves(curves.true, t, seed, umin, umax)

        ### visualization ####
        if (curves.visualize) {
            ### registration on distorted curves ####
            start <- Sys.time()
            pracma::fprintf("lambda = %.2f\n", srvf.lambda)
            curves.reg.gamma <- array(dim = c(N, n, n.f))
            for (k in 1:n.f) {
                srvf.res <- fdasrvf::multiple_align_functions(curves.random[, , k], t, gamma.true[, k],
                    lambda = srvf.lambda, showplot = showplot.srvf
                )
                curves.reg.gamma[, , k] <- srvf.res$fn
            }
            print(Sys.time() - start)

            ### plot ####
            par(mfrow = c(n.f, 4))
            for (k in 1:n.f) {
                matplot(t, curves[, , k],
                    type = "l", main = paste0("initial c_", k, "(t)"),
                    ylim = range(c(range(curves.true[, , k]), range(gamma.true[, k])))
                )
                lines(t, gamma.true[, k], lwd = 3, col = 1, lty = 1)

                matplot(t, curves.true[, , k],
                    type = "l", main = paste0("true c_", k, "(t)"),
                    ylim = range(c(range(curves.true[, , k]), range(gamma.true[, k])))
                )
                lines(t, gamma.true[, k], lwd = 3, col = 1, lty = 1)

                matplot(t, curves.random[, , k],
                    type = "l", main = paste0("distorted c_", k, "(t)"),
                    ylim = range(c(range(curves.true[, , k]), range(gamma.true[, k])))
                )
                lines(t, gamma.true[, k], lwd = 3, col = 1, lty = 1)

                matplot(t, curves.reg.gamma[, , k],
                    type = "l", main = paste0("register c_", k, "(t) on gamma(t)"),
                    ylim = range(c(range(curves.true[, , k]), range(gamma.true[, k])))
                )
                lines(t, gamma.true[, k], lwd = 3, col = 1, lty = 1)
            }
            legend("bottom",
                legend = c(expression(gamma(t))),
                col = 1:2, lty = 1:2, lwd = 3
            )
        }
    } else { # only scalar covariates
        theta2 <- beta.true <- bw <- 0
        curves <- curves.true <- curves.random <- gamma.true <- beta.true <- t <- NULL
    }


    #### (3) z & link function & Q #####
    z <- ax + bw
    fprintf("Range of z: (%.3G, %.3G)\n", range(z)[1], range(z)[2])

    #### link function ####
    if (is.null(z_mu)) {
        z_mu <- round(mean(z), 1)
    }
    if (is.null(z_sigma)) {
        z_sigma <- max(floor(sd(z)), 0.5)
    }
    print(paste0("z_mu = ", z_mu, ", z_sigma = ", z_sigma), digits = 4)

    if (fz.idx == 0) { # exp
        fz <- function(z) 1 - exp(-(z - z_mu) * z_sigma)
    } else if (fz.idx == 1) { # normal
        fz <- function(z) pnorm(z, mean = z_mu, sd = z_sigma)
    } else if (fz.idx == 2) { # Laplace
        fz <- function(z) ald::pALD(z, mu = z_mu, sigma = z_sigma / 2, p = 0.5)
        # } else if (fz.idx == 3) { # exp-log
        #     pexplog <- function(z, p, beta) {
        #         1 - log(1 - (1 - p) * exp(-beta * z)) / log(p)
        #     }
        #     fz <- function(z) pexplog(z - z_mu, p = 0.9, beta = z_sigma)
    } else if (fz.idx == 3) { # pareto cdf
        pareto.cdf <- function(z, scale, shape) {
            cdf <- numeric(length = length(z))
            cdf[z >= scale] <- 1 - (scale / z[z >= scale])^shape
            return(cdf)
        }
        fz <- function(z) pareto.cdf(z, scale = z_mu, shape = z_sigma)
    }
    Q <- fz(z + z0)
    fprintf("Range of Q: (%.3G, %.3G)\n", range(Q)[1], range(Q)[2])

    ### (4) ALD errors & response ####
    epsilons.range <- range(c(0, 1) - range(Q))
    epsilons <- NULL
    iters <- 0
    while (length(epsilons) < n && iters <= 100) {
        iters <- iters + 1
        epsilon <- ald::rALD(10 * n, mu = 0, sigma = 0.025, p = tau)
        accept <- epsilon[epsilon >= epsilons.range[1] & epsilon <= epsilons.range[2]]
        epsilons <- c(epsilons, accept)
    }
    epsilons <- epsilons[1:n]
    print(summary(epsilons))

    y <- as.numeric(Q) + epsilons
    if (response.plot) {
        par(mfrow = c(2, 2))
        plot(z, Q, ylim = c(0, 1), main = "f(z)")
        points(z, y, col = 2)
        abline(h = c(0, 1), lwd = 2, col = 2, lty = 2)
        hist(Q, breaks = 20, main = "hist of Q")
        hist(y, breaks = 20, main = "hist of y")
    }

    ### (5) save the variables ####
    data <- list(
        "x" = x, "curves.ini" = curves, "curves.true" = curves.true,
        "curves" = curves.random, "t" = t, "z" = z, "fz" = fz,
        "Q" = Q, "y" = y, "theta1" = theta1, "r1" = r1, "alpha" = alpha.true,
        "theta2" = theta2, "beta" = beta.true, "gamma" = gamma.true,
        "z_mu" = z_mu, "z_sigma" = z_sigma
    )
    return(data)
}
