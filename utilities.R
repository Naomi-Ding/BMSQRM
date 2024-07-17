#### Transformation ####
alpha.transform <- function(theta, r = 1) {
    ### one to one transformation between theta & alpha
    ### p >= 2
    p <- length(theta) + 1
    alpha <- numeric(p)
    alpha[1] <- sin(theta[1])
    tmp <- cos(theta[1])
    if (p == 2) { # polar system
        alpha[2] <- tmp
    } else if (p >= 3) {
        for (j in 2:(p - 1)) {
            alpha[j] <- sin(theta[j]) * tmp
            tmp <- cos(theta[j]) * tmp
        }
        alpha[p] <- tmp
    }
    return(r * alpha)
}

alpha.transform.inverse <- function(alpha) {
    ### one to one transformation between theta & alpha
    ### transform alpha to theta
    ### p >= 2
    p <- length(alpha)
    theta <- numeric(p - 1)
    theta[p - 1] <- atan2(alpha[p - 1], alpha[p])
    if (theta[p - 1] < 0) {
        theta[p - 1] <- theta[p - 1] + 2 * pi
    }
    tmp <- alpha[p]^2 + alpha[p - 1]^2
    if (p >= 3) {
        for (j in (p - 2):1) {
            theta[j] <- atan2(alpha[j], sqrt(tmp))
            tmp <- tmp + alpha[j]^2
        }
        # r <- tmp
    }
    return(theta)
}

pi.transform <- function(b) {
    ### one to one transformation between b & pi
    ### K >= 2
    K <- length(b) + 1
    pi <- numeric(K)
    pi[1] <- b[1]
    if (K == 2) {
        pi[K] <- 1 - b[1]
    } else if (K >= 3) {
        tmp <- 1 - b[1]
        for (j in 2:(K - 1)) {
            pi[j] <- b[j] * tmp
            tmp <- (1 - b[j]) * tmp
        }
        pi[K] <- 1 - sum(pi[1:(K - 1)])
    }
    return(pi)
}


### distort curves & apply B-spline decomposition ####
distort.curves.old <- function(
    curves, t = NULL, seed = 2023,
    umin = 1, umax = 1.5) {
    require(fdasrvf)

    dims <- dim(curves)
    N <- dims[1]
    n <- dims[2]
    d <- dims[3]
    if (is.na(d)) {
        d <- 1
        curves <- array(curves, dim = c(N, n, 1))
        print(dim(curves))
    }
    if (is.null(t)) {
        t <- seq(0, 1, length.out = N)
    }

    set.seed(seed)
    f.random <- array(dim = c(N, n, d))
    h <- array(dim = c(N, n, d))
    for (k in 1:d) {
        u <- runif(n, umin, umax)
        for (i in 1:n) {
            h[, i, k] <- t^u[i]
        }
        hk.mean <- rowMeans(h[, , k])
        hk.inv <- invertGamma(hk.mean)
        for (i in 1:n) {
            h[, i, k] <- warp_f_gamma(h[, i, k], t, hk.inv)
            f.random[, i, k] <- warp_f_gamma(curves[, i, k], t, h[, i, k])
        }
    }
    return(f.random)
}

distort.curves <- function(
    curves, t = NULL, seed = 2023,
    umin = 1, umax = 1.5) {
    require(fdasrvf)

    dims <- dim(curves)
    N <- dims[1]
    n <- dims[2]
    d <- dims[3]
    if (is.na(d)) {
        d <- 1
        curves <- array(curves, dim = c(N, n, 1))
        print(dim(curves))
    }
    if (is.null(t)) {
        t <- seq(0, 1, length.out = N)
    }

    set.seed(seed)
    f.random <- array(dim = c(N, n, d))
    h <- array(dim = c(N, n, d))
    for (k in 1:d) {
        u <- runif(n, umin, umax)
        for (i in 1:n) {
            f.random[, i, k] <- warp_f_gamma(curves[, i, k], t, t^u[i])
        }
    }
    return(f.random)
}


reconstruct.curves <- function(curves, t = NULL, norder = 4, nbasis = 7) {
    require(fda)
    dims <- dim(curves)
    N <- dims[1]
    n <- dims[2]
    d <- dims[3]
    if (is.na(d)) {
        d <- 1
        curves <- array(curves, dim = c(N, n, 1))
        print(dim(curves))
    }
    if (is.null(t)) {
        t <- seq(0, 1, length.out = N)
    }
    psi0 <- create.bspline.basis(rangeval = range(t), norder = norder, nbasis = nbasis)
    eval.psi0 <- eval.basis(t, psi0) # (T, L)

    ### B-spline decompose & new predictors
    w <- matrix(nrow = d * nbasis, ncol = n)
    for (k in 1:d) {
        idx <- (1:nbasis) + (k - 1) * nbasis
        w[idx, ] <- crossprod(eval.psi0, curves[, , k]) / N # (nbasis, n)
    }
    return(w)
}


reconstruct.curves.new <- function(curves, t = NULL, basis_type = "bsp",
                                   rangeval = NULL,
                                   norder = 4, nbasis = 7,
                                   FVEthreshold = 0.9,
                                   usergrid = TRUE) {
    require(fda)
    dims <- dim(curves)
    N <- dims[1]
    n <- dims[2]
    d <- dims[3]
    if (is.na(d)) {
        d <- 1
        curves <- array(curves, dim = c(N, n, 1))
        print(dim(curves))
    }
    if (is.null(t)) {
        t <- seq(0, 1, length.out = N)
    }
    if (is.null(rangeval)) {
        rangeval <- range(t)
    }
    if (basis_type == "pca") {
        require(fdapace)
        w <- NULL
        basis <- NULL
        for (k in 1:d) {
            fdata <- MakeFPCAInputs(
                IDs = rep(1:n, each = N), tVec = rep(t, n),
                curves[, , k]
            )
            fpca.dense <- FPCA(
                fdata$Ly, fdata$Lt,
                list(FVEthreshold = FVEthreshold, usergrid = usergrid)
            )
            phi <- fpca.dense$phi # (N, npc)
            xi <- t(fpca.dense$xiEst) # (npc, n)
            rownames(xi) <- colnames(phi) <- c(paste0("c", k, ".pc", 1:fpca.dense$selectK))
            basis <- cbind(basis, phi)
            w <- rbind(w, xi)
        }
        res <- list(w = w, basis = basis)
    } else {
        if (basis_type == "bsp") {
            psi0 <- create.bspline.basis(rangeval = rangeval, norder = norder, nbasis = nbasis)
        } else if (basis_type == "fourier") {
            psi0 <- create.fourier.basis(rangeval = rangeval, nbasis = nbasis)
        }
        eval.psi0 <- eval.basis(t, psi0) # (N, L)
        ### Basis decompose & new predictors
        w <- matrix(nrow = d * nbasis, ncol = n)
        for (k in 1:d) {
            idx <- (1:nbasis) + (k - 1) * nbasis
            w[idx, ] <- crossprod(eval.psi0, curves[, , k]) / N # (nbasis, n)
        }
        res <- list(w = w, basis = eval.psi0)
    }
    return(res)
}

######### Diagnosis ################
mcmc.diagnose <- function(samples, paras, Nsamples, burnin = 5000, thinning = 2,
                          psrf = TRUE, autoburnin = FALSE, multivariate = TRUE,
                          autocorr = TRUE, trace = TRUE, density = FALSE,
                          effsize = TRUE, est = c("median", "mean")) {
    require(coda)
    n.chains <- length(samples)
    Ind <- burnin + (1:((Nsamples - burnin) / thinning)) * thinning

    num.paras <- length(paras)
    all.chains <- list()
    all.estimators <- list()

    for (pp in 1:num.paras) {
        print(paste0("parameter: ", paras[pp]))
        para <- paste0(paras[pp], ".samples")
        chains <- mcmc.list()
        for (i in 1:n.chains) {
            if (paras[pp] == "sigma") {
                chains[[i]] <- mcmc(data = samples[[i]][[para]][Ind])
            } else {
                chains[[i]] <- mcmc(data = samples[[i]][[para]][Ind, ])
            }
        }
        all.chains[[paras[pp]]] <- chains

        ### convergence diagnostic
        if (psrf) {
            if (n.chains > 1) {
                print(gelman.diag(chains, autoburnin, multivariate))
            }
        }

        if (autocorr) {
            autocorr.plot(chains)
        }
        if (trace) {
            plot(chains, density = density)
        }
        if (effsize) {
            print(paste0("Effective sample size: ", list(effectiveSize(chains))))
            print(paste0("Total number of samples: ", n.chains * length(Ind)))
        }

        if (is.null(dim(chains[[1]]))) {
            estimator <- as.numeric(lapply(chains, est))
            names(estimator) <- paste0("Chain ", 1:n.chains)
        } else {
            estimator <- matrix(
                nrow = n.chains, ncol = dim(chains[[1]])[2],
                dimnames = list(paste0("Chain ", 1:n.chains), paste0(paras[pp], "_", 1:dim(chains[[1]])[2])) # nolint: seq_linter, line_length_linter.
            )
            for (i in 1:n.chains) {
                estimator[i, ] <- apply(chains[[i]], 2, est)
            }
        }
        all.estimators[[paras[pp]]] <- estimator
    }

    res <- list(all.chains = all.chains, all.estimators = all.estimators)
    return(res)
}


### link function
mcmc.linkf <- function(samples, x, w, Nsamples, K = 15, burnin = 5000, thinning = 2, psrf = TRUE, autoburnin = FALSE,
                       multivariate = TRUE, autocorr = FALSE, trace = FALSE, density = FALSE,
                       effsize = TRUE) {
    require(coda)
    n.chains <- length(samples)
    Ind <- burnin + (1:((Nsamples - burnin) / thinning)) * thinning
    nInd <- length(Ind)
    n <- dim(x)[2]
    fz.samples <- array(dim = c(nInd, n, n.chains))
    fz.sample.median <- matrix(nrow = n, ncol = n.chains)
    for (i in 1:n.chains) {
        alpha.sample <- samples[[i]][["alpha.samples"]][Ind, ] # (#Ind, p)
        if (!is.null(w)) {
            beta.sample <- samples[[i]][["beta.samples"]][Ind, ]
            z.sample <- alpha.sample %*% x + beta.sample %*% w
        } else {
            z.sample <- alpha.sample %*% x # (#Ind, n)
        }
        pi.sample <- samples[[i]][["pi.samples"]][Ind, ] # (#Ind, K)
        mu.sample <- samples[[i]][["mu.samples"]][Ind, ] # (#Ind, K)
        xi.sample <- samples[[i]][["xi.samples"]][Ind, ] # (#Ind, K)
        fz.sample <- matrix(nrow = nInd, ncol = n)
        for (j in 1:nInd) {
            for (ii in 1:n) {
                fz.sample[j, ii] <- sum(pi.sample[j, ] * pnorm(z.sample[j, ii], mu.sample[j, ], sqrt(xi.sample[j, ])))
            }
        }
        fz.sample.median[, i] <- apply(fz.sample, 2, median)
        fz.samples[, , i] <- fz.sample
    }
    fz.sample.median.mean <- rowMeans(fz.sample.median) # (n)

    if (trace) {
        subsamples <- 1:8
        fz.samples.mcmc <- mcmc.list()
        for (i in 1:n.chains) {
            fz.samples.mcmc[[i]] <- mcmc(fz.samples[, subsamples, i])
        }
        # fz.samples.mcmc <- mcmc.list(mcmc(fz.samples[, subsamples, 1]), mcmc(fz.samples[, subsamples, 2]))
        par(mfrow = c(4, 4))
        plot(fz.samples.mcmc, auto.layout = FALSE, density = density, main = expression(f(z[j])))
    }
    if (psrf || effsize) {
        fz.samples.mcmc <- mcmc.list()
        for (j in 1:n.chains) {
            fz.samples.mcmc[[j]] <- mcmc(fz.samples[, , j])
        }
        if (psrf && n.chains > 1) {
            psrf.fz <- gelman.diag(fz.samples.mcmc, autoburnin = FALSE)
            print(paste0("PSRF for f(z): ", psrf.fz$mpsrf))
        }
        if (effsize) {
            pracma::fprintf("Effective sample size = %s", list(effectiveSize(fz.samples.mcmc)))
            print(paste0("Total number of samples: ", n.chains * nInd))
        }
    }

    res <- list(
        fz.samples = fz.samples, fz.sample.median = fz.sample.median,
        fz.sample.median.mean = fz.sample.median.mean
    )
    return(res)
}


#### evaluate the errors #####
errors_evaluation <- function(res.est, paras.true, z, Q, func.p = TRUE, link = FALSE,
                              pred.response = TRUE, n.f = 2, norder = 4, nbasis = 7,
                              visualization = TRUE) {
    ### res.est: results by the function MCMC_estimators
    ### paras.true: a list of ground truth, (alpha, beta, gamma, fz)
    require(fda)
    require(fdasrvf)
    require(pracma)

    errors <- errors.mean <- NULL

    # alpha
    alpha <- paras.true$alpha
    p <- length(alpha)
    alpha.all.est <- res.est$alpha.all.est
    n.chains <- dim(alpha.all.est)[1]
    if (n.chains > 1) {
        sse.alpha <- rowSums((alpha.all.est[, 1:p] - pracma::repmat(alpha, n.chains, 1))^2)
        errors <- cbind(errors, sse.alpha)
    }
    sse.alpha.mean <- sum((res.est$alpha.est - alpha)^2)
    errors.mean <- c(errors.mean, sse.alpha.mean)

    if (func.p) {
        # beta
        beta <- paras.true$beta
        if (n.chains > 1) {
            beta.est <- alpha.all.est[, (1:(n.f * nbasis) + p)] # (n.chains, nbasis * n.f)
            sse.beta <- rowSums((beta.est - pracma::repmat(t(beta), n.chains, 1))^2)
            errors <- cbind(errors, sse.beta)
        }
        beta.est.mean <- res.est$beta.est # vector (n.f * nbasis)
        sse.beta.mean <- sum((beta.est.mean - beta)^2) # scalar

        # gamma(t)
        gamma <- paras.true$gamma
        N <- dim(gamma)[1]
        t <- seq(0, 1, length.out = N)
        psi0 <- fda::create.bspline.basis(norder = norder, nbasis = nbasis)
        eval.psi0 <- fda::eval.basis(t, psi0) # (T, L)
        if (n.chains > 1) {
            gamma.est <- array(dim = c(N, n.chains, n.f))
            gamma.est.mean <- matrix(nrow = N, ncol = n.f)
            ise.gamma <- matrix(nrow = n.chains, ncol = n.f)
            for (k in 1:n.f) {
                idx <- (1:nbasis) + (k - 1) * nbasis
                gamma.est[, , k] <- tcrossprod(eval.psi0, beta.est[, idx])
                gamma.est.mean[, k] <- eval.psi0 %*% beta.est.mean[idx]
                # ISE of gamma(t)
                ise.gamma[, k] <- colMeans((gamma.est[, , k] - gamma[, k])^2)
            }
            errors <- cbind(errors, ise.gamma)
        } else if (n.chains == 1) {
            gamma.est.mean <- matrix(nrow = N, ncol = n.f)
            for (k in 1:n.f) {
                idx <- (1:nbasis) + (k - 1) * nbasis
                gamma.est.mean[, k] <- eval.psi0 %*% beta.est.mean[idx]
            }
            gamma.est <- gamma.est.mean
        }
        ise.gamma.mean0 <- colMeans((gamma.est.mean - gamma)^2)
        errors.mean <- c(
            errors.mean, sse.beta.mean, # sse.beta.mean1,
            ise.gamma.mean0 # , ise.gamma.mean1
        )
        ## visualization
        if (visualization) {
            par(mfrow = c(1, n.f))
            if (n.chains > 1) {
                for (k in 1:n.f) {
                    matplot(t, gamma.est[, , k],
                        type = "l", lwd = 3, lty = 2, col = 2:(n.chains + 1),
                        main = paste0("gamma_", k, "(t)"),
                        ylim = range(c(gamma.est[, , k], gamma[, k]))
                    )
                    lines(t, gamma[, k], lwd = 3, col = 1, lty = 1)
                    lines(t, gamma.est.mean[, k], lwd = 4, col = n.chains + 2, lty = 3)
                    legend("bottomright",
                        legend = c(paste0("chain_", 1:n.chains), "true", "avg_chains"),
                        lwd = c(rep(3, (n.chains + 1)), 4), lty = c(rep(2, n.chains), 1, 3),
                        col = c(2:(n.chains + 1), 1, (n.chains + 2))
                    )
                }
            } else if (n.chains == 1) {
                for (k in 1:n.f) {
                    plot(t, gamma[, k],
                        type = "l", main = paste0("gamma_", k, "(t)"),
                        lwd = 3, col = 1, lty = 1, ylim = range(c(gamma.est.mean[, k], gamma[, k]))
                    )
                    lines(t, gamma.est.mean[, k], lwd = 4, col = 2, lty = 2)
                    legend("bottomright",
                        legend = c("true", "chain_1"),
                        lwd = 3:4, lty = 1:2, col = 1:2
                    )
                }
            }
        }
    } else {
        gamma.est <- gamma.est.mean <- NULL
    }

    ## f(z)
    if (!link) {
        pi.est <- res.est$fz.est.paras$pi.est
        mu.est <- res.est$fz.est.paras$mu.est
        xi.est <- res.est$fz.est.paras$xi.est
        fz.hat <- res.est$fz.hat.func

        res.f <- t(rbind(pi.est, mu.est, xi.est))
        colnames(res.f) <- paste0(rep(c("pi", "mu", "xi"), each = n.chains), "_chain_", rep(1:n.chains, 3))
        rownames(res.f) <- paste0("var_", 1:(dim(pi.est)[2]))
        print(res.f, digits = 4)

        n <- length(z)
        fz.est <- matrix(nrow = n.chains, ncol = n)
        z.sort <- sort(z)
        for (i in 1:n) {
            fz.est[, i] <- fz.hat(z.sort[i])
        }
        fz <- paras.true$fz
        diff.fz <- fz.est - pracma::repmat(fz(z.sort), n.chains, 1)
        mse.fz <- rowMeans(diff.fz^2)
        names(mse.fz) <- c(paste0("mse.fz_chain", 1:n.chains))
        mse.fz.mean <- mean((colMeans(fz.est) - fz(z.sort))^2)
        if (n.chains > 1) {
            errors <- cbind(errors, mse.fz)
        }
        errors.mean <- c(errors.mean, mse.fz.mean)

        ## visualization of f(z) & estimators
        if (visualization) {
            par(mfrow = c(1, 2))
            matplot(z.sort, t(fz.est),
                type = "l", lwd = 3,
                main = expression(fz),
                ylim = c(0, 1), col = 2:(n.chains + 1), lty = 2,
                xlab = "z", ylab = "f(z)"
            )
            lines(z.sort, fz(z.sort), lwd = 3, lty = 1, col = 1)
            if (n.chains > 1) {
                lines(z.sort, colMeans(fz.est), lwd = 4, lty = 3, col = n.chains + 2)
                legend("bottomright",
                    legend = c(paste0("chain_", 1:n.chains), "true", "avg"),
                    lty = c(rep(2, n.chains), 1, 3), col = c(2:(n.chains + 1), 1, n.chains + 2),
                    lwd = c(rep(3, n.chains + 1), 4)
                )
            } else if (n.chains == 1) {
                legend("bottomright", legend = c("chain_1", "true"), lty = c(2, 1), col = c(2, 1), lwd = 3)
            }
        }
    }

    ## prediction of response
    if (pred.response) {
        if (n.chains > 1) {
            z.diff <- res.est$z.est - t(pracma::repmat(z, 1, n.chains))
            mse.z <- rowMeans(z.diff^2)
            names(mse.z) <- c(paste0("mse.z_chain", 1:n.chains))
            q.hat <- res.est$q.hat
            diff.q <- pracma::repmat(t(Q), n.chains, 1) - q.hat # (n.chains, n)
            mse.q <- rowMeans(diff.q^2)
            names(mse.q) <- c(paste0("mse.q_chain", 1:n.chains))
            errors <- cbind(errors, mse.z, mse.q)
        }

        z.est.mean <- res.est$z.est.mean
        mse.z.mean <- mean((z.est.mean - t(z))^2)
        q.hat.mean <- res.est$q.hat.mean
        mse.q.mean <- mean((q.hat.mean - Q)^2)
        errors.mean <- c(errors.mean, mse.z.mean, mse.q.mean)

        ## visualization
        if (visualization) {
            plot(z, Q,
                col = 1, pch = 1, cex = 2, main = "Q_hat",
                xlim = range(c(res.est$z.est, z)), ylim = range(c(res.est$q.hat, Q))
            )
            for (j in 1:n.chains) {
                points(res.est$z.est[j, ], res.est$q.hat[j, ], col = j + 1, pch = j + 1)
            }
            if (n.chains > 1) {
                points(res.est$z.est.mean, q.hat.mean, col = n.chains + 2, pch = n.chains + 2)
                legend("bottomright",
                    legend = c("true", paste0("chain_", 1:n.chains), "avg"),
                    col = 1:(n.chains + 2), pch = 1:(n.chains + 2)
                )
            } else if (n.chains == 1) {
                legend("bottomright", legend = c("true", "chain_1"), col = 1:2, pch = 1:2)
            }
        }
    }

    ### errors ###
    errors.mean.names <- c("sse.alpha")
    if (func.p) {
        errors.mean.names <- c(errors.mean.names, "sse.beta", paste0("ise.gamma_", 1:n.f))
    }
    if (!link) {
        errors.mean.names <- c(errors.mean.names, "mse.f")
    }
    if (pred.response) {
        errors.mean.names <- c(errors.mean.names, "mse.z", "mse.q")
    }
    names(errors.mean) <- errors.mean.names
    if (n.chains > 1) {
        colnames(errors) <- errors.mean.names
        rownames(errors) <- c(paste0("chain_", 1:n.chains))
        print(errors, digits = 4)
    } else {
        errors <- NULL
    }
    print(errors.mean, digits = 4)

    res <- list(
        "errors" = errors, "errors.mean" = errors.mean,
        "gamma.est" = gamma.est, "gamma.mean.est" = gamma.est.mean
    )
    return(res)
}



#### load the real data #####
realdata_load <- function(result_dir, include_bsl = TRUE, response = "final",
                          func.p = TRUE, fun.pre_reg = "no", fun2 = "lose",
                          func.cut_warmup = FALSE, curve_center = FALSE,
                          concate = FALSE, first_signals = "win",
                          concate_shift = FALSE, basis_expansion = "original",
                          func.ratio = 1, basis_type = "bsp", norder = 4,
                          nbasis = 30, FVEthreshold = 0.95,
                          rescale_w = FALSE) {
    data <- readRDS(paste0(result_dir, "eeg33_orig_data.rds"))

    #### (1) scalar preds & response ####
    vars <- data$variables
    n <- dim(vars)[1] # sample size
    N <- 1024 # function length
    if (include_bsl) {
        x0 <- vars[, c(2:4, 6:8)]
        x0[, 3] <- scale(x0[, 3], center = TRUE, scale = TRUE)
    } else {
        x0 <- vars[, c(2, 3, 6:8)] # (n,p)
    }
    x0[, 1] <- scale(x0[, 1], center = TRUE, scale = TRUE)
    x0 <- t(as.matrix(x0)) # (p, n)
    p <- dim(x0)[1]

    if (response == "final") {
        y0 <- vars$Dys_C_W3
    } else if (response == "diff") {
        y0 <- vars$Dys_C_W3 - vars$Dys_C_W1
    }
    y <- (y0 - min(y0)) / (max(y0) - min(y0))

    #### (2) functional preds ####
    if (func.p) {
        require(fda)
        ##### (i) construct the functional pred ####
        t0 <- data$index
        winf <- data$winf
        losf <- data$losf
        if (fun.pre_reg == "no") { # use original curves
            curves <- array(dim = c(N, n, 2))
            curves[, , 1] <- winf
            if (fun2 == "lose") {
                curves[, , 2] <- losf
            } else if (fun2 == "diff") {
                curves[, , 2] <- winf - losf
            }
        } else { # use pre_registered curves
            if (fun.pre_reg == "registr") {
                win_file <- paste0(result_dir, "register_fpca_win_npc=0.9_Kt=20_Kh=4_fpca_max=500.rds")
                los_file <- paste0(result_dir, "register_fpca_los_npc=0.9_Kt=20_Kh=4_fpca_max=500.rds")
                Ywin.reg <- readRDS(win_file)
                Ylos.reg <- readRDS(los_file)
                t_win <- Ywin.reg$Y$t_hat
                t_los <- Ylos.reg$Y$t_hat
                winf.reg <- winf
                losf.reg <- losf
                for (i in 1:n) {
                    winf.reg[, i] <- spline(x = t_win[(i - 1) * N + (1:N)], y = winf[, i], xout = t0)$y
                    losf.reg[, i] <- spline(x = t_los[(i - 1) * N + (1:N)], y = losf[, i], xout = t0)$y
                }
                rm(Ywin.reg, Ylos.reg, t_win, t_los, win_file, los_file)
            } else if (fun.pre_reg == "fdasrvf") {
                win_file <- paste0(result_dir, "register_fpca_win_npc=0.9_Kt=20_Kh=4_fpca_max=500_fdasrvf.rds")
                los_file <- paste0(result_dir, "register_fpca_los_npc=0.9_Kt=20_Kh=4_fpca_max=500_fdasrvf.rds")
                Ywin.reg <- readRDS(win_file)
                Ylos.reg <- readRDS(los_file)
                winf.reg <- Ywin.reg$winf.srvf.0.1
                losf.reg <- Ylos.reg$losf.srvf.0.1
                rm(Ywin.reg, Ylos.reg, win_file, los_file)
            }
            curves <- array(dim = c(N, n, 2))
            curves[, , 1] <- winf.reg
            if (fun2 == "lose") {
                curves[, , 2] <- losf.reg
            } else if (fun2 == "diff") {
                curves[, , 2] <- winf.reg - losf.reg
            }
            rm(winf.reg, losf.reg)
        }
        rm(winf, losf, data)

        if (func.cut_warmup) {
            curves <- curves[-(1:200), , ]
            N <- N - 200
            t0 <- t0[-(1:200)]
        }
        if (curve_center) {
            curves.mean <- apply(curves, 1, colMeans) # (n.f, N)
            curves[, , 1] <- curves[, , 1] - curves.mean[1, ]
            curves[, , 2] <- curves[, , 2] - curves.mean[2, ]
        }
        if (concate) { # flip win signals & concate them together
            if (first_signals == "win") {
                if (concate_shift) {
                    c1 <- curves[, , 1]
                    c2 <- curves[N:1, , 2]
                    c2.shift <- c2 + repmat(t(c1[N, ] - c2[1, ]), N, 1)
                    curves <- rbind(c1, c2.shift)
                } else {
                    curves <- rbind(curves[N:1, , 1], curves[, , 2])
                }
            } else if (first_signals == "lose") {
                if (concate_shift) {
                    c1 <- curves[, , 2]
                    c2 <- curves[N:1, , 1]
                    c2.shift <- c2 + repmat(t(c1[N, ] - c2[1, ]), N, 1)
                    curves <- rbind(c1, c2.shift)
                } else {
                    curves <- rbind(curves[, , 2], curves[N:1, , 1])
                }
            }
            t0 <- seq(0, 2, length.out = N * 2)
            n.f <- 1
            N <- 2 * N
        } else {
            n.f <- 2
        }

        ##### (ii) construct the predictors using B-splines ####
        # source("utilities.R")
        # t0 <- seq(0, 1, length.out = N)
        if (basis_expansion == "original") {
            rangeval <- NULL
        } else if (basis_expansion == "ext_bs") {
            rangeval <- c(min(t0) - 0.01 * max(t0), max(t0) + 0.01 * max(t0))
        }
        print(rangeval)
        w.basis <- reconstruct.curves.new(curves / func.ratio, t0,
            basis_type = basis_type, rangeval = rangeval,
            norder = norder, nbasis = nbasis,
            FVEthreshold = FVEthreshold,
            usergrid = TRUE
        )
        w <- w.basis$w
        basis <- w.basis$basis
        if (rescale_w) {
            print("rescale w")
            w.mu <- rowMeans(w)
            w.sd <- apply(w, 1, sd)
            w <- w / t(repmat(w.sd, n, 1))
            if (basis_type != "pca") {
                # basis0 = basis
                basis.rescale <- NULL
                for (k in 1:n.f) {
                    idx <- (1:nbasis) + (k - 1) * nbasis
                    basis.rescale <- cbind(
                        basis.rescale,
                        basis / repmat(w.sd[idx], N, 1)
                    )
                }
            } else {
                basis.rescale <- basis / repmat(w.sd, N, 1)
                if (!curve_center) {
                    w <- w + t(repmat(w.mu / w.sd, n, 1))
                }
            }
            basis0 <- basis
            basis <- basis.rescale
        }
        print(range(w))
        nbasis <- dim(w)[1]
        print(paste0("Basis type = ", basis_type, ", with nbasis = ", nbasis))
        x <- x0[-p, ] # remove RewP
        p <- p - 1
    } else {
        curves <- t0 <- w <- NULL
        x <- x0
    }

    res.list <- list(x = x, w = w, y = y, curves = curves, t = t0, basis = basis)
    return(res.list)
}
