source("utilities.R")
source("MCMC_sampling.R")

############ Main Process for Simulation ##########
main_process <- function(x, y, curves, t = NULL, tau = 0.5, registration = TRUE,
                         srvf.lambda = 0, showplot.srvf = FALSE,
                         curves.template, norder = 4, nbasis = 7,
                         K = 10, seed = 2023, n.chains = 2, Nsamples = 1e5,
                         slice.w = 1, slice.m = 1, a.sigma = 1, b.sigma = 1,
                         a.DP = 2, mean.mu = 0, sd.mu = 10, shape.xi = 2, scale.xi = 3,
                         scale.eta = 2, burnin = 5e3, thinning = 300,
                         psrf = FALSE, autocorr = FALSE, trace = FALSE, density = FALSE,
                         effsize = TRUE, max_iters = 5,
                         effsize.fz = FALSE, save.samples = FALSE,
                         r = 1, p.kl = 0.5) {
    ### Descriptions ####
    #### Input ####
    ### x: scalar predictors, (p,n)
    ### y: response, (n)
    ### curves: functional predictors, (N, n, nf)
    ### t: index of the curves, (N)
    ### tau: the quantile level to be examined
    ### registration: bool, if conduct the registration iteratively together with MCMC
    ### srvf.lambda: control the level of penalty in the elastic registration
    ### showplot.srvf: bool, visualize the detailed process of registration
    ### curves.template: the template used in the elastic registration,
    ###         a matrix (N, nf),
    ###         or vector (nf) of integers denoting the index of curves as the template
    ###         or string "median"/"mean" denoting the Karcher Median/Mean as the template
    ### norder: the order of B-Splines for the functional coeffficients
    ### nbasis: the number of basis of B-Splines for the functional coefficients
    ### K: the number of gaussian cdfs used in estimating the link function
    ### seed: for random number generator
    ### n.chains: the number of MCMC chains to obtain the estimators
    ### Nsamples: the number of samples obtained in each chain
    ### slice.w: the window size in the slice sampler, a scalar or a list of window size for each parameter
    ### slice.m: the max number of steps within the slice sampler

    #### Hyperparamters: ####
    ### a.sigma, b.sigma: the prior of sigma
    ### a.DP:
    ### mean.mu, sd.mu: the prior of each mu_j in mixture gaussian cdf
    ### shape.xi, scale.xi: the prior of each xi_j in mixture gaussian cdf
    ### scale.eta: the parameter in the prior of eta_{k,b}

    ### burnin & thinning: number of samples to discard and thinning in each chain
    ### psrf, autocorr, trace, density, effsize, effsize.fz: bool, for chain diagnosis, check R package "coda" for more details
    ### max_iters: the max number of iterations for the whole process
    ### save.samples: bool, whether to save all the samples in each chain, (will always save the effective samples)
    ### r: the constraint (2-norm) of scalar coefficients
    ### p.kl: in (0,1), control the sparsity level of beta

    ### n: sample size
    ### N: number of grids for curves
    ### p: number of scalar predictors
    ### nf: number of functional predictors

    #### Output:####
    # **a list containing the following attributes:**
    # - **res.list0: a list containing the final estimates after convergence**
    #   - res0: a list containing the effective samples and Bayesian estimates for each chain
    #     - all.chains: a list containing the mcmc chains after burnin & thinning for each chain
    #     - all.estimators: a list containing the Bayesian estimates for each chain
    #   - paras.init: a list containing the initial values for each parameter
    #   - sigma.est
    #   - theta.est
    #   - alpha.est
    #   - beta.est
    #   - z.est
    #   - fz.est
    #   - b.est
    #   - mu.est
    #   - xi.est
    # - iters: the number of iterations running (registration + MCMC) procedure
    # - all.alpha.est: a (iters, p) matrix, each row contains $\hat{\boldsymbol{\alpha}}$ after one iteration of (registration + MCMC) procedure
    # - all.beta.est: a (iters, nbasis*nf) matrix, each row contains $\hat{\boldsymbol{\beta}}$ after one iteration of (registration + MCMC) procedure
    # - all.gamma.est: a (N, nf, iters) array, each slice contains $(\hat{\gamma}_1(t), \hat{\gamma}_2(t))$ after one iteration of (registration + MCMC) procedure
    # - all.curves.reg: a list containing the registered curves for each iteration 
    # - all.z.est: a (iters, n) matrix, each row contains $\hat{z}_i, i=1,\ldots,n$ after one iteration of (registration + MCMC) procedure
    # - all.fz.est: a (iters, n) matrix, each row contains $\hat{f}(\hat{z}_i), i=1,\ldots, n$ after one iteration of (registration + MCMC) procedure



    ### Decide if the model contains functional predictors ####
    ## if no functional predictors, directly obtain the estimators from MCMC;
    ## else, conduct the registration & MCMC iteratively until convergency
    if (is.null(curves)) { # no functional predictors
        res.list0 <- MCMC_estimators_spike_slab(x, y, NULL, nbasis, 2, tau,
            r = r, K = K, seed = seed, p.kl = p.kl, lam = 0.5, n.chains = n.chains,
            Nsamples = Nsamples, slice.w = slice.w, slice.m = slice.m,
            a.sigma = a.sigma, b.sigma = b.sigma,
            a.DP = a.DP, mean.mu = mean.mu, sd.mu = sd.mu, shape.xi = shape.xi,
            scale.xi = scale.xi, scale.eta = scale.eta, burnin = burnin,
            thinning = thinning, psrf = psrf, autocorr = autocorr,
            trace = trace, density = density, effsize = effsize,
            effsize.fz = effsize.fz, save.samples = save.samples
        )
        res.list <- list(res.list0 = res.list0)
    } else {
        require(fda)
        n <- length(y)
        N <- dim(curves)[1]
        if (length(dim(curves)) == 3) {
            nf <- dim(curves)[3]
        } else {
            nf <- 1
            curves <- array(curves, dim = c(N, n, nf))
        }
        if (is.null(t)) {
            t <- seq(0, 1, length.out = N)
        }
        if (!registration) { # no registration
            print("Obtain the estimators based on original functional predictors")
            #### Obtain MCMC estimators based on non-registered curves ####
            w <- reconstruct.curves(curves, t, norder = norder, nbasis = nbasis)
            res.list0 <- MCMC_estimators_spike_slab(x, y, w, nbasis, nf, tau,
                r = r, K = K, seed = seed, p.kl = p.kl, lam = 0.5, n.chains = n.chains,
                Nsamples = Nsamples, slice.w = slice.w, slice.m = slice.m,
                a.sigma = a.sigma, b.sigma = b.sigma,
                a.DP = a.DP, mean.mu = mean.mu, sd.mu = sd.mu, shape.xi = shape.xi,
                scale.xi = scale.xi, scale.eta = scale.eta, burnin = burnin,
                thinning = thinning, psrf = psrf, autocorr = autocorr,
                trace = trace, density = density, effsize = effsize,
                effsize.fz = effsize.fz, save.samples = save.samples
            )

            #### construct gamma(t) ####
            beta.est0 <- res.list0$beta.est
            psi <- create.bspline.basis(norder = norder, nbasis = nbasis)
            eval.psi <- eval.basis(t, psi)
            gamma.est0 <- matrix(nrow = N, ncol = nf)
            for (k in 1:nf) {
                idx <- (1:nbasis) + (k - 1) * nbasis
                gamma.est0[, k] <- eval.psi %*% beta.est0[idx]
            }
            res.list <- list(res.list0 = res.list0, gamma.est = gamma.est0)
        } else {
            require(fdasrvf)
            print("Iteratively conduct Registration & MCMC until convergency")
            ### Iteratively conduct Registration & MCMC until convergency ####
            # all.theta1.est <- NULL
            # all.theta2.est <- NULL
            all.alpha.est <- NULL
            all.beta.est <- NULL
            all.z.est <- NULL
            all.fz.est <- NULL
            all.gamma.est <- array(dim = c(N, nf, max_iters))
            all.curves.reg <- list()
            alpha.diff.norm <- beta.diff.norm <- 1e4
            iters <- 1
            paras.init <- NULL

            # while ((alpha.diff.norm > 0.1) && (beta.diff.norm > 1) && (iters <= max_iters)) {
            while ((beta.diff.norm > 1) && (iters <= max_iters)) {
                print(paste("iter", iters))
                #### (1) registration based on the template ####
                print("Registration based on the template")
                if (!is.null(curves.template)) {
                    if (is.integer(curves.template)) {
                        ## the index of the curve to be used as the template
                        curves.template <- curves[, curves.template, ] # (N, nf)
                    }
                    curves.reg0 <- array(dim = c(N, n, nf))
                    for (k in 1:nf) {
                        srvf.res <- multiple_align_functions(curves[, , k], t,
                            curves.template[, k],
                            lambda = srvf.lambda, showplot = showplot.srvf
                        )
                        curves.reg0[, , k] <- srvf.res$fn
                    }
                } else { # no template, use Karcher Median as the template,
                    # but it's very time-consuming
                    # instead, no registration in the first step
                    curves.reg0 <- curves
                }
                all.curves.reg[[iters]] <- curves.reg0

                #### (2) construct the new predictors by chosen B-Splines ####
                print("Construct the new predictors")
                w.reg0 <- reconstruct.curves(curves.reg0, t, norder = norder, nbasis = nbasis)
                print(dim(w.reg0))
                #### (3) Obtain MCMC estimators based on registered curves ####
                print("Obtain MCMC estimators")
                res.list0 <- MCMC_estimators_spike_slab(x, y, w.reg0, nbasis, nf, tau,
                    r = r, K = K, seed = seed, p.kl = p.kl, lam = 0.5, n.chains = n.chains,
                    Nsamples = Nsamples, slice.w = slice.w, slice.m = slice.m,
                    a.sigma = a.sigma, b.sigma = b.sigma,
                    a.DP = a.DP, mean.mu = mean.mu, sd.mu = sd.mu, shape.xi = shape.xi,
                    scale.xi = scale.xi, scale.eta = scale.eta, burnin = burnin,
                    thinning = thinning, psrf = psrf, autocorr = autocorr,
                    trace = trace, density = density, effsize = effsize,
                    effsize.fz = effsize.fz, save.samples = save.samples,
                    paras.init = paras.init
                )
                alpha.est0 <- res.list0$alpha.est
                beta.est0 <- res.list0$beta.est
                z.est0 <- res.list0$z.est
                fz.est0 <- res.list0$fz.est

                #### (4) check the change of estimators d ####
                if (iters > 1) {
                    alpha.diff.norm <- norm(alpha.est0 - alpha.est1, "2")
                    beta.diff.norm <- norm(beta.est0 - beta.est1, "2")
                    pracma::fprintf(
                        "iter = %d, alpha.diff.norm = %.3G, beta.diff.norm = %.3G\n",
                        iters, alpha.diff.norm, beta.diff.norm
                    )
                }
                all.alpha.est <- rbind(all.alpha.est, alpha.est0)
                all.beta.est <- rbind(all.beta.est, beta.est0)
                all.z.est <- rbind(all.z.est, z.est0)
                all.fz.est <- rbind(all.fz.est, fz.est0)

                #### (5) construct gamma(t) & use it as the new template ####
                psi <- create.bspline.basis(rangeval = range(t), norder = norder, nbasis = nbasis)
                eval.psi <- eval.basis(t, psi)
                gamma.est0 <- matrix(nrow = N, ncol = nf)
                for (k in 1:nf) {
                    idx <- (1:nbasis) + (k - 1) * nbasis
                    gamma.est0[, k] <- eval.psi %*% beta.est0[idx]
                }
                all.gamma.est[, , iters] <- gamma.est0
                curves.template0 <- curves.template
                curves.template <- gamma.est0

                #### (6) update the other parameters to check the convergence ####
                alpha.est1 <- alpha.est0
                beta.est1 <- beta.est0
                iters <- iters + 1

                #### update the initial values in MCMC ####
                paras.init <- res.list0$paras.init
                # paras.init$theta.t <- res.list0$theta.est
                paras.init$theta.t <- alpha.transform.inverse(alpha.est0)
                paras.init$alpha.t <- alpha.est0
                paras.init$beta.t <- beta.est0
                paras.init$sigma.t <- res.list0$sigma.est
                paras.init$b.t <- res.list0$b.est
                paras.init$pi.t <- pi.transform(paras.init$b.t)
                paras.init$mu.t <- res.list0$mu.est
                paras.init$xi.t <- res.list0$xi.est
            }

            ### Obtain the final registered curves & estimators ####
            res.list <- list(
                res.list0 = res.list0, # all.theta1.est = all.theta1.est,
                # all.theta2.est = all.theta2.est,
                all.alpha.est = all.alpha.est,
                all.beta.est = all.beta.est, all.gamma.est = all.gamma.est,
                all.curves.reg = all.curves.reg, iters = iters,
                curves.template = curves.template,
                curves.template0 = curves.template0,
                all.z.est = all.z.est, all.fz.est = all.fz.est
            )
        }
        print("Finished Implementing the Framework!")
    }
    return(res.list)
}
