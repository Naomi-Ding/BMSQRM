############ Run Multiple Simulations on Server ##########
### Input: simulation settings ####
args <- commandArgs(TRUE)
n <- (eval(parse(text = args[[1]]))) # sample size, n=200
N <- (eval(parse(text = args[[2]]))) # number of grids
tau <- (eval(parse(text = args[[3]]))) # quantile level, 0.25, 0.5, 0.75
func.p <- (eval(parse(text = args[[4]]))) # true or false
fz.idx <- (eval(parse(text = args[[5]]))) # choice of link function, 0: exp, 1: normal, 2: laplace
K <- (eval(parse(text = args[[6]]))) # number of mixture gaussian cdf
n.chains <- (eval(parse(text = args[[7]]))) # number of chains
Nsamples <- (eval(parse(text = args[[8]]))) # number of samples from each chain
burnin <- (eval(parse(text = args[[9]])))
thinning <- (eval(parse(text = args[[10]])))
max_iters <- (eval(parse(text = args[[11]]))) # max number of iterations in (registration, MCMC)
jobid <- (eval(parse(text = args[[12]])))
# range from 1-200, index for replications, used as the seed in data generation
seed <- 2023 + jobid
save.samples <- (eval(parse(text = args[[13]])))
visualize.results <- (eval(parse(text = args[[14]])))

if (func.p) {
    nbasis <- (eval(parse(text = args[[15]])))
    print(paste0("nbasis for gamma_k(t) = ", nbasis))
    nf <- (eval(parse(text = args[[16]])))
    p.kl <- (eval(parse(text = args[[17]])))
    print(paste0("p.kl for sparsity = ", p.kl))

    registration <- (eval(parse(text = args[[18]]))) # registration or not
    if (registration) {
        srvf.lambda <- (eval(parse(text = args[[19]])))
        print(paste0("Registration on the curves: ", registration, " using lambda = ", srvf.lambda))
    }
} else {
    registration <- FALSE
}

spike <- TRUE
beta.constrain <- FALSE
beta.sparse <- TRUE
curve_center <- FALSE
slice.m <- 1

print(paste0(
    "Sample size=", n, ", number of grids=", N, ", tau=", tau,
    ", functional predictors=", func.p,
    ", number of functional predictors =", nf,
    ", choice of link function=", fz.idx,
    ", number of mixture gaussian=", K, ", number of chains=", n.chains,
    ", number of samples in each chain=", Nsamples, ", burnin=", burnin,
    ", thinning=", thinning,
    ", max number of iterations for (registration, MCMC)=", max_iters,
    ", jobid=", jobid, ", seed=", seed, ", save samples=", save.samples
))

### Settings ####
scale.eta <- 2
theta2.case <- 2
p <- 4
r <- 1
r2 <- 10
norder <- 4
# nbasis <- 7
umin <- 0.5
umax <- 1.5
z0 <- 0
if (func.p) { # contain functional predictors
    if (nf == 1) {
        if (theta2.case == 1) {
            if (fz.idx == 2) {
                z_mu <- 8
                z_sigma <- 1
            } else if (fz.idx == 3) {
                z_mu <- 7
                z_sigma <- 5
            }
        } else if (theta2.case == 2) {
            if (fz.idx == 3) { # pareto
                z_mu <- 3
                z_sigma <- 4
            } else if (fz.idx == 2) { # laplace
                z_mu <- 4
                z_sigma <- 0.8
            }
        } else if (theta2.case == 3) {
            if (fz.idx == 2) { # laplace
                z_mu <- 11
                z_sigma <- 1.2
            } else if (fz.idx == 3) {
                z_mu <- 9
                z_sigma <- 5
            }
        }
    } else if (nf == 2) {
        if (beta.constrain) {
            if (fz.idx == 0) { # for exp
                z_mu <- 6
                z_sigma <- 0.2
            } else if (fz.idx == 3) { # for pareto
                z_mu <- 9
                z_sigma <- 5
            } else { # for normal & laplace
                z_mu <- 10.9
                z_sigma <- 1
            }
        } else {
            if (fz.idx == 1) {
                z_mu <- 1.6
                z_sigma <- 0.9
            } else if (fz.idx == 2) { # laplace cdf
                if (nbasis == 7) {
                    z_mu <- 2
                } else if (nbasis == 9) {
                    z_mu <- 1.6
                }
                z_sigma <- 0.6
            } else if (fz.idx == 3) {
                if (nbasis == 9) {
                    z_mu <- 1
                    z_sigma <- 3
                } else if (nbasis == 7) {
                    z_mu <- 1.2
                    z_sigma <- 3
                }
            }
        }
    }
} else { # only scalar predictors
    if (fz.idx == 0) { # for exp
        z_mu <- -2
        z_sigma <- 0.8
    } else if (fz.idx == 3) { # for pareto
        z0 <- 2
        z_mu <- 1.2
        z_sigma <- 2.5
    } else { # for normal & laplace
        z_mu <- 0
        z_sigma <- 0.5
    }
}

#### hyperparamters ####
a.sigma <- 1
b.sigma <- 1
a.DP <- 2
mean.mu <- 0
sd.mu <- 10
shape.xi <- 2
scale.xi <- 3
slice.w <- 1
psrf <- FALSE

showplot.srvf <- FALSE
curves.visualize <- FALSE
response.plot <- FALSE
bss.sd1 <- bss.sd2 <- 0.5


### STEP 0. set up path for saving results ####
results.name.full <- paste0(
    "func=", func.p, "_regist=", registration,
    "_nf=", nf, "_nbasis=", nbasis,
    "_pkl=", p.kl, "_case=", fz.idx, "_n=", n, "_N=", N, "_tau=", tau, ".rds"
)
print(results.name.full)


### STEP 1. Data generation #####
print("Simulating Data...")
start <- Sys.time()
source("simu_data_gen.R")
data <- simu.data.gen(
    n, N, p, r, r2, func.p, tau, fz.idx,
    nf, norder, nbasis, seed, umin, umax, z_mu, z_sigma, z0,
    srvf.lambda = srvf.lambda, showplot.srvf,
    curves.visualize, response.plot,
    bss.sd1 = bss.sd1, bss.sd2 = bss.sd2,
    beta.constrain = beta.constrain, beta.sparse = beta.sparse,
    curve_center = curve_center
)
x <- data$x
if (func.p) {
    curves <- data$curves
    t <- data$t
} else {
    curves <- t <- NULL
}
z <- data$z
fz <- data$fz
Q <- data$Q
y <- data$y
theta1 <- data$theta1
alpha <- data$alpha
gamma <- data$gamma
theta2 <- data$theta2
beta <- data$beta

#### check the range of y ####
print(paste0("z_mu = ", z_mu, ", z_sigma = ", z_sigma), digits = 4)
fprintf("Range of z: (%.3G, %.3G)\n", range(z)[1], range(z)[2])
fprintf("Range of Q: (%.3G, %.3G)\n", range(Q)[1], range(Q)[2])
fprintf("Range of y: (%.3G, %.3G)\n", range(y)[1], range(y)[2])
if (min(y) < 0 || max(y) > 1) {
    print("Improper range of y! Adjust the link function!")
} else {
    ### STEP2. CONDUCT OVERALL PROCEDURE ####
    print("Implementing the Framework...")
    source("main_process.R")
    res.list <- main_process(x, y, curves, t, tau, registration,
        srvf.lambda, showplot.srvf, gamma, norder, nbasis,
        K = K, seed = seed, n.chains = n.chains, Nsamples = Nsamples,
        slice.w = slice.w, slice.m = slice.m, a.sigma = a.sigma, b.sigma = b.sigma,
        a.DP = a.DP, mean.mu = mean.mu, sd.mu = sd.mu, shape.xi = shape.xi,
        scale.xi = scale.xi, , scale.eta = scale.eta,
        burnin = burnin, thinning = thinning, psrf = psrf,
        autocorr = FALSE, trace = FALSE, density = FALSE, effsize = TRUE,
        max_iters = max_iters, est.fz = TRUE, effsize.fz = FALSE,
        save.samples = save.samples,
        r = r, p.kl = p.kl
    )
    end <- Sys.time()
    print(paste0("Total Elapsed Time: ", end - start))

    ### STEP3. Performance Evaluation ####
    print("Performance Evaluation")
    res.list0 <- res.list$res.list0
    alpha.est0 <- res.list0$alpha.est
    sse.alpha <- sum((alpha.est0 - alpha)^2)
    errors <- c(sse.alpha)

    if (func.p) {
        require(fda)
        beta.est0 <- res.list0$beta.est
        psi <- create.bspline.basis(norder = norder, nbasis = nbasis)
        eval.psi <- eval.basis(t, psi)
        gamma.est0 <- matrix(nrow = N, ncol = nf)
        for (k in 1:nf) {
            idx <- (1:nbasis) + (k - 1) * nbasis
            gamma.est0[, k] <- eval.psi %*% beta.est0[idx]
        }
        # iters <- res.list$iters - 1
        # gamma.est0 <- res.list$all.gamma.est[, , iters]
        sse.beta <- sum((beta.est0 - beta)^2)
        ise.gamma <- colMeans((gamma.est0 - gamma)^2)
        errors <- c(errors, sse.beta, ise.gamma)
        errors.name <- c("sse.alpha", "sse.beta", paste0("ise.gamma", 1:nf), "mse.z", "mse.f")
    } else {
        errors.name <- c("sse.alpha", "mse.z", "mse.f")
    }

    z.est0 <- res.list0$z.est
    fz.est0 <- res.list0$fz.est
    mse.z <- mean((z.est0 - as.numeric(z))^2)
    mse.f <- mean((fz.est0 - Q)^2)
    errors <- c(errors, mse.z, mse.f)
    names(errors) <- errors.name
    print(errors)


    ### STEP 4. Visualization of gamma(t) & f(z) ####
    if (visualize.results) {
        require(ggplot2)
        print("Visualization of gamma(t) & f(z)")
        if (func.p) {
            ##### (1) gamma(t)
            q.gamma <- ggplot(data = data.frame(
                t = rep(t, 2), gamma = c(gamma),
                gamma.hat = c(gamma.est0),
                group = factor(rep(1:2, each = N))
            )) +
                geom_line(aes(x = t, y = gamma, group = group), colour = "black", linewidth = 1.1) +
                geom_line(aes(x = t, y = gamma.hat, group = group), colour = "red", linewidth = 1.1) +
                ggtitle(expression(gamma(t))) +
                ylab("") +
                facet_wrap(~group, scales = "free")
            ggsave(sprintf("fz%d_tau%.2f_gamma.png", fz, tau), width = 5, height = 5)
        }

        ##### (2) f(z)
        require(KernSmooth)
        df.q <- cbind(z, Q)[order(z), ]
        fz.lp <- KernSmooth::locpoly(z.est0, fz.est0, bandwidth = 0.1)
        qz <- ggplot(data = data.frame(z = z, Q = Q, y = y)) +
            geom_point(aes(x = z, y = y), color = "gray50", show.legend = T) +
            ylim(0, 1) +
            ggtitle(expression(f(z))) +
            xlab("z") +
            geom_line(aes(z, Q),
                data = data.frame(df.q),
                linewidth = 1.1, show.legend = TRUE
            ) +
            geom_line(aes(x, y),
                data = data.frame(fz.lp), linewidth = 1.1,
                color = "red"
            )
        ggsave(sprintf("fz%d_tau%.2f_fz.png", fz, tau), width = 8, height = 6)
    }

    ### STEP 4. Save the results ####
    print("Saving the Results...")
    results <- list(args = args, data = data, res.list = res.list, errors = errors)
    saveRDS(results, file = results.name.full)
}
