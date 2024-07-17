#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double logLikelihood_cpp(rowvec alpha, rowvec beta, double sigma, vec e, 
        bool link, int fz_idx, double z_mu, double z_sigma, vec pi, vec mu, 
        vec xi, vec y, mat preds_all, int k1, int k2){
            int n = y.n_rows; 
            int p = alpha.n_cols;
            int p_all = preds_all.n_rows;
            arma::mat x = preds_all.rows(0,p-1); //(p, n)
            arma::rowvec z = alpha * x;
            arma::vec f_z(n);
            double logL; 
            if (p < p_all) {
                arma::mat w = preds_all.rows(p, p_all - 1); // (L, n)
                z = alpha * x + beta * w; //(1, n)
            } else {
                z = alpha * x;
                }
            if (!link){
                for (int i=0; i < n; i++){
                    arma::vec pdfs = arma::normcdf(z(i), mu, sqrt(xi));
                    f_z(i) = sum(pi % pdfs);
                }
            } else {
                if (fz_idx == 0){
                    f_z = 1 - exp(-(z.st() + z_mu) * z_sigma);
                } else if (fz_idx == 1){
                    f_z = arma::normcdf(z.st(), z_mu, z_sigma);
                }
            }
            arma::vec kernel = square(f_z) / e + 2 * (k1 - y / e) % f_z;
            logL = - sum(kernel) / (2 * k2 * sigma);
            return logL;
        }