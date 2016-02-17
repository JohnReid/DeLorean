functions {
    #
    # Periodic function
    #
    real
    periodise(real r, real period) {
        return period / 2 * sin(r * pi() / period);
    }
    #
    # Squared exponential covariance function
    #
    real
    se_cov(real r, real l) {
        return exp(- pow(r / l, 2) / 2);
    }
    #
    # Matern nu=3/2 covariance function
    #
    real
    matern32_cov(real r, real l) {
        real x;
        x <- sqrt(3) * fabs(r / l);
        return (1 + x) * exp(- x);
    }
    #
    # Matern nu=5/2 covariance function
    #
    real
    matern52_cov(real r, real l) {
        real x;
        x <- sqrt(5) * fabs(r / l);
        return (1 + x + pow(x, 2) / 3) * exp(- x);
    }
    #
    # Covariance function
    #
    real
    cov_fn(real r, int periodic, real period, real l) {
        real rprime;
        if (periodic) {
            rprime <- periodise(r, period);
        } else {
            rprime <- r;
        }
        return matern32_cov(rprime, l);
        # return se_cov(rprime, l);
        # return se_cov(rprime, l);
        # return matern52_cov(rprime, l);
        # return matern52_cov(rprime, l);
    }
    #
    # Calculate symmetric covariance matrix
    #
    matrix
    cov_symmetric(row_vector tau, int periodic, real period, real l) {
        matrix[cols(tau),cols(tau)] result;
        for (c1 in 1:cols(tau)) {
            for (c2 in 1:c1) {
                result[c1,c2] <- cov_fn(tau[c2] - tau[c1], periodic,
                                        period, l);
                if(c1 != c2) {
                    result[c2,c1] <- result[c1,c2];
                }
            }
        }
        return result;
    }
    #
    # Calculate covariance matrix
    #
    matrix
    cov(row_vector tau1,
        row_vector tau2,
        int periodic,
        real period,
        real l
    ) {
        matrix[cols(tau1),cols(tau2)] result;
        for (c1 in 1:cols(tau1)) {
            for (c2 in 1:cols(tau2)) {
                result[c1,c2] <- cov_fn(tau2[c2] - tau1[c1], periodic, period, l);
            }
        }
        return result;
    }
    #
    # Calculate approximation to precision matrix
    #
    matrix
    calc_approx_prec(
        matrix Autau,
        real psi_g,
        real omega_g
    ) {
        vector[cols(Autau)] Binv;
        vector[cols(Autau)] Binvsqrt;
        Binv <- psi_g * (diagonal(crossprod(Autau) - 1) + omega_g);
        for(c in 1:cols(Autau)) {
            Binvsqrt[c] <- sqrt(Binv[c]);
        }
        Binv <- psi_g * (diagonal(crossprod(Autau) - 1) + omega_g);
        return
            diag_matrix(Binv)
            - psi_g * tcrossprod(mdivide_right_tri_low(
                diag_pre_multiply(Binv, Autau'),
                cholesky_decompose(
                    diag_matrix(rep_vector(1, rows(Autau)))
                    + psi_g*tcrossprod(diag_post_multiply(Autau,
                                                          Binvsqrt)))));
    }
    void
    pretty_print_tri_lower(matrix x) {
        if (rows(x) == 0) {
            print("empty matrix");
            return;
        }
        print("rows=", rows(x), " cols=", cols(x));
        for (m in 1:rows(x))
            for (n in 1:m)
                print("[", m, ",", n, "]=", x[m,n]);
    }
}
data {
    #
    # Dimensions
    #
    int<lower=2> C;  # Number of cells
    int<lower=2> G;  # Number of genes
    int<lower=0> H;  # Number of held out genes
    int<lower=0> M;  # Number of inducing pseudotimes
    #
    # Data
    #
    # Time data
    int periodic; # Are the expression patterns periodic?
    real period;  # Cyclic period
    row_vector<lower=0>[C] time;  # Time index for cell c
    # Expression data
    vector[C] expr[G+H];
    row_vector[G+H] phi;  # gene mean
    #
    # Inducing pseudotimes
    #
    row_vector[M] u; # pseudotimes
    #
    # Hyperparameters
    #
    real mu_psi;  # Mean of log between time variation, psi
    real<lower=0> sigma_psi;  # S.d. of log between time variation, psi
    real mu_omega;  # Mean of log within time variation, omega
    real<lower=0> sigma_omega;  # S.d. of log within time variation, omega
    real l;  # Length scale
    real<lower=0> sigma_tau;  # Standard deviation for pseudotime
    #
    # Held out parameters
    #
    row_vector<lower=0>[H] heldout_psi;    # Between time variance
    row_vector<lower=0>[H] heldout_omega;  # Within time variance
    #
    # Generated quantities
    #
    int numtest; # Number of test inputs for predictions
    row_vector[numtest] testinput; # Test inputs for predictions
}
transformed data {
    vector[C] expradj[G+H];  # Mean adjusted expression
    cholesky_factor_cov[M] KuuChol; # pseudotime covariance Cholesky
    KuuChol <- cholesky_decompose(cov_symmetric(u, periodic, period, l));
    # Calculate adjusted genes
    for (g in 1:(G+H)) {
        expradj[g] <- expr[g] - phi[g];
    }
}
parameters {
    row_vector[C] tau;    # Pseudotime
    row_vector<lower=0>[G] psi;    # Between time variance
    row_vector<lower=0>[G] omega;  # Within time variance
}
model {
    #
    # Approximation to pseudotime covariance matrix
    matrix[M,C] Autau;  # Sqrt of approximation
    vector[C] Qtautaudiag;   # Approximation to covariance
    Autau <- mdivide_left_tri_low(KuuChol, cov(u, tau, periodic, period, l));
    Qtautaudiag <- columns_dot_self(Autau)';
    # Check that diag(Ktautau - Qtautau) is positive
    for (c in 1:C) {
        if (Qtautaudiag[c] > 1) {
            reject("Qtautaudiag must be less than 1.", Qtautaudiag[c]);
        }
    }
    #
    # Sample gene-specific factors
    psi ~ lognormal(mu_psi, sigma_psi);  # Temporal variation
    omega ~ lognormal(mu_omega, sigma_omega);  # Noise
    #
    # Sample pseudotime
    tau ~ normal(time, sigma_tau);  # Pseudotime
    #
    # Expression values for each gene
    for (g in 1:G) {
        vector[C] Binv;
        vector[C] Binvy;
        vector[M] b;  # Autau * B-1 * y
        matrix[M,M] Vchol;  # Low dimension decomposition
        real det_cov;  # The determinant of the covariance
        #
        # Approximation to pseudotime cross-covariance matrix
        #
        # Sample covariance parameters
        psi ~ lognormal(mu_psi, sigma_psi);  # Temporal variation
        omega ~ lognormal(mu_omega, sigma_omega);  # Noise
        #
        # Inverse of high dimensional covariance diagonal
        # Here we are assuming that the diagonal of Ktautau is 1
        Binv <- 1 ./ (omega[g] + psi[g] * (1 - Qtautaudiag));
        Binvy <- Binv .* expradj[g];
        #
        # Invert low dimensional matrix
        Vchol <- cholesky_decompose(
            diag_matrix(rep_vector(psi[g], M))
            + diag_post_multiply(Autau, Binv) * Autau');
        #
        # Calculate term in quadratic form part of likelihood
        b <- mdivide_left_tri_low(Vchol, Autau * Binvy);
        #
        # Calculate determinant of the covariance
        det_cov <- square(prod(diagonal(Vchol))) / prod(Binv);
        #
        # Increment log probability with multivariate normal log likelihood
        increment_log_prob(-.5 * (log(det_cov)
                                    + dot_product(expradj[g], Binvy)
                                    - psi[g] * dot_self(b)));
    }
}
generated quantities {
}
