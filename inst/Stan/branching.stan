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
        x = sqrt(3.) * fabs(r / l);
        return (1 + x) * exp(- x);
    }
    #
    # Matern nu=5/2 covariance function
    #
    real
    matern52_cov(real r, real l) {
        real x;
        x = sqrt(5.) * fabs(r / l);
        return (1 + x + pow(x, 2) / 3) * exp(- x);
    }
    #
    # Covariance function
    #
    real
    cov_fn(real rtau, real rz, int periodic, real period, real l) {
      real rprime;
      real r;
      # Check if we need to make periodic
      if (periodic) {
        rprime = periodise(rtau, period);
      } else {
        rprime = rtau;
      }
      # Length scale only applies to pseudotime dimension
      # Need fabs inside sqrt to make gradient calculation behave
      r = sqrt(fabs(square(rprime / l) + square(rz)));
      return matern32_cov(r, 1);
      # return se_cov(r, 1);
      # return se_cov(r, 1);
      # return matern52_cov(r, 1);
      # return matern52_cov(r, 1);
    }
    #
    # Calculate symmetric covariance matrix
    #
    matrix
    cov_symmetric(row_vector tau, row_vector z, int periodic, real period, real l) {
      matrix[cols(tau),cols(tau)] result;
      for (c1 in 1:cols(tau)) {
        for (c2 in 1:c1) {
          result[c1,c2] = cov_fn(tau[c2] - tau[c1],
                                 z[c2] - z[c1],
                                 periodic,
                                 period, l);
          if(c1 != c2) {
            result[c2,c1] = result[c1,c2];
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
        # row_vector x1,  # TODO: how to integrate latent dimension into test input prediction?
        # row_vector x2,
        int periodic,
        real period,
        real l
    ) {
      matrix[cols(tau1),cols(tau2)] result;
      for (c1 in 1:cols(tau1)) {
        for (c2 in 1:cols(tau2)) {
          # TODO: fix zero that ignores latent dimension
          result[c1,c2] = cov_fn(tau2[c2] - tau1[c1], 0, periodic, period, l);
        }
      }
      return result;
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
    # Hyperparameters
    #
    real mu_psi;  # Mean of log between time variation, psi
    real<lower=0> sigma_psi;  # S.d. of log between time variation, psi
    real mu_omega;  # Mean of log within time variation, omega
    real<lower=0> sigma_omega;  # S.d. of log within time variation, omega
    real l;  # Length scale
    real<lower=0> sigma_tau;  # Standard deviation for pseudotime
}
transformed data {
    #
    # Unit diagonal covariance matrix
    cov_matrix[C] identity;
    #
    # Transformations of expression
    #
    # Unit diagonal covariance matrix
    identity = diag_matrix(rep_vector(1, C));
}
parameters {
    row_vector[C] tau;    # Pseudotime
    row_vector[C] z;      # Unconstrained latent dimension
    row_vector<lower=0>[G] psi;    # Between time variance
    row_vector<lower=0>[G] omega;  # Within time variance
}
model {
    #
    # Sample gene-specific factors
    target += lognormal_lpdf(psi|mu_psi, sigma_psi);
    target += lognormal_lpdf(omega|mu_omega, sigma_omega);
    #
    # Sample pseudotime
    target += normal_lpdf(tau|time, sigma_tau);  # Pseudotime
    #
    # Sample latent z
    target += normal_lpdf(z|0, 1);  # z
    #
    # Expression values for each gene
    for (g in 1:G) {
        target += multi_normal_lpdf(expr[g]|
                    rep_vector(phi[g], C),
                    psi[g] * cov_symmetric(tau, z, periodic, period, l)
                        + omega[g] * identity);
    }
}
