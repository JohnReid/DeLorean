functions {
    #
    # Squared exponential covariance function
    #
    real
    se_cov(real r) {
        return exp(- pow(r, 2) / 2);
    }
    #
    # Matern nu=3/2 covariance function
    #
    real
    matern32_cov(real r) {
        real x;
        x = sqrt(3.) * fabs(r);
        return (1 + x) * exp(- x);
    }
    #
    # Matern nu=5/2 covariance function
    #
    real
    matern52_cov(real r) {
        real x;
        x = sqrt(5.) * fabs(r);
        return (1 + x + pow(x, 2) / 3) * exp(- x);
    }
    #
    # Covariance function
    #
    real
    cov_fn(real r) {
      return matern32_cov(r);
      # return se_cov(r);
      # return matern52_cov(r);
    }
    #
    # Calculate symmetric covariance matrix
    #
    matrix
    cov_symmetric(matrix latent, vector lengthscales) {
      matrix[cols(latent),cols(latent)] result;
      for (c1 in 1:cols(latent)) {
        for (c2 in 1:c1) {
          real r;
          if (c1 == c2) {
            r = 0;
          } else {
            r = distance(col(latent, c1) ./ lengthscales, col(latent, c2) ./ lengthscales);
          }
          result[c1,c2] = cov_fn(r);
          if (c1 != c2) {
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
    cov(matrix latent1, matrix latent2, vector lengthscales) {
      matrix[cols(latent1),cols(latent2)] result;
      for (c1 in 1:cols(latent1)) {
        for (c2 in 1:cols(latent2)) {
          result[c1,c2] = cov_fn(distance(col(latent1, c1) ./ lengthscales, col(latent2, c2) ./ lengthscales));
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
    #
    # Data
    #
    # Time data
    row_vector<lower=0>[C] time;  # Time index for cell c
    vector[C] expr[G];            # Expression data
    row_vector[G] phi;            # Gene mean
    #
    # Hyperparameters
    #
    real mu_psi;                  # Mean of log between time variation, psi
    real<lower=0> sigma_psi;      # S.d. of log between time variation, psi
    real mu_omega;                # Mean of log within time variation, omega
    real<lower=0> sigma_omega;    # S.d. of log within time variation, omega
    real<lower=0> sigma_tau;      # Standard deviation for z
    vector[2] lengthscales;       # Length scales for pseudotime and z
}
transformed data {
    cov_matrix[C] identity;  # Unit diagonal covariance matrix
    identity = diag_matrix(rep_vector(1, C));
}
parameters {
    row_vector[C] tauoffset;       # Pseudotime offset
    row_vector[C] z;               # Latent dimension
    row_vector<lower=0>[G] psi;    # Between time variance
    row_vector<lower=0>[G] omega;  # Within time variance
}
transformed parameters {
    row_vector[C] tau;             # Pseudotime
    tau = time + tauoffset;
}
model {
    matrix[2, C] latent;           # Latent variables
    latent[1] = tau;
    latent[2] = z;
    #
    # Sample gene-specific temporal variance and noise
    target += lognormal_lpdf(psi|mu_psi, sigma_psi);
    target += lognormal_lpdf(omega|mu_omega, sigma_omega);
    #
    # Sample latent variables
    target += normal_lpdf(tauoffset|0, sigma_tau);  # Pseudotime offsets
    target += normal_lpdf(z|0, 1);                  # Z
    #
    # Expression values for each gene
    for (g in 1:G) {
        target += multi_normal_lpdf(expr[g]|
                    rep_vector(phi[g], C),
                    psi[g] * cov_symmetric(latent, lengthscales)
                        + omega[g] * identity);
    }
}
