functions {
    //
    // Squared exponential covariance function
    //
    real
    se_cov(real r) {
        return exp(- pow(r, 2) / 2);
    }
    //
    // Matern nu=3/2 covariance function
    //
    real
    matern32_cov(real r) {
        real x;
        x = sqrt(3.) * fabs(r);
        return (1 + x) * exp(- x);
    }
    //
    // Matern nu=5/2 covariance function
    //
    real
    matern52_cov(real r) {
        real x;
        x = sqrt(5.) * fabs(r);
        return (1 + x + pow(x, 2) / 3) * exp(- x);
    }
    //
    // Covariance function
    //
    real
    cov_fn(real r) {
      return matern32_cov(r);
      // return se_cov(r);
      // return matern52_cov(r);
    }
    //
    // Calculate symmetric covariance matrix
    //
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
    //
    // Calculate covariance matrix
    //
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
    //
    // Calculate approximation to precision matrix
    //
    matrix
    calc_approx_prec(
        matrix Aut,
        real psi_g,
        real omega_g
    ) {
        vector[cols(Aut)] Binv;
        vector[cols(Aut)] Binvsqrt;
        Binv = psi_g * (diagonal(crossprod(Aut) - 1) + omega_g);
        for(c in 1:cols(Aut)) {
            Binvsqrt[c] = sqrt(Binv[c]);
        }
        Binv = psi_g * (diagonal(crossprod(Aut) - 1) + omega_g);
        return
            diag_matrix(Binv)
            - psi_g * tcrossprod(mdivide_right_tri_low(
                diag_pre_multiply(Binv, Aut'),
                cholesky_decompose(
                    diag_matrix(rep_vector(1, rows(Aut)))
                    + psi_g*tcrossprod(diag_post_multiply(Aut,
                                                          Binvsqrt)))));
    }
    //
    // Calculate diagonal of Q
    //
    vector
    calc_Q_diag(matrix KuuChol, matrix Kux) {
      return columns_dot_self(mdivide_left_tri_low(KuuChol, Kux))';
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
    //
    // Dimensions
    //
    int<lower=2> C;               // Number of cells
    int<lower=2> G;               // Number of genes
    int<lower=0> M;               // Number of inducing pseudo-inputs

    //
    // Data
    //
    // Time data
    row_vector<lower=0>[C] time;  // Time index for cell c
    vector[C] expr[G];            // Expression data
    row_vector[G] phi;            // Gene mean
    //
    // Inducing pseudo-inputs
    //
    matrix[2, M] u;               // pseudo-inputs

    //
    // Hyperparameters
    //
    real mu_psi;                  // Mean of log between time variation, psi
    real<lower=0> sigma_psi;      // S.d. of log between time variation, psi
    real mu_omega;                // Mean of log within time variation, omega
    real<lower=0> sigma_omega;    // S.d. of log within time variation, omega
    real<lower=0> sigma_tau;      // Standard deviation for z
    vector[2] lengthscales;       // Length scales for pseudotime and z
}
transformed data {
    cov_matrix[C] identity;         // Unit diagonal covariance matrix
    vector[C] expradj[G];           // Mean adjusted expression
    matrix[M,M] Kuu;                // Covariance between inducing points
    cholesky_factor_cov[M] KuuChol; // pseudotime covariance Cholesky
    identity = diag_matrix(rep_vector(1, C));
    Kuu = cov_symmetric(u, lengthscales);
    KuuChol = cholesky_decompose(Kuu);
    for (g in 1:G) {
      expradj[g] = expr[g] - phi[g];
    }
}
parameters {
    row_vector[C] tauoffset;       // Pseudotime offset
    row_vector[C] z;               // Latent dimension
    row_vector<lower=0>[G] psi;    // Between time variance
    row_vector<lower=0>[G] omega;  // Within time variance
}
transformed parameters {
    row_vector[C] tau;             // Pseudotime
    tau = time + tauoffset;
}
model {
    matrix[2, C] latent;           // Latent variables
    //
    // Approximation to pseudotime covariance matrix
    matrix[M,C] Aut;  // Sqrt of approximation
    vector[C] KminusQdiag;   // Diagonal adjustment to approximate covariance
    //
    // Make a matrix of the pseudotime and latent variables
    latent[1] = tau;
    latent[2] = z;
    //
    // Calculate the square root of Qtautau. Complexity: O(CM^2)
    Aut = mdivide_left_tri_low(KuuChol, cov(u, latent, lengthscales));
    KminusQdiag = 1 - columns_dot_self(Aut)';
    // Check that diag(Ktautau - Qtt) is positive
    for (c in 1:C) {
        if (KminusQdiag[c] < 0) {
            reject("KminusQdiag must be greater or equal to 0. : ",
                    KminusQdiag[c]);
        }
    }
    //
    // Sample gene-specific temporal variance and noise
    target += lognormal_lpdf(psi|mu_psi, sigma_psi);
    target += lognormal_lpdf(omega|mu_omega, sigma_omega);
    //
    // Sample latent variables
    target += normal_lpdf(tauoffset|0, sigma_tau);  // Pseudotime offsets
    target += normal_lpdf(z|0, 1);                  // Z
    //
    // Expression values for each gene
    for (g in 1:G) {
      vector[C] Binv;
      vector[C] Binvy;
      vector[M] b;           // Aut * B-1 * y
      matrix[M,M] Vchol;     // Low dimension decomposition
      real log_det_cov;      // The logarithm of the determinant of the covariance
      //
      // Inverse of high dimensional covariance diagonal
      // Here we are assuming that the diagonal of Ktautau is 1
      Binv = 1 ./ (omega[g] + psi[g] * KminusQdiag);
      Binvy = Binv .* expradj[g];
      //
      // Invert low dimensional matrix
      // Complexity: O(GM^3)
      Vchol = cholesky_decompose(
          diag_matrix(rep_vector(1/psi[g], M))
          + diag_post_multiply(Aut, Binv) * Aut');
      //
      // Calculate term in quadratic form part of likelihood
      // Complexity: O(GM^2 + GMC)
      b = mdivide_left_tri_low(Vchol, Aut * Binvy);
      //
      // Calculate determinant of the covariance
      log_det_cov =
        2*sum(log(diagonal(Vchol)))  // Twice log determinant of Cholesky
        + M * log(psi[g])            // Add log determinant of psi_g I
        - sum(log(Binv));            // Add log determinant of B
      //
      // Increment log probability with multivariate normal log likelihood
      // Complexity: O(GC)
      target += -.5 * (log_det_cov + dot_product(expradj[g], Binvy) - dot_self(b));
    }
}
