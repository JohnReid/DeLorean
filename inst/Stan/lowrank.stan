functions {
    //
    // Periodic function
    //
    real
    periodise(real r, real period) {
        return period / 2 * sin(r * pi() / period);
    }
    //
    // Squared exponential covariance function
    //
    real
    se_cov(real r, real l) {
        return exp(- pow(r / l, 2) / 2);
    }
    //
    // Matern nu=3/2 covariance function
    //
    real
    matern32_cov(real r, real l) {
        real x;
        x = sqrt(3.) * fabs(r / l);
        return (1 + x) * exp(- x);
    }
    //
    // Matern nu=5/2 covariance function
    //
    real
    matern52_cov(real r, real l) {
        real x;
        x = sqrt(5.) * fabs(r / l);
        return (1 + x + pow(x, 2) / 3) * exp(- x);
    }
    //
    // Covariance function
    //
    real
    cov_fn(real r, int periodic, real period, real l) {
        real rprime;
        if (periodic) {
            rprime = periodise(r, period);
        } else {
            rprime = r;
        }
        return matern32_cov(rprime, l);
        // return se_cov(rprime, l);
        // return se_cov(rprime, l);
        // return matern52_cov(rprime, l);
        // return matern52_cov(rprime, l);
    }
    //
    // Calculate symmetric covariance matrix
    //
    matrix
    cov_symmetric(row_vector tau, int periodic, real period, real l) {
        matrix[cols(tau),cols(tau)] result;
        for (c1 in 1:cols(tau)) {
            for (c2 in 1:c1) {
                result[c1,c2] = cov_fn(tau[c2] - tau[c1], periodic,
                                        period, l);
                if(c1 != c2) {
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
    cov(row_vector tau1,
        row_vector tau2,
        int periodic,
        real period,
        real l
    ) {
        matrix[cols(tau1),cols(tau2)] result;
        for (c1 in 1:cols(tau1)) {
            for (c2 in 1:cols(tau2)) {
                result[c1,c2] = cov_fn(tau2[c2] - tau1[c1], periodic, period, l);
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
  int<lower=2> C;  // Number of cells
  int<lower=2> G;  // Number of genes
  int<lower=0> H;  // Number of held out genes
  int<lower=0> M;  // Number of inducing pseudotimes
  //
  // Data
  //
  // Time data
  int periodic; // Are the expression patterns periodic?
  real period;  // Cyclic period
  row_vector<lower=0>[C] time;  // Time index for cell c
  // Expression data
  vector[C] expr[G+H];
  row_vector[G+H] phi;  // gene mean
  //
  // Inducing pseudotimes
  //
  row_vector[M] u; // pseudotimes
  //
  // Hyperparameters
  //
  real mu_psi;  // Mean of log between time variation, psi
  real<lower=0> sigma_psi;  // S.d. of log between time variation, psi
  real mu_omega;  // Mean of log within time variation, omega
  real<lower=0> sigma_omega;  // S.d. of log within time variation, omega
  real l;  // Length scale
  real<lower=0> sigma_tau;  // Standard deviation for pseudotime
  //
  // Held out parameters
  //
  row_vector<lower=0>[H] heldout_psi;    // Between time variance
  row_vector<lower=0>[H] heldout_omega;  // Within time variance
  //
  // Generated quantities
  //
  int numtest; // Number of test inputs for predictions
  row_vector[numtest] testinput; // Test inputs for predictions
}
transformed data {
  //
  // Subscripts:
  //   s - test inputs (star)
  //   u - inducing inputs
  //   t - pseudotime (observed/training) inputs
  vector[C] expradj[G+H];  // Mean adjusted expression
  matrix[M,M] Kuu; // Covariance between inducing points
  matrix[M,numtest] Kus; // Covariance between inducing and test inputs
  cholesky_factor_cov[M] KuuChol; // pseudotime covariance Cholesky
  Kuu = cov_symmetric(u, periodic, period, l);
  Kus = cov(u, testinput, periodic, period, l);
  KuuChol = cholesky_decompose(Kuu);
  // Calculate adjusted gene expression
  for (g in 1:(G+H)) {
    expradj[g] = expr[g] - phi[g];
  }
}
parameters {
  row_vector[C] tauoffsets;    // Pseudotime
  row_vector<lower=0>[G] psi;    // Between time variance
  row_vector<lower=0>[G] omega;  // Within time variance
}
transformed parameters {
  row_vector[C] tau;    // Pseudotime
  tau = time + tauoffsets;
}
model {
  //
  // Approximation to pseudotime covariance matrix
  matrix[M,C] Aut;  // Sqrt of approximation
  vector[C] KminusQdiag;   // Diagonal adjustment to approximate covariance
  //
  // Calculate the square root of Qtautau. Complexity: O(CM^2)
  Aut = mdivide_left_tri_low(KuuChol, cov(u, tau, periodic, period, l));
  KminusQdiag = 1 - columns_dot_self(Aut)';
  // Check that diag(Ktautau - Qtt) is positive
  for (c in 1:C) {
      if (KminusQdiag[c] < 0) {
          reject("KminusQdiag must be greater or equal to 0. : ",
                  KminusQdiag[c]);
      }
  }
  //
  // Sample gene-specific factors
  target += lognormal_lpdf(psi|mu_psi, sigma_psi);
  target += lognormal_lpdf(omega|mu_omega, sigma_omega);
  //
  // Sample pseudotime
  target += normal_lpdf(tauoffsets|0, sigma_tau);  // Pseudotime
  //
  // Expression values for each gene
  for (g in 1:G) {
    vector[C] Binv;
    vector[C] Binvy;
    vector[M] b;  // Aut * B-1 * y
    matrix[M,M] Vchol;  // Low dimension decomposition
    real log_det_cov;  // The logarithm of the determinant of the covariance
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
      + M * log(psi[g]) // Add log determinant of psi_g I
      - sum(log(Binv)); // Add log determinant of B
    //
    // Increment log probability with multivariate normal log likelihood
    // Complexity: O(GC)
    target += -.5 * (log_det_cov + dot_product(expradj[g], Binvy) - dot_self(b));
  }
}
generated quantities {
  vector[numtest] predictedmean[G+H];
  vector[numtest] predictedvar[G+H];
  real logmarglike[G+H];
  {
    matrix[M,C] Kut;
    Kut = cov(u, tau, periodic, period, l);
    //
    // For each gene (including held out genes)
    for (g in 1:(G+H)) {
      matrix[numtest,M] KsuSigmag;
      real psi_g;
      real omega_g;
      real ratio;
      //
      // Get the temporal variance and noise parameters
      if (g <= G) { // Sampled parameters
          psi_g = psi[g];
          omega_g = omega[g];
      } else { // Parameters for held out genes
          psi_g = heldout_psi[g-G];
          omega_g = heldout_omega[g-G];
      }
      //
      // Pre-calculate some useful terms
      ratio = psi_g/omega_g;
      KsuSigmag = Kus' / (ratio*tcrossprod(Kut) + Kuu);
      //
      // Calculate the predicted mean
      predictedmean[g] = ratio*KsuSigmag*Kut*expradj[g];
      //
      // Calculate the predicted variance
      predictedvar[g] = psi_g * (
        1
        - calc_Q_diag(KuuChol, cov(u, testinput, periodic, period, l))
        + diagonal(KsuSigmag * Kus));
      //
      // Calculate log marginal likelihood
      logmarglike[g] = 0;
    }
  }
}
