functions {
  real add_two(real x) {
    return x + 2;
  }
  
  // ------------------------------------------------------------
    // [NEW] SPD logdet with symmetrization + jitter (stability)
  // ------------------------------------------------------------
    real logdet_spd_jitter(matrix M, real jitter) {
      int d = rows(M);
      matrix[d,d] S = 0.5 * (M + M') + jitter * diag_matrix(rep_vector(1.0, d));
    return log_determinant_spd(S);
  }

  // ------------------------------------------------------------
  // [OLD] eigen-based inverse sqrt (numerically fragile in ADVI)
  // ------------------------------------------------------------
  /*
  matrix inv_sqrt_mat(matrix M) {
    int d = rows(M);
    vector[d] lambda = eigenvalues_sym(M);
    matrix[d,d] P = eigenvectors_sym(M);
    vector[d] inv_sqrt_lambda = inv_sqrt(lambda);
    return P * diag_matrix(inv_sqrt_lambda) * P';
    }
  */
    
    // ------------------------------------------------------------
    // [NEW] chol-based orthonormalization: CA*(CA'CA)^(-1/2) = CA * (chol(CA'CA)')^{-1}
  // This avoids eigenvalues/eigenvectors and is much more stable.
  // ------------------------------------------------------------
  matrix orthonormalize_from_C(matrix C, real jitter) {
    int k = cols(C);
    matrix[k,k] CtC = crossprod(C);
    // L is lower-triangular Cholesky of symmetrized CtC + jitter*I
    matrix[k,k] L = cholesky_decompose(0.5 * (CtC + CtC') + jitter * diag_matrix(rep_vector(1.0, k)));
// compute (L')^{-1} via triangular solves:
    matrix[k,k] Linv = mdivide_left_tri_low(L, diag_matrix(rep_vector(1.0, k))); // L^{-1}
    matrix[k,k] LinvT = Linv';                                                 // (L')^{-1}
    return C * LinvT; // C*(L')^{-1}
}

matrix get_gamma(matrix A) {
  int u = cols(A);
  int r_minus_u = rows(A);
  int r = r_minus_u + u;
  
  matrix[u,u] I_u = diag_matrix(rep_vector(1.0, u));
  matrix[r,u] CA = append_row(I_u, A);
  
  // --- [NEW] chol-based version ---
    return orthonormalize_from_C(CA, 1e-8);
  // --- [OLD] eigen-based version ---
    // matrix[u,u] CAtCA_minus_half;
  // CAtCA_minus_half = inv_sqrt_mat(crossprod(CA));
  // return CA * CAtCA_minus_half;
}

matrix get_gamma0(matrix A) {
  int u = cols(A);
  int r_minus_u = rows(A);
  int r = r_minus_u + u;
  
  matrix[r_minus_u,r_minus_u] I_r_minus_u = diag_matrix(rep_vector(1.0, r_minus_u));
  matrix[r,r_minus_u] DA = append_row(-A', I_r_minus_u);

    // --- [NEW] chol-based version ---
    return orthonormalize_from_C(DA, 1e-8);
    // --- [OLD] eigen-based version ---
    // matrix[r_minus_u,r_minus_u] DAtDA_minus_half;
    // DAtDA_minus_half = inv_sqrt_mat(crossprod(DA));
    // return DA * DAtDA_minus_half;
  }
}

data {
  int<lower=0> n;
  int<lower=0> r;
  int<lower=0> u;
  int<lower=0> r_minus_u;
  int<lower=0> nu;
  int<lower=0> nu0;
  matrix[u,u] Psi;
  matrix[r_minus_u,r_minus_u] Psi0;
  matrix[r,r] G_tilde;
  matrix[r,r] YctYc;
  matrix[r_minus_u,u] A0;
  matrix[r_minus_u,r_minus_u] K_inv;
  matrix[u,u] L_inv;
}

parameters {
  matrix[r_minus_u,u] A;
}

transformed parameters {
  matrix[r, u] gamma_A  = get_gamma(A);
  matrix[r, r_minus_u] gamma0_A = get_gamma0(A);
}

model {
  // --- [NEW] symmetrize inputs once (helps when supplied matrices are nearly symmetric) ---
  matrix[r,r] Gsym = 0.5 * (G_tilde + G_tilde');
  matrix[r,r] Ysym = 0.5 * (YctYc   + YctYc');

  // --- term1: logdet SPD with jitter (stable) ---
  // real term1 = (-0.5) * (n+nu-1) *
  //   log_determinant((gamma_A' * G_tilde * gamma_A) + Psi);              // [OLD]
real term1 = (-0.5) * (n + nu - 1) *
  logdet_spd_jitter((gamma_A' * Gsym * gamma_A) + Psi, 1e-8);            // [NEW]

  // --- term2: logdet SPD with jitter (stable) ---
  // real term2 = (-0.5) * (n+nu0-1) *
  //   log_determinant((gamma0_A' * YctYc * gamma0_A) + Psi0);              // [OLD]
real term2 = (-0.5) * (n + nu0 - 1) *
  logdet_spd_jitter((gamma0_A' * Ysym * gamma0_A) + Psi0, 1e-8);          // [NEW]

  // --- term3: prior term (kept as-is) ---
  real term3 = (-0.5) * sum(K_inv .* ((A - A0) * L_inv * (A - A0)'));

target += (term1 + term2 + term3);
}
