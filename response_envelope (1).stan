functions {
  real add_two(real x) {
    return x + 2;
  }

  real logdet_spd_jitter(matrix M, real jitter){
    int d = rows(M);
    matrix[d,d] S = 0.5*(M + M') + jitter * diag_matrix(rep_vector(1.0, d));
    return log_determinant_spd(S);
  }
  
  matrix inv_sqrt_mat(matrix M) {
    int d = rows(M);
    vector[d] lambda = eigenvalues_sym(M + 1e-6 * diag_matrix(rep_vector(1.0, d)));
    matrix[d,d] P = eigenvectors_sym(M + 1e-6 * diag_matrix(rep_vector(1.0, d)));
    vector[d] inv_sqrt_lambda = inv_sqrt(lambda);
    return P * diag_matrix(inv_sqrt_lambda) * P';
  }

  matrix get_gamma(matrix A) {
    int u = cols(A);
    int r_minus_u = rows(A);
    int r = r_minus_u + u;
    vector[u] one_u;
    matrix[u,u] I_u;
    matrix[r,u] CA;
    matrix[u,u] CAtCA_minus_half;

    u = cols(A);
    r_minus_u = rows(A);

    for (i in 1:u) {
      one_u[i] = 1.0;
    }

    I_u = diag_matrix(one_u);
    CA = append_row(I_u, A);
    CAtCA_minus_half = inv_sqrt_mat(crossprod(CA));

    return CA * CAtCA_minus_half;
  }


  matrix get_gamma0(matrix A) {
    int u = cols(A);
    int r_minus_u = rows(A);
    int r = r_minus_u + u;
    vector[r_minus_u] one_r_minus_u;
    matrix[r_minus_u,r_minus_u] I_r_minus_u;
    matrix[r,r_minus_u] DA;
    matrix[r_minus_u,r_minus_u] DAtDA_minus_half;

    u = cols(A);
    r_minus_u = rows(A);

    for (i in 1:r_minus_u) {
      one_r_minus_u[i] = 1.0;
    }

    I_r_minus_u = diag_matrix(one_r_minus_u);
    DA = append_row(-A', I_r_minus_u);
    DAtDA_minus_half = inv_sqrt_mat(crossprod(DA));

    return DA * DAtDA_minus_half;
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


model{
  matrix[u,u] Iu = diag_matrix(rep_vector(1.0, u));
  matrix[r_minus_u,r_minus_u] Iru = diag_matrix(rep_vector(1.0, r_minus_u));

  matrix[r,u] CA = append_row(Iu, A);
  matrix[r,r_minus_u] DA = append_row(-A', Iru);

  matrix[r,r] Gsym = 0.5*(G_tilde + G_tilde');
  matrix[r,r] Ysym = 0.5*(YctYc + YctYc');

  // Psi = psi * Iu, Psi0 = psi0 * Iru 라고 가정
  real psi  = Psi[1,1];
  real psi0 = Psi0[1,1];

  matrix[u,u] CAtCA = crossprod(CA);
  matrix[r_minus_u,r_minus_u] DAtDA = crossprod(DA);

  matrix[u,u] S1  = crossprod(CA, Gsym)  * CA + psi  * CAtCA;
  matrix[r_minus_u,r_minus_u] S0 = crossprod(DA, Ysym) * DA + psi0 * DAtDA;

  real ld1 = logdet_spd_jitter(S1,  1e-8) - logdet_spd_jitter(CAtCA, 1e-8);
  real ld0 = logdet_spd_jitter(S0,  1e-8) - logdet_spd_jitter(DAtDA, 1e-8);

  real term1 = -0.5 * (n + nu  - 1) * ld1;
  real term2 = -0.5 * (n + nu0 - 1) * ld0;

  // MN prior (정확한 형태)
  real term3 = -0.5 * sum( (K_inv * (A - A0) * L_inv) .* (A - A0) );

  target += term1 + term2 + term3;
}
