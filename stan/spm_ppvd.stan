functions {
  vector calzetti(vector lambda, real tauv) {
    int nl = rows(lambda);
    vector[nl] lt = 10000. * inv(lambda);
    vector[nl] fk;
    fk = -0.2688+0.7958*lt-4.785e-2*(lt .* lt)-6.033e-3*(lt .* lt .* lt)
            +7.163e-4*(lt .* lt .* lt .* lt);

    return exp(-tauv * fk);
  }

  vector convolve(vector x, vector kernel) {
    int nr = rows(x);
    int kl = rows(kernel);
    int kl1 = kl-1;
    vector[nr-kl1] y;
    for (i in 1:(nr-kl1)) {
        y[i] = dot_product(x[i:(i+kl1)], kernel);
    }
    return y;
  }
  matrix emprof(vector lpl, vector lp_em, real sigma_em, real voff_em, real[] h) {
    int nr = rows(lpl);
    int n_em = rows(lp_em);
    matrix[nr, n_em] lprof;
    vector[nr] y;
    vector[nr] H3;
    vector[nr] H4;
    
    for (j in 1:n_em) {
      y = (lpl - lp_em[j] - voff_em)/sigma_em;
      H3 = (2.0*(y .* y .* y) - 3.0*y)/sqrt(3.);       //3rd Hermite polynomial
      H4 = (4.0*(y .* y .* y .* y) - 12.0*(y .* y) + 3.)/2./sqrt(6.); //4th Hermite polynomial
      lprof[ , j] = exp(-(y .* y)/2).*(1.+h[1]*H3+h[2]*H4)/sigma_em/sqrt(2.*pi());
    }
    return lprof;
  }
}
data {
  int<lower=1> nf;  //number of flux values
  int<lower=1> nr; //rows in pc matrix
  int<lower=1> nt;  //number of stellar ages
  int<lower=1> nz;  //number of metallicity bins
  int<lower=1> n_em; //number of emission lines
  int<lower=1> kl; //kernel length
  
  vector[nf] lambda;
  vector[nf] gflux;
  vector<lower=0>[nf] g_std;
  int<lower=1> ins[nf];  //indexes of the non-missing flux values
  matrix[nr, nt*nz] sp_st;     //predictors
  vector[nt*nz] dT;        // width of age bins
  vector[n_em] lambda_em;
  real norm_em;  //normalizing constant to make emission line contributions <~ 1
  
}
transformed data {
  vector[nf] lpl = (10000./log(10.)) * log(lambda);
  vector[n_em] lp_em = (10000./log(10.)) * log(lambda_em);
}
parameters {
  simplex[kl] kernel;
  vector<lower=0>[nt*nz] b_st;
  vector[n_em] b_em;
  real<lower=1.> sigma_em;
  real voff_em;
  real h[2];
  real<lower=0> tauv;
}
model {
  vector[nr-kl+1] gmod;
  matrix[nf, n_em] lprof;

  b_st ~ cauchy(0., 4. * dT);
  b_em ~ cauchy(0., 1.);
  sigma_em ~ normal(0., 10.);
  voff_em ~ normal(0., 10.);
  h ~ normal(0., 0.25);
  tauv ~ normal(0., 1.);
  
  gmod = convolve(sp_st * b_st, kernel);
  lprof = norm_em * emprof(lpl, lp_em, sigma_em, voff_em, h); 
  gflux ~ normal(gmod[ins] .* calzetti(lambda, tauv)
              + lprof * b_em, g_std);
}
generated quantities {
  vector[nf] mu_st;
  vector[nf] mu_g;
  vector[nf] gflux_rep;
  vector[nf] log_lik;
  
  mu_st = (convolve(sp_st * b_st, kernel)[ins]) .* calzetti(lambda, tauv);
  mu_g = mu_st + norm_em * emprof(lpl, lp_em, sigma_em, voff_em, h) * b_em;
  for (i in 1:nf) {
    gflux_rep[i] = normal_rng(mu_g[i], g_std[i]);
    log_lik[i] = normal_lpdf(gflux[i] | mu_g[i], g_std[i]);
  }
}

