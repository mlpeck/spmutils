// calzetti attenuation
// precomputed emission lines

functions {
  vector calzetti_mod(vector lambda, real tauv, real delta) {
    int nl = rows(lambda);
    real lt;
    vector[nl] fk;
    for (i in 1:nl) {
      lt = 5500./lambda[i];
      fk[i] = -0.10177 + 0.549882*lt + 1.393039*pow(lt, 2) - 1.098615*pow(lt, 3)
            +0.260618*pow(lt, 4);
      fk[i] *= pow(lt, delta);
    }

    return exp(-tauv * fk);
  }
  
  real sum_ll(int[] ind, int start, int end, matrix sp_st, matrix sp_em, 
              vector gflux, vector g_std, vector lambda, real a, real tauv, real delta, 
              vector b_st_s, vector b_em) {
    return normal_lpdf(gflux[start:end] | a * (sp_st[start:end, :]*b_st_s) .* calzetti_mod(lambda[start:end], tauv, delta) 
                    + sp_em[start:end, :]*b_em, g_std[start:end]);
  }

}
data {
    int<lower=1> nt;  //number of stellar ages
    int<lower=1> nz;  //number of metallicity bins
    int<lower=1> nl;  //number of data points
    int<lower=1> n_em; //number of emission lines
    vector[nl] lambda;
    vector[nl] gflux;
    vector[nl] g_std;
    real norm_g;
    vector[nt*nz] norm_st;
    real norm_em;
    matrix[nl, nt*nz] sp_st; //the stellar library
    vector[nt*nz] dT;        // width of age bins
    matrix[nl, n_em] sp_em; //emission line profiles
}
transformed data {
  int grainsize = 1;
  int ind[nl] = rep_array(1, nl);
}
parameters {
    real a;
    simplex[nt*nz] b_st_s;
    vector<lower=0>[n_em] b_em;
    real<lower=0> tauv;
    real<lower= -0.5> delta;
}
model {
    b_em ~ normal(0, 100.);
    a ~ normal(1, 10.);
    tauv ~ normal(0, 1.);
    delta ~ normal(0., 0.1);
    target += reduce_sum(sum_ll, ind, grainsize, sp_st, sp_em, gflux, g_std, lambda,
                          a, tauv, delta, b_st_s, b_em);
    }

// put in for posterior predictive checking and log_lik

generated quantities {
    vector[nt*nz] b_st;
    vector[nl] gflux_rep;
    vector[nl] mu_st;
    vector[nl] mu_g;
    vector[nl] log_lik;
    real ll;
    
    b_st = a * b_st_s;
    
    mu_st = (sp_st*b_st) .* calzetti_mod(lambda, tauv, delta);
    mu_g = mu_st + sp_em*b_em;
    
    for (i in 1:nl) {
        gflux_rep[i] = norm_g*normal_rng(mu_g[i] , g_std[i]);
        log_lik[i] = normal_lpdf(gflux[i] | mu_g[i], g_std[i]);
        mu_g[i] = norm_g * mu_g[i];
    }
    ll = sum(log_lik);
}

