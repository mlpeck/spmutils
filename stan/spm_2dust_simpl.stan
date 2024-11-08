// 2 component dust model 
// calzetti attenuation
// precomputed emission lines

functions {
  vector calzetti(vector lambda, real tauv) {
    int nl = rows(lambda);
    vector[nl] lt = 10000. * inv(lambda);
    vector[nl] fk;
    fk = -0.2688+0.7958*lt-4.785e-2*(lt .* lt)-6.033e-3*(lt .* lt .* lt)
            +7.163e-4*(lt .* lt .* lt .* lt);

    return exp(-tauv * fk);
  }

}
data {
    int<lower=1> nt;  //number of stellar ages
    int<lower=1> nz;  //number of metallicity bins
    int<lower=1> nl;  //number of data points
    int<lower=1> n_em; //number of emission lines
    int ind_young[nz];  //indexes of young ssp's
    vector[nl] lambda;
    vector[nl] gflux;
    vector[nl] g_std;
    real norm_g;
    vector[nt*nz] norm_st;
    real norm_em;
    matrix[nl, nt*nz] sp_st; //the stellar library
    matrix[nl, n_em] sp_em; //emission line profiles
}
parameters {
    real a;
    simplex[nt*nz] b_st_s;
    vector<lower=0>[nz] b_st_young;
    vector<lower=0>[n_em] b_em;
    real<lower=0> tauv;
}
model {
    b_em ~ normal(0, 100.);
    a ~ normal(1., 10.);
    b_st_young ~ normal(0, 1.);
    tauv ~ normal(0, 1);
    gflux ~ normal(a * (sp_st*b_st_s) .* calzetti(lambda, tauv) +
                   (sp_st[:, ind_young] * b_st_young) .* calzetti(lambda, tauv+1.) +                    
                   sp_em*b_em, g_std);
}

// put in for posterior predictive checking and log_lik

generated quantities {
    vector[nt*nz] b_st;
    vector[nl] gflux_rep;
    vector[nl] mu_g;
    vector[nl] log_lik;
    real ll;
    
    mu_g = a * (sp_st*b_st_s) .* calzetti(lambda, tauv) + 
           (sp_st[:, ind_young] * b_st_young) .* calzetti(lambda, tauv+1.) +
            sp_em*b_em;
    b_st = a*b_st_s;
    b_st[ind_young] += b_st_young;
    
    for (i in 1:nl) {
        gflux_rep[i] = norm_g * normal_rng(mu_g[i] , g_std[i]);
        log_lik[i] = normal_lpdf(gflux[i] | mu_g[i], g_std[i]);
        mu_g[i] = norm_g * mu_g[i];
    }
    ll = sum(log_lik);
}

