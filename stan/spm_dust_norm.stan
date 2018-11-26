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
    vector[nl] lambda;
    vector[nl] gflux;
    vector[nl] g_std;
    matrix[nl, nt*nz] sp_st; //the stellar library
    vector[nt*nz] dT;        // width of age bins
    matrix[nl, n_em] sp_em; //emission line profiles
}
parameters {
    vector<lower=0>[nt*nz] b_st;
    vector<lower=0>[n_em] b_em;
    real<lower=0> tauv;
}
model {
    b_em ~ cauchy(0, 1);
    b_st ~ cauchy(0, 4.*dT);
    tauv ~ normal(0, 1);
    gflux ~ normal((sp_st*b_st) .* calzetti(lambda, tauv) 
                    + sp_em*b_em, g_std);
}

// put in for posterior predictive checking and log_lik

generated quantities {
    vector[nl] gflux_rep;
    vector[nl] mu_st;
    vector[nl] mu_g;
    vector[nl] log_lik;
    
    mu_st = (sp_st*b_st) .* calzetti(lambda, tauv);
    mu_g = mu_st + sp_em*b_em;
    
    for (i in 1:nl) {
        gflux_rep[i] = normal_rng(mu_g[i] , g_std[i]);
        log_lik[i] = normal_lpdf(gflux[i] | mu_g[i], g_std[i]);
        
    }
}

