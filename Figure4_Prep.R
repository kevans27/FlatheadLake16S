sink("satAndSine.stan")
cat("
      data {
      int <lower=1> N; //number of data points
      int <lower=1> J; //number of depths
      int <lower=1> K; //number of samples per depth
      vector [K] Days; //input days between samples
      matrix[K,J] Jaccs; //jaccard dissimilarity model matrix
      }
      
      parameters {
      vector <lower=0, upper=1> [J] alph;
      real <lower=0> mu_alph;
      real <lower=0> sd_alph;
      vector <lower=0.000001> [J] alphPB;
      real <lower=0> mu_alphPB;
      real <lower=0> sd_alphPB;
      vector <lower=0.000001> [J] PBS;
      real <lower=0> mu_PBS;
      real <lower=0> sd_PBS;
      vector <lower=-3.1, upper=2.8> [J] omic;
      real <lower=-3.1, upper=2.8> mu_omic;
      real <lower=0> sd_omic;
      vector <lower=0> [J] sigm;
      real <lower=0> mu_sigm;
      real <lower=0> sd_sigm;
      vector [J] interc;
      real mu_interc;
      real <lower=0> sd_interc;
      }
      
      model { 
      //priors. 
      alph ~ normal(mu_alph, sd_alph);
      mu_alph ~ normal(0.25, 0.1);
      sd_alph ~ normal(0, 0.2);
      alphPB ~ normal(mu_alphPB, sd_alphPB);
      mu_alphPB ~ normal(0.02, 0.01);
      sd_alphPB ~ normal(0.2, 0.1);
      PBS ~ normal(mu_PBS, sd_PBS);
      mu_PBS ~ normal(0.6, 0.2);
      sd_PBS ~ normal(0.2, 0.1);
      omic ~ normal(mu_omic, sd_omic);
      mu_omic ~ normal(-0.5, 1);
      sd_omic ~ normal(0.5, 0.5);
      sigm ~ normal(mu_sigm, sd_sigm);
      mu_sigm ~ normal(0.6, 0.1);
      sd_sigm ~ normal(0.1, 0.05);
      interc ~ normal(0.2, 0.2);
      mu_interc ~ normal(0.2, 0.2);
      sd_interc ~ normal(0.2, 0.2);

      //likelihood    	
      for (j in 1:J) {
        for (k in 1:K) {
            Jaccs[k,j] ~ normal(PBS[j]*(1-exp(-alphPB[j]*Days[k]/PBS[j]))+alph[j]*sin(6.283*Days[k]/365+omic[j])+interc[j], sigm[j]);
        }
      }
      
      }",
    fill=TRUE)
sink()



sink("sineOnly.stan")
cat("
      data {
      int <lower=1> N; //number of data points
      int <lower=1> J; //number of depths
      int <lower=1> K; //number of samples per depth
      vector [K] Days; //input days between samples
      matrix[K,J] Jaccs; //jaccard dissimilarity model matrix
      }
      
      parameters {
      vector <lower=0> [J] alph;
      real <lower=0> mu_alph;
      real <lower=0> sd_alph;
      vector <lower=-3.1, upper=2.8> [J] omic;
      real mu_omic;
      real <lower=0> sd_omic;
      vector <lower=0> [J] sigm;
      real <lower=0> mu_sigm;
      real <lower=0> sd_sigm;
      vector [J] interc;
      real mu_interc;
      real <lower=0> sd_interc;
      }
      
      model { 
      //priors. 
      alph ~ normal(mu_alph, sd_alph);
      mu_alph ~ normal(0.1, 0.1);
      sd_alph ~ normal(0, 0.2);
      omic ~ normal(mu_omic, sd_omic);
      mu_omic ~ normal(-2, 1.5);
      sd_omic ~ normal(0.5, 1);
      sigm ~ normal(mu_sigm, sd_sigm);
      mu_sigm ~ normal(0.06, 0.1);
      sd_sigm ~ normal(0.05, 0.05);
      interc ~ normal(mu_interc, sd_interc);
      mu_interc ~ normal(0.75, 0.1);
      sd_interc ~ normal(0.05, 0.05);

      //likelihood    	
      for (j in 1:J) {
        for (k in 1:K) {
            Jaccs[k,j] ~ normal((alph[j]*sin((6.283*Days[k])/365+omic[j]))+interc[j], sigm[j]);
        }
      }
      
      }",
    fill=TRUE)
sink()