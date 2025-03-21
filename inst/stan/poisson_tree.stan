functions{

    real x_to_t2(real[] x, int[] parentidx, int[] xidx,real[] tip_min_age, real[] t, int i);
    //real[] xx_to_t2(real[] x, int[] parentidx, int[] xidx,real[] tip_min_age, real[] t, int i);
    //TODO (See above). The following is a transcription of the C code that supported pass-by-reference
    //pass-by-ref is not supported in stan - so we can set each t multiple times. Should be possible to fix
    //by passing back a vector.
    real x_to_t2(real[] x, int[] parentidx, int[] xidx,real[] tip_min_age, real[] t, int i){
      real factor;
      if(t[i]>=0){
        return t[i];
      }
      if(parentidx[i]<1){
        if(xidx[i]<1){
          //tip is directly ancestral to root..
          return tip_min_age[i];
        }else{
          //t[i]=x[xidx[i]]*tip_min_age[i];
          return x[xidx[i]]*tip_min_age[i];
        }
      }else{
        real tmp=0;
        int k=i;
        while(parentidx[k]>=1){
          k=parentidx[k];
          tmp+=x_to_t2(x,parentidx,xidx,tip_min_age,t,k);
        }
        if(xidx[i]>=1){
          //interior
          factor=tip_min_age[i]-tmp;
          //t[i]=x[xidx[i]]*factor;
          return x[xidx[i]]*factor;
        }else{
          //t[i]=tip_min_age[i]-tmp;
          return tip_min_age[i]-tmp;
        }
        //return t[i];
      }
    }
    //Need to make sure all indices are 1+
      vector x_to_t(real[] x, int[] parentidx, int[] xidx, real[] tip_min_age, int n){
        int N=n;
        real t[N];
        for(i in 1:N){
          t[i]=-1;
        }
        for(i  in 1:N){
          t[i]=x_to_t2(x,parentidx,xidx,tip_min_age,t,i);
        }
        return to_vector(t);
      }

      vector logisticMean(real L,real k,real midpoint,vector a,vector b){
        return ((L/k)*(log(1+exp(k*(b-midpoint)))-log(1+exp(k*(a-midpoint))))) ./ (b-a);
      }

      vector getExtraLambdaRates(vector t,int[] parentidx,int N){
        real extralambda=149.1968140;
        real kg=-50.0161248;
        real mp=0.2246393;
        real ct=0;//current time
        int k;
        vector[N] t0;

        for(i in 1:N){
          t0[i]=-1;
        }
        for(i in 1:N){
          k=parentidx[i];
          ct=0.0;
          while(k>0){
            if(t0[k]>0){
              ct=ct+t0[k]+t[k];
              k=-1;
            }else{
              ct=ct+t[k];
              k=parentidx[k];
            }

          }
          t0[i]=ct;
        }
        for(i in 1:N){
          if (t0[i]<0){
          reject("x must not be negative; found x=", t0[i]);
          }

        }
        return logisticMean(extralambda,kg,mp,t0,t0+t);
      }



  }

data{
  int N; //num branches
  int NINT; //num internal branches
  int NLAMBDA;
  int parentidx[N];
  int xidx[N];
  real tip_min_age[N];
  int rates[N];
  int ratesp[N];
  vector[N] s; //sensitivity
  int m[N];
  int nh[NINT];
  vector[NINT] q;
  vector[NINT] concentration;
  int  idxcrossover[NLAMBDA-1];
  real lambda_est;
  real early_growth_model_on;
  vector[NLAMBDA-1] lb_xover;
  vector[NLAMBDA-1] ub_xover;
  
}

parameters {
  real<lower=0.0001,upper=0.9999> x[NINT];
  vector<lower=0.001,upper=0.999>[N] S;
  vector<lower=1,upper=200>[NLAMBDA] lambda;
  real<lower=0.05,upper=0.999> p;  //nb 1/overdispersion
  vector<lower=0.00001,upper=0.99999>[NLAMBDA-1] x0_raw; //fractional crossover..
}

transformed parameters {
  vector[NLAMBDA-1] x0 = lb_xover + (ub_xover-lb_xover) .* x0_raw;
}


model {
  vector[N] t;
  vector[N] t0;
  vector[N] lambda_per_branch;
  vector[N] tmp;
  vector[N] x0v;
  x0v=rep_vector(0.0,N);
  x ~ beta(concentration .* q ./ (1-q),concentration);
  S ~ beta(100,100*(1-s) ./ s);
  
  for(i in 1:(NLAMBDA-1)){
    //print("x0[",i,"]=",x0[i],":lb=",lb_xover[i],":ub=",ub_xover[i]);
    x0v[idxcrossover[i]]=x0[i];
  }
  lambda ~ normal(lambda_est,0.25*lambda_est);
  
  t=x_to_t(x, parentidx, xidx, tip_min_age, N);
  lambda_per_branch=(1-x0v) .* lambda[rates]+x0v .* lambda[ratesp]+early_growth_model_on*getExtraLambdaRates(t,parentidx,N);
  m ~ poisson(t .* lambda_per_branch .* S);
}

generated quantities {
  vector[N] ta;
  //vector[N] elb;
  ta=x_to_t(x, parentidx, xidx, tip_min_age, N);
  ///elb=ta .* ((1-alpha) .* lambda[rates]+alpha .* lambda[ratesp] + getExtraLambdaRates(ta,parentidx,N)) .* S;

}
