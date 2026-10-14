functions{

    // Transforms stick-breaking fractions, x, to branch durations, t, in a single O(N) pass.
    // Equivalent to x_to_t in poisson_tree.stan, but REQUIRES that every branch's parent appears
    // before it, i.e. parentidx[i] < i (checked in transformed data; ape "cladewise" order satisfies this).
    // start[i] is the time at which branch i begins (= sum of the durations of all its ancestors).
      vector x_to_t(array[] real x, array[] int parentidx, array[] int xidx, array[] real tip_min_age, int n){
        vector[n] t;
        vector[n] start;
        for(i in 1:n){
          if(parentidx[i]>=1){
            start[i]=start[parentidx[i]]+t[parentidx[i]];
          }else{
            start[i]=0;
          }
          if(xidx[i]>=1){
            //interior
            t[i]=x[xidx[i]]*(tip_min_age[i]-start[i]);
          }else{
            //tip
            t[i]=tip_min_age[i]-start[i];
          }
        }
        return t;
      }

      vector logisticMean(real L,real k,real midpoint,vector a,vector b){
        return ((L/k)*(log(1+exp(k*(b-midpoint)))-log(1+exp(k*(a-midpoint))))) ./ (b-a);
      }

      vector getExtraLambdaRates(vector t,array[] int parentidx,int N){
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
  array[N] int parentidx;
  array[N] int xidx;
  array[N] real tip_min_age;
  array[N] int rates;
  array[N] int ratesp;
  vector[N] s; //sensitivity
  array[N] int m;
  array[NINT] int nh;
  vector[NINT] q;
  vector[NINT] concentration;
  array[NLAMBDA-1] int idxcrossover;
  real lambda_est;
  real early_growth_model_on;
  vector[NLAMBDA-1] lb_xover;
  vector[NLAMBDA-1] ub_xover;
  
}

transformed data {
  // x_to_t relies on parents being listed before children.
  for(i in 1:N){
    if(parentidx[i]>=i){
      reject("Edges must be ordered with parents before children (parentidx[i] < i); failed at i=", i);
    }
  }
}

parameters {
  array[NINT] real<lower=0.0001,upper=0.9999> x;
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
