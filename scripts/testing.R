source('scripts/funcs.R')

rho_fitness_single = function(lambda, g, pars, b, k0){
  n1 = 1:k0
  n2 = 0:(k0-1)
  out1 = sum(exp(cumsum(log(1-lambda/(lambda+g(n1-1,pars))))))
  out2 = exp(sum(log(1-lambda/(lambda+g(n2,pars)))))
  return(out1 + g(k0,pars)*out2/(lambda-b))
}
rho_fitness_dist = Vectorize(function(lambda, g, fitsample, b, k0, n=10){
  smps = fitsample(n)
  rhos = apply(smps, 1, rho_fitness_single, lambda=lambda, g=g, b=b, k0=k0)
  return(mean(rhos)-1)
}, vectorize.args = 'lambda')
find_lambda_fitness = function(g, fitsample, b, k0, n=1e4){
  rt = uniroot(rho_fitness_dist, c(0.1,10), g=g, fitsample=fitsample, b=b, k0=k0, n=n, extendInt = 'yes')$root
  return(rt)
}

S_fitness_single = Vectorize(function(x, g,pars, lambda, k0, b){
  if(x<0){
    return(1)
  }
  if(x>=k0){
    k0s = 0:(k0-1)
    
    phi_g = lgamma((lambda+g(k0,pars))/b)  + sum(log(g(k0s,pars)/(lambda+g(k0s,pars))))-
      (lgamma(lambda/b)+lgamma(g(k0,pars)/b))
    
    return(exp(phi_g+lbeta(g(k0,pars)/b + x-k0+1, lambda/b)))
  }
  if(x>0){
    ks = 0:(x)
    return(prod(g(ks,pars)/(lambda+g(ks,pars))))
  }
  if(x==0){
    return(g(x,pars)/(lambda+g(x,pars)))
  }
  
},vectorize.args = 'x')

S_fitness = Vectorize(function(x, g, fitsample, lambda, k0, b, n=1000){
  smps = fitsample(n)
  return(mean(apply(smps, 1, S_fitness_single, x=x, g=g, lambda=lambda, k0=k0, b=b)))
}, vectorize.args = 'x')

sample_gpa_tree_fitness = function(n, g,fitsample, b, k0,quiet=T){
  root_fitness = fitsample(1)
  df = data.frame(deg = c(0),
                  pref = c(g(0, root_fitness)),
                  fitness = root_fitness)
  
  for(i in 1:(n-1)){
    if(!quiet) message(i)
    selected = sample(1:nrow(df),1, prob = df$pref)
    df$deg[selected] = df$deg[selected] + 1
    df$pref[selected]  = ifelse(df$deg[selected]<=k0, g(df$deg[selected], df[selected,-(1:2)]), 
                                g(df$deg[selected], df[selected,-(1:2)])+ b*(df$deg[selected]-k0))
    newfit = fitsample(1)
    df = rbind(df, c(0, g(0,newfit),newfit))
  }
  return(df)
}


# -------------------------------------------------------------------------
fitsample = function(n){
  pars1 = runif(n, 0,1.5)
  pars2 = runif(n)
  return(cbind(pars1, pars2))
}

g = function(x,pars){
  return(x^pars[1] + pars[2])
}


b = 1
k0 = 25
lambda = find_lambda_fitness(g, fitsample, b, k0)

x = 1:100
y = S_fitness(x, g, fitsample, lambda, k0, b)
plot(x, y, log='xy', ylim=c(1e-4, 1), pch=20)

df = sample_gpa_tree_fitness(10000, g, fitsample, b, k0, quiet=F)
points(twbfn::deg_surv(df$deg), pch=20, col='red')













