f = Vectorize(function(x, g, lambda, k0, b){
  if(x>=k0){
    k0s = 0:(k0-1)
    phi_g = lambda*gamma((lambda+g(k0))/b)*prod(g(k0s)/(lambda+g(k0s)))/
      (b*gamma(lambda/b + 1)*gamma(g(k0)/b))
    return(phi_g*beta(g(k0)/b + x-k0, lambda/b + 1))
  }
  if(x>0){
    ks = 0:(x-1)
    return(lambda/(lambda+g(x))*prod(g(ks)/(lambda+g(ks))))
  }
  if(x==0){
    return(lambda/(lambda+g(x)))
  }
}, vectorize.args = 'x')
S = Vectorize(function(x, g, lambda, k0, b){
  if(x<0){
    return(1)
  }
  if(x>=k0){
    k0s = 0:(k0-1)
    
    phi_g = lgamma((lambda+g(k0))/b)  + sum(log(g(k0s)/(lambda+g(k0s))))-
      (lgamma(lambda/b)+lgamma(g(k0)/b))
    
    return(exp(phi_g)*beta(g(k0)/b + x-k0+1, lambda/b))
  }
  if(x>0){
    ks = 0:(x)
    return(prod(g(ks)/(lambda+g(ks))))
  }
  if(x==0){
    return(g(x)/(lambda+g(x)))
  }
  
},vectorize.args = 'x')
sample_gpa_tree = function(n, g,b,k0,m=1){
  degs = c(0)
  prefs = c(g(0))
  for(i in 1:n){
    # message(i)
    selected = sample(1:length(degs), min(m, length(degs)), prob=prefs)
    degs = c(degs, 0)
    prefs = c(prefs, g(0))
    degs[selected] = degs[selected]+1
    prefs[selected] = ifelse(degs[selected]<k0, g(degs[selected]), g(k0) + b*(degs[selected]-k0))
  }
  return(degs)
}
rho = Vectorize(function(lambda, g, b, k0){
  n1 = 1:k0
  n2 = 0:(k0-1)
  out1 = sum(exp(cumsum(log(1-lambda/(lambda+g(n1-1))))))
  out2 = exp(sum(log(1-lambda/(lambda+g(n2)))))
  return(out1 + g(k0)*out2/(lambda-b))
}, vectorize.args = 'lambda')
rho_optim = function(lambda, g, b, k0){
  return(rho(lambda, g, b, k0)-1)
}
find_lambda2 = function(g, b, k0){
  out = try(uniroot(rho_optim, c(b,100), g=g,b=b,k0=k0, extendInt = 'yes')$root,silent=F)
  if(is.character(out)){
    return(NULL)
  }
  return(out)
}
samplepolylin = function(n,a,eps,b,k0,m=1){
  g01 = function(x){
    return((x)^a+eps)
  }
  degs = sample_gpa_tree(n,g01,b,k0,m)
  return(degs)
}
polylin = function(a,eps){
  return(function(x){
    return((x)^a+eps)
  })
}
llh_x = Vectorize(function(x, pars,lambda,minx = 0){
  a = pars[1]
  eps = pars[2]
  k0 = pars[3]
  b = pars[4]
  if(x==0){
    return(log(lambda) - log(lambda + eps))
  }else if(x<k0){
    return(log(lambda) - log(lambda + eps + (x)^a) +
             sum(log(1 - lambda/(lambda + eps + ((minx:(x-1)))^a))))
  }else{
    return(sum(log(1 - lambda/(lambda + eps + ((minx:(k0-1)) )^a))) +
             lbeta(x-k0 + ((k0)^a + eps)/b,1 + lambda/b) - lbeta(((k0)^a + eps)/b,lambda/b))
  }
},vectorize.args='x')
llh = function(dat, pars, lambda,minx=0){
  if(any(pars<0)){
    return(-Inf)
  }
  return(sum(dat[,2]*llh_x(dat[,1], pars, lambda,minx)))
}
prior = function(pars,xmax = 1e3){
  return(
    dgamma(pars[1],1,0.01,log=T)+
      dgamma(pars[2],1,0.01,log=T)+
      dunif(pars[3],1,1e4,log=T)+
      dgamma(pars[4],1,0.01,log=T)
  )
}
posterior = function(dat, pars, lambda,minx=min(dat[,1])){
  return(llh(dat, pars, lambda,minx) + prior(pars,xmax = max(dat[,1])))
}
posterior_k0 = Vectorize(function(dat, k0, pars, lambda,minx=min(dat[,1])){
  if(k0<=min(dat[,1])){
    return(-Inf)
  }
  if(k0>=max(dat[,1])){
    return(-Inf)
  }
  return(posterior(dat, c(pars[1],pars[2],k0,pars[4]), lambda,minx))
},vectorize.args = 'k0')
polylin_mcmc = function(N, dat,init =c(a = 1,eps=1,k0=dat[min(10,nrow(dat)-2),1],b=1),sig = matrix(c(5e-4, 5e-5, 3e-5, 5e-5, 2.5e-5, 4e-5, 3e-5, 4e-5, 1.5e-3),ncol=3,nrow=3),k0_jump = 2,stop_at = 1e4,  update_period = 2e2,quiet=T,burn.in=5e3,cov_scale = 0.1){
  par_mat = matrix(as.matrix(init),nrow=1)
  lambda_vec = c(2)
  post_vec = c(posterior(dat, par_mat[1,],lambda_vec[1]))
  lock = FALSE
  old_pars = par_mat[1,]
  accepted=numeric(N)
  st = 1
  for(i in 1:N){
    if(length(post_vec)%%update_period == 0 & !quiet){
      message(i,': ',sum(accepted[1:update_period + (i - update_period)])/update_period,'A ',' | ', paste0(signif(old_pars,2),sep='||'))
      lock=TRUE
    }
    #update covariance matrix
    if(length(post_vec)>update_period  & length(post_vec)>=burn.in){
      acc_rate  = sum(accepted[1:update_period + (i - update_period)])/update_period
      st = st*exp(cov_scale*(acc_rate-0.234))
      sig = st*(cov(par_mat[,c(1,2,4)]))
      # sig = (cov(par_mat[,c(1,2,4)]))*(2.38^2)/4
    }
    #propose new values
    cur_k0 = tail(par_mat[,3],1)
    # possk0  = max((cur_k0-k0_jump),2):(cur_k0+k0_jump)
    possk0 = dat[,1][max(1,(which(dat[,1]==cur_k0) - k0_jump)):min(nrow(dat),(which(dat[,1]==cur_k0) + k0_jump))]
    y = posterior_k0(dat, possk0, old_pars,tail(lambda_vec,1))
    # print(exp(y - max(y)))
    
    
    possk0_probs = exp(y - max(y))
    if(any(is.na(possk0_probs))){
      new_k0 = -1
    }else{
      new_k0 = sample(possk0,1,prob = possk0_probs)
    }
    new_pars0 = mvtnorm::rmvnorm(1, c(old_pars[-3]), sig)
    new_pars = c(new_pars0[1],new_pars0[2],new_k0,new_pars0[3])
    if(!any(new_pars<0)){
      #obtain lambda
      new_lambda = find_lambda2(polylin(new_pars[1], new_pars[2]), new_pars[4], new_pars[3])
      if(is.null(new_lambda)){
        par_mat = rbind(par_mat, old_pars)
        lambda_vec = c(lambda_vec, tail(lambda_vec,1))
        post_vec = c(post_vec, tail(post_vec, 1))
      }else{
        #obtain new_posterior
        new_post = posterior(dat,new_pars,new_lambda)
        logA = min(0,new_post - tail(post_vec,1))
        #accept/reject
        if(log(runif(1))<logA | runif(1)<0.05){
          old_pars=new_pars
          par_mat = rbind(par_mat, new_pars)
          lambda_vec = c(lambda_vec, new_lambda)
          post_vec = c(post_vec, new_post)
          accepted[i] = 1
          lock=FALSE
        }
        else{
          par_mat = rbind(par_mat, old_pars)
          lambda_vec = c(lambda_vec, tail(lambda_vec,1))
          post_vec = c(post_vec, tail(post_vec, 1))
        }
      }
      
    }
    else{
      par_mat = rbind(par_mat, old_pars)
      lambda_vec = c(lambda_vec, tail(lambda_vec,1))
      post_vec = c(post_vec, tail(post_vec, 1))
    }
  }
  par_mat = as.data.frame(par_mat)
  names(par_mat) = c('a','eps','k0','b')
  return(list(pars = par_mat, lambdas = lambda_vec, post = post_vec, sig=sig))
}
recover_params = function(a,eps,k0,b,stop_at=1e4,quiet=T,network.size=1e4,N=1e5,m=1){
  message('Simulating ... ')
  degs = samplepolylin(network.size, a, eps, b, k0,m=m)
  dat = twbfn::deg_count(degs)
  message('Fitting ...')
  out = polylin_mcmc(N, dat,init = c(a=1, eps=1,k0=10,b=1),
                     stop_at = stop_at,sig = diag(c(0.0001,0.0001,0.0001)),k0_jump = 1,
                     update_period = 1e2,quiet=quiet)
  return(list(output = out,degs  = degs))
}

counts_to_degs = function(dat){
  return(rep(dat[,1],dat[,2]))
}
fit_model = function(dat,N,burn.in=5e3,thin.by=5,trunc.at = min(dat[,1]),interval.size = 0.95,k0_jump=2,quiet=T){
  dat = dat[dat[,1]>=trunc.at,]
  out = polylin_mcmc(N, dat,quiet=quiet,burn.in = burn.in,k0_jump = k0_jump)
  smps = out$pars[seq(burn.in*1.2, nrow(out$pars),by=thin.by),]
  names(dat) = c('x','count')
  x = dat$x
  y_mat = numeric(length(x))
  pref_mat = numeric(length(x))
  
  pref = function(x,a,eps, k0, b){
    return(ifelse(x<k0, (x)^a + eps, (k0)^a  +eps + b*(x-k0)))
  }
  
  for(i in 1:nrow(smps)){
    est = as.numeric(smps[i,])
    lambda = find_lambda2(polylin(est[1], est[2]),est[4],est[3])
    if(!is.null(lambda)){
      ls = c(ls,lambda)
      y_mat = rbind(y_mat,S(x, polylin(est[1], est[2]),lambda, est[3],est[4])/S(min(dat[,1])-1, polylin(est[1], est[2]),lambda, est[3],est[4]))
      pref_mat = rbind(pref_mat,pref(x,est[1], est[2], est[3], est[4])/sum(dat[,2]*pref(x,est[1], est[2], est[3], est[4])))
    }
  }
  y_mat =y_mat[-1,]
  pref_mat = pref_mat[-1,]
  
  probs = c((1-interval.size)/2,1-(1-interval.size)/2)
  
  pref_CI = apply(pref_mat, 2, quantile, prob=probs,na.rm=T) 
  pref_50 = apply(pref_mat, 2, median,na.rm=T) 
  y_CI = apply(y_mat, 2, quantile, prob=probs,na.rm=T) 
  y_50 = apply(y_mat, 2, median,na.rm=T)
  return(list(mcmc = out, smps=smps, PA = list(CI = pref_CI, est = pref_50), surv = list(CI=y_CI, est=y_50),dat=dat))
}


recover_params = function(a, eps, b, k0,network.size=1e4,n.iter=3e4,quiet=T){
  degs = samplepolylin(network.size, a, eps, b, k0)
  dat = twbfn::deg_count(degs)
  out = fit_model(dat, n.iter, burn.in=n.iter/4, thin.by=5,quiet=quiet)
  out$smps$true_a = a
  out$smps$true_eps = eps
  out$smps$true_b = b
  out$smps$true_k0 = k0
  return(list(mcmc = out, degs = degs))
}


logdiff <- function(l1, l2) { l1 + VGAM::log1mexp(l1-l2); }
digp_log = function(x, shape, scale, threshold){
  out=x
  v=threshold
  sig=scale
  xi=shape
  p1 = pmax(x[x>v]*0, 1 + xi * (x[x>v]+1-v)/(sig))
  p2 = pmax(x[x>v]*0, 1 + xi * (x[x>v]-v)/(sig))
  pows2 = pows1 = x[x>v]*0 + 1
  pows1[p1>0] = -1/xi
  pows2[p2>0] = -1/xi
  # out[x>v] =  ( p2^pows2 - p1^pows1)
  out[x>v] = logdiff(p2^pows2,p1^pows1)
  out[x<=v] = 0
  return(out)
}

igp_ll_neg = function(pars,dat,threshold=0){
  shape = pars[1]
  scale = pars[2]
  if(shape==0){
    return(1e9)
  }
  if(shape<0 & max(dat[,1]+1)<(threshold - scale/shape)){
    return(1e9)
  }
  out = -sum(dat[,2]*log(digp(dat[,1], shape, scale, threshold)))
  if(is.infinite(out)){
    return(1e9)
  }else{
    return(out)
  }
}

fit.igp = function(dat, threshold=0){
  return(optim(c(0.5, 1), fn = igp_ll_neg,dat=dat[dat[,1]>threshold,], threshold=threshold,
               lower=c(-Inf,0.001), upper=c(Inf, Inf),method='L-BFGS-B'))
}














