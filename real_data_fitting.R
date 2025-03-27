source('funcs.R')
library(ggplot2)
data_raw = read.csv('degs.csv')
nms = unique(data_raw$name)
selected = c(2,3,5,6,7,8,9,11,12,13,17,22)
fits = list()
plots = list()
for(i in 1:length(selected)){
  print(i)
  fits[[i]] = fit_model(data_raw[data_raw$name==nms[selected[i]],2:3],N=5e4,quiet=F,
                        trunc.at = 5,k0_jump = 1)
}


saveRDS(fits, 'real_mcmc.rds')


for(i in 1:length(fits)){
  plots[[i]] = ggplot() + scale_x_log10() + scale_y_log10() +
    geom_point(data = twbfn::deg_surv(counts_to_degs(fits[[i]]$dat)),aes(x=degree, y=surv))+
    geom_line(aes(x=!!fits[[i]]$dat[,1],y=!!fits[[i]]$surv$est),colour='red') + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[1,]),linetype=2, colour='red') + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[2,]),,linetype=2,colour='red')+
    ggtitle(nms[selected[i]])
}


real_survs = ggpubr::ggarrange(plotlist = plots)

saveRDS(real_survs,'real_survs.rds')




real_survs
# igpd --------------------------------------------------------------------


fits = readRDS('real_mcmc.rds')

gpfit = mev::fit.gpd(counts_to_degs(fits[[1]]$dat),threshold=6)

scale = gpfit$estimate[1]
shape = gpfit$estimate[2]

x = fits[[1]]$dat[,1]
y = mev::pgp(x,loc=6,scale=scale,shape=shape,lower.tail = F)
degs = counts_to_degs(fits[[1]]$dat)

plot(x,y,log='xy',type='l')
points(twbfn::deg_surv(degs[degs>6]))

# -------------------------------------------------------------------------

data_raw = read.csv('degs.csv')
nms = unique(data_raw$name)
selected = c(2,3,5,6,7,8,9,11,12,13,17,22)

nm <- nms[selected[1]]

dat = data_raw[data_raw$name==nm,][,2:3]


degs = counts_to_degs(dat)


plot(twbfn::deg_surv(degs), log='xy')


fit = fit.gpd(degs, threshold=15)


# -------------------------------------------------------------------------


fits = readRDS('real_mcmc.rds')







# -------------------------------------------------------------------------



plts = list()

for(j in 1:length(fits)){
  print(j)
  fit = fits[[j]]
  shapes<- numeric()
  xs <- numeric()
  for(i in 1:(nrow(fit$dat)-3)){
    if(fit$dat[i,1]>max(fit$mcmc$pars$k0[seq(2e4, 5e4,by=10)])){
      break;
    }
    igp = mev::fit.gpd(counts_to_degs(fit$dat),threshold=fit$dat[i,1])
    shapes = c(shapes,igp$estimate[2])
    xs = c(xs,fit$dat[i,1])
  }
  plts[[j]] = ggplot() +
    geom_point(aes(x=!!fit$mcmc$pars$k0[seq(2e4, 5e4,by=10)],y=!!(fit$mcmc$pars$b[seq(2e4, 5e4,by=10)]/fit$mcmc$lambdas[seq(2e4, 5e4,by=10)])),colour='blue', alpha=0.0025) + 
    geom_point(aes(x=!!xs[xs<=max(fit$mcmc$pars$k0[seq(2e4, 5e4,by=10)])],y=!!shapes[xs<=max(fit$mcmc$pars$k0[seq(2e4, 5e4,by=10)])]), colour='red') +
    ylab('Shape') + xlab('k0') + ggtitle(ggtitle(nms[selected[j]]))
}


igp_plt = ggpubr::ggarrange(plotlist = plts,nrow=3,ncol=4)


saveRDS(igp_plt, 'gp_plt.rds')


plot(igp_plt)

df = data.frame(k0 = fit$mcmc$pars$k0,shape = fit$mcmc$pars$b/fit$mcmc$lambdas)

df$real = shapes[match(df$k0,xs)]


plt = readRDS('real_survs.rds')






























