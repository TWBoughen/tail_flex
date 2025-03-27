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




fits = readRDS('real_mcmc.rds')


plots = list()
for(i in 1:length(fits)){
  plots[[i]] = ggplot() + scale_x_log10() + scale_y_log10() +
    geom_point(data = twbfn::deg_surv(counts_to_degs(fits[[i]]$dat)),aes(x=degree, y=surv))+
    geom_line(aes(x=!!fits[[i]]$dat[,1],y=!!fits[[i]]$surv$est),colour='red') + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[1,]),linetype=2, colour='red') + 
    geom_line(aes(x=!!fits[[i]]$dat[,1],y= !!fits[[i]]$surv$CI[2,]),linetype=2,colour='red')+
    geom_vline(aes(xintercept = !!fits[[i]]$smps$k0), colour='blue',alpha=0.1)+
    ggtitle(nms[selected[i]])
}


real_survs = ggpubr::ggarrange(plotlist = plots)

real_survs


saveRDS(real_survs,'real_survs.rds')







max()







fits[[9]]$surv$CI[2,]






