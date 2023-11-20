# this script analysis the MCMC output.

chain_name = 'test'

################################################################################
dir='~/Desktop/vsp-hierarchy/output/'
load(paste0(dir,chain_name,'.RData'))

source("~/Desktop/vsp-hierarchy/code/vsp_hierarchy_func.R")
library(igraph)
library(Rgraphviz)
library(graph)
library(ggplot2)
library(ggpubr)
library(kableExtra)
library(RColorBrewer)
library(reshape)

N = sum(output$loglk<0);N
num_actors=nrow(output$U0[[1]])
num_assessors=length(output$U[[1]])
burn_in = 2000

# traceplots
par(mfrow=c(2,1))
## paramters - main
plot(output$loglkd[1:N],type='l',main='loglikelihood',ylab='loglikelihood')
abline(v = burn_in, col='red')
plot(output$rho[1:N],type='l',main='latent matrix depth control parameter',ylab='rho')
abline(v = burn_in, col='red')
plot(output$tau[1:N],type='l',main='latent matrix dispersion parameter',ylab='tau')
abline(v=burn_in, col='red')
## paramters - noise model
plot(output$p[1:N],type='l',main='error probability for the queue-jumping noise model',ylab='p')
abline(v=burn_in, col='red')
plot(output$theta[1:N],type='l',main="dispersion parameter for the Mallow's noise model",ylab='theta')
abline(v=burn_in, col='red')
## latent matrices
U0_means=do.call(rbind,lapply(output$U0[1:N], function(m){apply(m,1,mean)}))
U0_means=data.frame(index=1:N,U0_means)
colnames(U0_means)=c('index',paste0('actor',1:num_actors))

df=melt(U0_means,id='index')
ggplot(df,aes(x=index,y=value,group=variable,color=variable))+
  geom_line()+
  ylab('latent path mean')+
  ggtitle('latent preference weight matrix for the global partial order')+
  theme_classic()

for (i in 1:num_assessors){
  U_assessor=lapply(output$U[1:N], function(U){U[[i]]})
  U_means=do.call(rbind,lapply(U_assessor, function(m){apply(m,1,mean)}))
  U_means=data.frame(index=1:N,U_means)
  colnames(U_means)=c('index',paste0('actor',1:num_actors))
  
  df=melt(U_means,id='index')
  g=ggplot(df,aes(x=index,y=value,group=variable,color=variable))+
    geom_line()+
    ylab('latent path mean')+
    ggtitle(paste0('latent preference weight matrix assessor ',i))+
    theme_classic()
  print(g)
}

# partial orders

## global vsp
global_pos=lapply(output$U0[1:N],u2po)
global_pos=lapply(global_pos,transitive.closure,mat=TRUE,loops=FALSE)
PO_unlist=table(sapply(global_pos[burn_in:N],paste,collapse = " "))
### mode
sort(PO_unlist,decreasing = TRUE)[1:5]
PO_mp=str2po(names(which(PO_unlist==max(PO_unlist))),num_actors)
par(mfrow=c(1,1))
showDAG(transitive.reduction(PO_mp))
showDAG(transitive.reduction(str2po(names(sort(PO_unlist,decreasing = TRUE)[2]),num_actors)))
### concensus order
concensus_po = Reduce('+',global_pos[burn_in:N])/length(global_pos[burn_in:N]) # po_posterior transitive closure
showDAGcon(concensus_po,threshold =.3)

## individual vsps
for (i in 1:num_assessors){
  U_assessor=lapply(output$U[1:N], function(U){U[[i]]})
  individual_pos=lapply(U_assessor,u2po)
  PO_unlist=table(sapply(individual_pos[burn_in:N],paste,collapse = " "))
  ### mode
  sort(PO_unlist,decreasing = TRUE)[1]
  PO_mp=str2po(names(which(PO_unlist==max(PO_unlist))),num_actors)
  par(mfrow=c(1,1))
  showDAG(transitive.reduction(PO_mp))
  showDAG(transitive.reduction(str2po(names(sort(PO_unlist,decreasing = TRUE)[2]),num_actors)))
  ### concensus order
  concensus_po = Reduce('+',individual_pos[burn_in:N])/length(individual_pos[burn_in:N]) # po_posterior transitive closure
  showDAGcon(concensus_po,threshold =.3)
}

# posterior distributions
## rho
ggplot(melt(data.frame(prior=rbeta(N-burn_in,1,RHO_HYPERPAMA),posterior = output$rho[(burn_in+1):N])),aes(x=value,y=..density..,group=variable,color=variable,linetype=variable))+
  geom_density(alpha=0.4,adjust = 3)+
  labs(title="Prior/Posterior Distributiosn for rho",
       y="", x = "rho")+
  scale_color_manual(values=c("#999999", "#0072B2"))+
  scale_linetype_manual(values=c('dashed','solid'))+
  theme_classic()+
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title=element_blank(),text=element_text(family='Times New Roman'))
mean(output$rho[burn_in:N])
## tau
ggplot(melt(data.frame(prior=1/(1+exp(-rnorm(N-burn_in,sd=T_HYPERPAMA))),posterior = output$tau[(burn_in+1):N])),aes(x=value,y=..density..,group=variable,color=variable,linetype=variable))+
  geom_density(alpha=0.4,adjust = 3)+
  labs(title="Prior/Posterior distributiosn for tau",
       y="", x = "tau")+
  scale_color_manual(values=c("#999999", "#0072B2"))+
  scale_linetype_manual(values=c('dashed','solid'))+
  theme_classic()+
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title=element_blank(),text=element_text(family='Times New Roman'))
mean(output$tau[(burn_in+1):N])
## p
ggplot(melt(data.frame(prior=1/(1+exp(-rnorm(N-burn_in,sd=R_HYPERPAMA))),posterior = output$p[(burn_in+1):N])),aes(x=value,y=..density..,group=variable,color=variable,linetype=variable))+
  geom_density(alpha=0.4,adjust = 3)+
  labs(title="Prior/Posterior distributiosn for Error Probability p",y="", x = "p")+
  scale_color_manual(values=c("#999999", "#0072B2"))+
  scale_linetype_manual(values=c('dashed','solid'))+
  theme_classic()+
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title=element_blank(),text=element_text(family='Times New Roman'))
mean(output$p[(burn_in+1):N])
## theta
ggplot(melt(data.frame(prior=rgamma(N-burn_in,shape=THETA_HYPERPAMA,rate=1),posterior = output$theta[(burn_in+1):N])),aes(x=value,y=..density..,group=variable,color=variable,linetype=variable))+
  geom_density(alpha=0.4,adjust = 3)+
  labs(title="Prior/Posterior distributiosn for theta",
       y="", x = "theta")+
  scale_color_manual(values=c("#999999", "#0072B2"))+
  scale_linetype_manual(values=c('dashed','solid'))+
  theme_classic()+
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title=element_blank(),text=element_text(family='Times New Roman'))
mean(Result$theta[(burn_in+1):N])



