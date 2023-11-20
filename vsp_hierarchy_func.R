# This script contains helper functions for the hierarchical vsp model, including
# - essential helper functions 
# - relevant distributions
# - data structure helper functions
# - visualisation helper functions

# required libaries
library(foreach)
library(doParallel)
library(truncnorm)
library(mvtnorm)
library(condMVNorm)

# essential helper functions 
po2tree<-function(h,tree=vector('list',0)) {
  # tests whether a partial order is a vsp, and generates a binary decomposition 
  # tree from a partial order if so.
  # 
  # Parameters
  # ----------
  # h: matrix. 
  #   The transitive closure of the partial order. 
  # tree: list, default tree = vector('list',0). 
  #   The working tree (vsp). 
  #
  # Returns
  # -------
  #   list, with attributes
  #     h: matrix.  
  #       The unstandardised partial order (transitive closure).
  #     tree: list. 
  #       The binary decomposition tree corresponding to the input partial order. 
  #     is_vsp: bool.
  #       Whether the partial order is a vsp.
  # 
  # example
  # -------
  # N=6; actors=1:N; tree=rvsp(actors,0.5); h=tree2po(tree)
  # par(mfrow=c(1,2)); showDAG(transitive.reduction(h),edge.arrow.size=3/N); 
  # showTREE(tree,edge.arrow.size=3/N)
  n=dim(h)[1]
  is_vsp=TRUE
  if (length(tree)==0) {
    rn=rownames(h)
    tree=lapply(1:n,function(x) node(actor=as.numeric(rn[x])))
    rn=1:n
    rownames(h)<-colnames(h)<-rn 
  }
  rn=sapply(rownames(h),as.numeric)
  
  new.node.i=length(tree)+1
  if (n>2) {
    i=0; finished_1=FALSE
    while (!finished_1) {
      i=i+1
      if (i==n) {return(list(h=h,tree=tree,is_vsp=FALSE))}
      j=i; finished_2=FALSE
      while (!finished_2 & j < n) {
        #ij can be a cherry if they have same relns to all other nodes
        j=j+1
        finished_1=finished_2=is.cherry(h,i,j)
      }
    }
  } else {i=1;j=2}
  
  i.i=rn[i]
  j.i=rn[j]
  tree[[i.i]]$parent=new.node.i
  tree[[j.i]]$parent=new.node.i
  if (h[i,j]==1) {
    tree[[new.node.i]]=node(child=c(i.i,j.i),type='S',order='');
    tree[[i.i]]$order='+'
    tree[[j.i]]$order='-'
  }
  if (h[j,i]==1) {
    tree[[new.node.i]]=node(child=c(j.i,i.i),type='S',order='');
    tree[[j.i]]$order='+'
    tree[[i.i]]$order='-'
  }
  if (h[i,j]==0 & h[j,i]==0 ) {
    tree[[new.node.i]]=node(child=c(j.i,i.i),type='P',order='');
  }
  
  h<-h[-j,-j,drop=FALSE]; 
  rn[i]<-new.node.i; rn<-rn[-j] 
  rownames(h)<-colnames(h)<-rn
  
  if (dim(h)[1]>1) {
    out=po2tree(h,tree)
    if (!out$is_vsp) {return(list(h=out$h,tree=out$tree,is_vsp=out$is_vsp))}
    h=out$h
    tree=out$tree
  }
  
  # add children counts
  v=which(sapply(tree,function(x) is.root(x)))
  cc=get.child.count(tree,v)
  n=dim(cc)[2]
  for (i in 1:n) tree[[cc[1,i]]]$nc=cc[2,i]
  
  return(list(h=h,tree=tree,is_vsp=is_vsp))
}

u2po = function(U){
  # converts the U-matrix to partial order (transitive closure).
  ## to add covariates (!)
  ## to adapt to varied k (!!)
  Z =invcdf_gumbel(pnorm(U))
  tc=transitive.closure(order2partial(latent2order(Z),n=NUM_ACTORS),
                         mat=TRUE,loops=FALSE)
  return(tc)
}

# relevant distributions
invcdf_gumbel = function(v) return(-log(-log(v)))

pU0 <- function(U,Sigma,log=TRUE){
  # The distribution over U0 given Sigma. 
  if (isTRUE(log)){
    return(sum(apply(U,1,dmvnorm,log=TRUE,sigma=Sigma)))
  } else {
    return(prod(apply(U,1,dmvnorm,log=FALSE,sigma=Sigma)))
  }
}

pU <- function(U,U0,tau,Sigma,log=TRUE){
  # The distribution over U given U0, tau and Sigma.
  if (isTRUE(log)){
    log_probs=sapply(U,function(u) sum(apply(u-tau*U0,1,dmvnorm,log=TRUE,
                                             sigma=(1-tau^2)*Sigma)))
    return(sum(log_probs))
  } else {
    probs = sapply(U,function(u) prod(apply(u-tau*U0,1,dmvnorm,log=FALSE,
                                            sigma=(1-tau^2)*Sigma)))
    return(prod(probs))
  }
}

loglikQJ <- function(trees,Y,p,theta=NULL,no_cores=2,phi=0,model='lkdup'){
  
  # The log-likelihood for the queue-jumping observation model. This function uses
  # parallelism to speed up calculation.
  # 
  # Parameters
  # ----------
  # trees: list. 
  #   The binary-decomposition-trees for each U-matrix, presented of the same 
  #     order as Us.
  # Y: list
  #   The nested data list. Each data entry contains attributes 'assessor' and 
  #     'order'. 
  # p: float in [0,1].
  #   The error probability for the queue-jumping noise model. 
  # theta: float, default = NULL.
  #   Placeholder for the dispersion parameter in the Mallow's model. Not used
  #     in this function. 
  # no_cores: int, default=2.
  #   The number of cores to use for parallel computing.
  # phi: float in [0,1], default=0.
  #   The jumping direction probability in the bidirectional queue-jumping model.
  # model: string, in 'lkdup', 'lkddown' and 'bi-directn'.
  #   The queue-jumping to use. 'model' takes value 'lkdup' - queue-jumping-down, 
  #     'lkddown' - queue-jumping-up and 'bi-directn' for bi-directional queue-
  #     jumping.
  #
  # Returns
  # -------
  #   float.
  #     The log-likelihood on the data list based on a specified qj model.
  
  Q.LEProb.vsp <- function(tree,le,p,q=0) {
    
    #if tree is a vsp and le is one list then calculate the likeligood for the list
    #given the vsp-PO tree. Here p (err prob) and q (prob choose top down at a given le entry insertion)
    #adapts to the suborder so le and actors.in.tree have same content
    
    tree = sub.tree(tree,le)
    
    n=length(le)
    if (n==1) return(1)
    
    leaf.nodes=which(sapply(tree,is.leaf))
    actors.in.tree=sapply(leaf.nodes, function(x) tree[[x]]$actor)
    
    ntc=NA
    if (q>0) {
      top.i=leaf.nodes[which(actors.in.tree==le[1])]
      if (length(top.i)!=1) stop('err in QP.LEProb.vsp length(top.i)!=1')
      le.not=le[-1]
      tree.not=delete(tree,top.i)
      prob.top=p/n
      if (is.top.bot(tree,top.i,'top')) {
        ntc=nle.tree(tree); 
        prob.top=prob.top+(1-p)*nle.tree(tree.not)/ntc
      }
      top.fac=q*prob.top*QP.LEProb.vsp(tree.not,le.not,p,q)
    } else {top.fac=0}
    
    if (q<1) {
      bot.i=leaf.nodes[which(actors.in.tree==le[n])]
      if (length(bot.i)!=1) stop('err in QP.LEProb.vsp length(bot.i)!=1')
      le.nob=le[-n]
      tree.nob=delete(tree,bot.i)
      prob.bot=p/n
      if (is.top.bot(tree,bot.i,'bot')) {
        if (is.na(ntc)) {ntc=nle.tree(tree)}
        prob.bot=prob.bot+(1-p)*nle.tree(tree.nob)/ntc
      }
      bot.fac=(1-q)*prob.bot*QP.LEProb.vsp(tree.nob,le.nob,p,q)
    } else {bot.fac=0}
    
    return(top.fac+bot.fac)
  }
  
  registerDoParallel(cores=no_cores)
  export_functions = c('nle.tree','sub.tree','find.parents','find.actor','is.root',
                       'is.top.bot','is.zombie','QP.LEProb.vsp','delete','is.leaf')
  
  if (model=='lkddown') {
    llkda = foreach(le=Y,.combine=c,.export=export_functions) %dopar% 
      Q.LEProb.vsp(tree=trees[[le$assessor]],le=le$order,p=p,q=1)
    }
  if (model=='lkdup') {
    llkda = foreach(le=Y,.combine=c,.export=export_functions) %dopar% 
      Q.LEProb.vsp(tree=trees[[le$assessor]],le=le$order,p=p,q=0)
    }
  if (model=='bi-directn') {
    q = 1/(1+exp(-phi))
    llkda = foreach(le=Y,.combine=c,.export=export_functions) %dopar% 
      Q.LEProb.vsp(tree=trees[[le$assessor]],le=le$order,p=p,q=q)
  }
  
  stopImplicitCluster()
  
  return(sum(log(llkda)))
}

loglikM <- function(trees,Y,p,theta){
  # the mallow's log-likelihood model - need to update.
  tcs = lapply(U, u2po)
  potree = lapply(tcs, po2tree)
  is_vsp = sapply(potree,'[[','is_vsp')
  log_probs = c()
  if (any(is_vsp)){
    vsps = which(is_vsp)
    Y_vsp = Y[sapply(Y,'[[','assessor') %in% vsps]
    trees = lapply(potree, '[[', 'tree')
    nles = sapply(trees, nle.tree)
    log_probs = c(log_probs, sapply(Y_vsp, function(y) 
      log(sum(sapply(allLE(tcs[[y$assessor]]),d.kendall,theta=theta,y=y$order,
                     log=FALSE)))))
  }
  if (prod(is_vsp)==0){
    non_vsps = which(!is_vsp)
    Y_nonvsp = Y[sapply(Y,'[[','assessor') %in% non_vsps]
    log_probs = c(log_probs, sapply(Y_nonvsp, function(y) 
      mallows(trs[y$assessor],y$order,theta)$lik))
  }
  return(sum(log_probs))
}

# data structure helper functions

initialisation = function(num_actors, num_assessors){
  # Simulates the initial states to each parameter from their priors.
  # 
  # Parameters
  # ----------
  # num_actors: int. 
  #   The number of actors. 
  # num_assessors: int. 
  #   The number of assessors. 
  # 
  # Returns
  # -------
  #   list, with attributes
  #     rho: float, in [0,1].
  #       The depth control parameter for the latent preference weight matrices.
  #     tau:   float=0.5.
  #       The dispersion parameter for latent matrices U. 
  #     theta: float, in [0,inf).
  #       The dispersion parameter for the Mallow's noise model. 
  #     U0:    matrix.
  #       The latent preference weight matrix for the global partial order.
  #     U:     list [matrix].
  #       The latent preference weight matrices for the partial orders 
  #         corresponding to each assessor.
  #     p:     float=0.5.
  #       The error probability for the queue-jumping noise model. 
  rho=rbeta(1,1,RHO_HYPERPAMA)
  Sigma=matrix(rho,nrow=2,ncol=2); diag(Sigma)=1
  U0 = rmvnorm(num_actors,sigma=Sigma)
  while (!po2tree(u2po(U0))$is_vsp){
    rho=rbeta(1,1,RHO_HYPERPAMA)
    Sigma=matrix(rho,nrow=2,ncol=2); diag(Sigma)=1
    U0 = rmvnorm(num_actors,sigma=Sigma)
  }
  tau=0.5
  U_error = lapply(rep(num_actors,num_assessors),function(n) rmvnorm(n,sigma=(1-tau^2)*Sigma))
  U = lapply(U_error, function(ue) tau*U0+ue)
  while (!all(sapply(lapply(U,u2po),function(po){po2tree(po)$is_vsp}))){
    U_error = lapply(rep(num_actors,num_assessors),function(n) rmvnorm(n,sigma=(1-tau^2)*Sigma))
    U = lapply(U_error, function(ue) tau*U0+ue)
  }
  theta=rgamma(1,shape=THETA_HYPERPAMA,rate=1)
  return(list(rho=rho,tau=tau,theta=theta,U0=U0,U=U,p=0.5))
}

data_simulation = function(num_actors, num_assessors, num_orders,U0=NULL,U=NULL){
  # simulates a data list with entries 'assessor' and 'order'.
  # 
  # Parameters
  # ----------
  # num_actors: int. 
  #   The number of actors. 
  # num_assessors: int. 
  #   The number of assessors. 
  # num_orders: vector. 
  #   The number of orders from each assessor. 
  # U0: matrix, default=NULL.
  #   The latent preference weight matrix for the global partial order. If null,
  #     simulate from prior. 
  # U: list, default=NULL.
  #   The latent preference weight matrices for the partial orders corresponding 
  #     to each assessor. If null, simulate from prior.
  # 
  # Returns
  # -------
  #   list.
  #     The data list - a nested list each with entries 'assessor' and 'order'.
  if(length(num_orders)!=num_assessors) 
    stop("The number of assessors and the number of different orders doesn't match! ")
  
  if (is.null(U)){
    initial_states = initialisation(num_actors, num_assessors)
    U = initial_states$U
    U0 = initial_states$U0
  }
  PO0 = u2po(initial_states$U0)
  trs = lapply(U, u2po)
  tcs = lapply(trs, transitive.closure, mat=TRUE, loops=FALSE)
  
  N = sum(num_orders)
  assessors = rep(1:3, times=num_orders)
  orders = lapply(assessors,function(i) unifLE(tcs[[i]]))
  
  Y = mapply(function(assessor, order) list(assessor=assessor, order=order), 
             assessors, orders, SIMPLIFY=FALSE)
  return(list(PO0=PO0, POs=trs, Y=Y))
}

## helper function for visualization

str2po=function(x,n){
  # turns a string to an adjacency matrix
  # 
  # Parameters
  # ----------
  # x: string.
  #   The string of an adjacency matrix.
  # n: int. 
  #   The number of items. 
  # 
  # Returns
  # -------
  #   matrix.
  #     The adjacency matrix. 
  matrix(as.numeric(strsplit(x," ")[[1]]),nrow=n)
}

showDAGcon <- function(con_po,threshold=.5,label_threshold=.9,size = 10){
  # plots the concensus order given a concensus adjacency matrix.
  # 
  # Parameters
  # ----------
  # con_po: matrix.
  #   The weighted concensus adjacency matrix (CAM).
  # threshold: float in [0,1], default=0.5.
  #   The threshold to plot an edge. If the weight in CAM is less than threshold, 
  #     such edge is not plotted.
  # label_threshold: float in [0,1], default=0.9.
  #   Order relations with weights larger than label_threshold will be colored red. 
  # size: int, default=10.
  #   The node size.
  con_po1 = (con_po>threshold)+0.0
  con_po2 = (con_po>label_threshold)+0.0
  
  con_po1 = transitive.reduction(con_po1)
  con_po2 = transitive.reduction(con_po2)
  
  g = graph_from_adjacency_matrix(con_po1,mode ="directed")#as(m,"graphNEL")
  
  g_label = graph_from_adjacency_matrix(con_po2)
  el1 <- apply(get.edgelist(g), 1, paste, collapse="-")
  el2 <- apply(get.edgelist(g_label), 1, paste, collapse="-")
  E(g)$color <- ifelse(el1 %in% el2, "red", "darkgrey")
  
  # tkplot(g,layout=layout_with_sugiyama(g)$layout,vertex.size=8,
  #     edge.arrow.size=.5,vertex.label=NA)
  plot(g,layout=layout_with_sugiyama(g)$layout,vertex.size=size,
       edge.arrow.size=.2)#,vertex.label=NA)
}
