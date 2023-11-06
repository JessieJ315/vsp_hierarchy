# essential helper functions 

u2po = function(U){
  # converts the U-matrix to partial order (transitive reduction)
  ## to add covariates (!)
  ## to adapt to varied k (!!)
  Z = invcdf_gumbel(pnorm(U))
  PO= transitive.reduction(order2partial(latent2order(Z),n=NUM_ACTORS))
  return(PO)
}

# distributions

## priors
# rho ~ Beta(1,1/6)
# tau ~ Gamma(shape,rate=1)
# theta ~ Gamma(shape,rate=1)

invcdf_gumbel = function(v) return(-log(-log(v)))

pU0 <- function(U,Sigma,log=TRUE){
  if (isTRUE(log)){
    return(sum(apply(U,1,dmvnorm,log=TRUE,sigma=Sigma)))
  } else {
    return(prod(apply(U,1,dmvnorm,log=FALSE,sigma=Sigma)))
  }
}

pU <- function(U,U0,tau,Sigma,log=TRUE){
  if (isTRUE(log)){
    log_probs=sapply(U,function(u) sum(apply(u-tau*U0,1,dmvnorm,log=TRUE,sigma=(1-tau^2)*Sigma)))
    return(sum(log_probs))
  } else {
    probs = sapply(U,function(u) prod(apply(u-tau*U0,1,dmvnorm,log=FALSE,sigma=(1-tau^2)*Sigma)))
    return(prod(probs))
  }
}

loglik <- function(U,theta,Y){
  # check!!!
  trs = lapply(U, u2po)
  tcs = lapply(trs, transitive.closure, mat=TRUE,loops=FALSE)
  potree = lapply(tcs, po2tree)
  is_vsp = sapply(potree,'[[','is_vsp')
  log_probs = c()
  if (any(is_vsp)){
    vsps = which(is_vsp)
    Y_vsp = Y[sapply(Y,'[[','assessor') %in% vsps]
    trees = lapply(potree, '[[', 'tree')
    nles = sapply(trees, nle.tree)
    log_probs = c(log_probs, sapply(Y_vsp, function(y) log(sum(sapply(allLE(tcs[[y$assessor]]),d.kendall,theta=theta,y=y$order,log=FALSE)))))
  }
  if (prod(is_vsp)==0){
    non_vsps = which(!is_vsp)
    Y_nonvsp = Y[sapply(Y,'[[','assessor') %in% non_vsps]
    log_probs = c(log_probs, sapply(Y_nonvsp, function(y) mallows(trs[y$assessor],y$order,theta)$lik))
  }
  return(sum(log_probs))
}

# data structure helper functions

initialisation = function(num_actors, num_assessors){
  # Simulates the initial states to each parameter from their priors.
  rho=rbeta(1,1,RHO_HYPERPAMA)
  tau=rgamma(1,shape=TAU_HYPERPAMA,rate=1)
  theta=rgamma(1,shape=THETA_HYPERPAMA,rate=1)
  Sigma=matrix(rho,nrow=2,ncol=2); diag(Sigma)=1
  U0 = rmvnorm(num_actors,sigma=Sigma)
  U_error = lapply(rep(num_actors,num_assessors),function(n) rmvnorm(n,sigma=(1-tau^2)*Sigma))
  U = lapply(U_error, function(ue) tau*U0+ue)
  return(list(rho=rho,tau=tau,theta=theta,U0=U0,U=U))
}

data_simulation = function(num_actors, num_assessors, num_orders){
  
  if(length(num_orders)!=num_assessors) stop("The number of assessors and the number of different orders doesn't match! ")
  
  initial_states = initialisation(num_actors, num_assessors)
  U = initial_states$U
  PO0 = u2po(initial_states$U0)
  trs = lapply(U, u2po)
  tcs = lapply(trs, transitive.closure, mat=TRUE, loops=FALSE)
  
  N = sum(num_orders)
  assessors = rep(1:3, times=num_orders)
  orders = lapply(assessors,function(i) unifLE(tcs[[i]]))
  
  Y = mapply(function(assessor, order) list(assessor=assessor, order=order), assessors, orders, SIMPLIFY=FALSE)
  return(list(PO0=PO0, POs=trs, Y=Y))
}
