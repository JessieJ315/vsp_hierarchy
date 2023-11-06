
vsp_hierarchy = function(initial_states,chain.name,n.itr=100,n.record=NULL){
  # mcmc function
  n = nrow(Z)
  if(is.null(n.record)) {n.record = 2*n}
  # load initial states
  U0 = initial_states$U0
  U = initial_states$U
  rho = initial_states$rho
  tau = initial_states$tau
  theta = initial_states$theta
  Sigma=matrix(rho,nrow=2,ncol=2); diag(Sigma)=1
  # initiate storage lists/vectors
  U0_PATH = vector(mode= "list", length = n.itr/n.record)
  U_PATH = vector(mode="list", length = n.itr/n.record)
  RHO_PATH = vector(mode="numeric", length = n.itr/n.record)
  TAU_PATH = vector(mode= "numeric", length = n.itr/n.record)
  THETA_PATH = vector(mode= "numeric", length = (n.itr/n.record))
  LOGLKD_PATH = vector(mode="numeric", length = n.itr/n.record)
  
  Sigma = matrix(rho, 2, 2); diag(Sigma) = 1
  loglkd = loglik(U, theta, Y)
  
  for(itr in 1:n.itr){
    # update on U0
    U0_temp = U0
    c = sample(NUM_ACTORS,1); i = sample(2,1)
    U0_temp[c,i] = rcmvnorm(n=1, mean=U0[c,], sigma=Sigma, dependent.ind=i, given.ind=(1:2)[-i], X.given=U0[c,-i])
    log_accept_rate = pU0(U0_temp,Sigma)-pU0(U0,Sigma)+pU(U,U0=U0_temp,tau,Sigma)-pU(U,U0=U0,tau,Sigma)
    if (log_accept_rate > log(runif(1))) {U0=U0_temp}
    # update on U
    for (j in 1:NUM_ASSESSORS){
      U_temp=U
      u_temp = U[[j]]
      c = sample(NUM_ACTORS,1); i = sample(2,1)
      u_temp[c,i] = rcmvnorm(n=1, mean=u_temp[c,], sigma=Sigma, dependent.ind=i, given.ind=(1:2)[-i], X.given=u_temp[c,-i])
      U_temp[[j]]=u_temp
      loglkd_temp = loglik(U_temp,theta,Y)
      log_accept_rate = loglkd_temp - loglkd + pU(U_temp,U0,tau,Sigma) - pU(U,U0,tau,Sigma)
      if (log_accept_rate > log(runif(1))) {U=U_temp; loglkd=loglkd_temp}
    }
    # update on rho
    delta = runif(1, RHO_PROPOSAL, 1/RHO_PROPOSAL)
    rho_temp = 1 - delta*(1-rho)
    Sigma_temp = matrix(rho_temp,nrow=2,ncol=2); diag(Sigma_temp)=1
    log_accept_rate = dbeta(rho_temp,1,RHO_HYPERPAMA,log=TRUE)-dbeta(rho,1,RHO_HYPERPAMA,log=TRUE)+pU0(U0,Sigma_temp)-pU0(U0,Sigma)+pU(U,U0,tau,Sigma_temp)-pU(U,U0,tau,Sigma)-log(delta)
    if(log_accept_rate > log(runif(1))){rho=rho_temp;Sigma=Sigma_temp}
    # update on theta
    theta_temp = rnorm(1,theta,0.5)
    if(theta_temp>0){
      loglkd_temp = loglik(U,theta_temp,Y)
      log_accept_rate = loglkd_temp+dgamma(theta_temp,shape=THETA_HYPERPAMA,rate=1,log=TRUE)-loglkd-dgamma(theta,shape=THETA_HYPERPAMA,rate=1,log=TRUE)
      if(log_accept_rate > log(runif(1))){theta=theta_temp; loglkd=loglkd_temp}
    }
    # update on tau
    tau_temp = rnorm(1,tau,0.5)
    if(tau_temp>0){
      log_accept_rate=dgamma(tau_temp,shape=TAU_HYPERPAMA,rate=1,log=TRUE)-dgamma(tau,shape=TAU_HYPERPAMA,rate=1,log=TRUE)+pU(U,U0,tau_temp,Sigma)-pU(U,U0,tau,Sigma)
      if(log_accept_rate > log(runif(1))){tau=tau_temp}
    }
    if (itr %% n.record == 0){
      U0_PATH[[itr/n.record]] = U0
      U_PATH[[itr/n.record]] = U
      RHO_PATH[itr/n.record] = rho
      THETA_PATH[itr/n.record] = theta
      TAU_PATH[itr/n.record] = tau
      LOGLKD_PATH[itr/n.record] = loglkd
    }
    if (itr %% (n.itr/10) == 0){
      print(itr)
      Result = list(U0=U0_PATH,U=U_PATH,rho=RHO_PATH,theta=THETA_PATH,tau=TAU_PATH,loglkd=LOGLKD_PATH)
      save(Result, file=paste0(chain.index,".RData"))
    }
  }
  return(Result)
}
