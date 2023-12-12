# This script contains the main function for hierarchical vsp mcmc update.

vsp_hierarchy = function(initial_states,Y,chain_name,n_itr=100,n_record=NULL,
                         noise_model='queue-jumping'){
  
  # MCMC update for the Hierarchical-VSP model. 
  # This function performs Metropolis Hasting MCMC for the paramters of interest
  # for the hierarchical vsp model. The output is stored in the working directory.
  # 
  # Parameters
  # ----------
  # initial_states: list.
  #   The initial states of the model parameters, including:
  #
  #     U0:    matrix.
  #       The latent preference weight matrix for the global partial order.
  #     U:     list [matrix].
  #       The latent preference weight matrices for the partial orders 
  #         corresponding to each assessor.
  #     rho:   float, in [0,1].
  #       The depth control parameter for the latent preference weight matrices.
  #     tau:   float, in [0,1].
  #       The dispersion parameter for latent matrices U. When tau=0, U0 and Us 
  #         are independent. 
  #     p:     float, in [0,1].
  #       The error probability for the queue-jumping noise model. 
  #     theta: float, in [0,inf).
  #       The dispersion parameter for the Mallow's noise model. 
  # 
  # Y: list. 
  #   The nested data list. Each data entry contains attributes 'assessor' and 
  #     'order'. 
  # 
  # chain_name: string. 
  #   The name of the output file. 
  #
  # n_itr: int.
  #   The number of MCMC iterations. 
  #
  # n_record: int. 
  #   The step between each path record. We only record every n_record iterations. 
  #
  # noise_model: string.
  #   The choice of noise model, either 'queue-jumping' (vsp) for 'mallows' 
  #     (general po). 
  # 
  # Returns
  # -------
  #   List, with attributes:
  #     U0, U, rho, p, tau: same as above. 
  
  # load initial states
  U0 = initial_states$U0
  U = initial_states$U
  rho = initial_states$rho
  tau = initial_states$tau
  p = initial_states$p
  theta = initial_states$theta
  
  Sigma=matrix(rho,nrow=2,ncol=2); diag(Sigma)=1
  pos = lapply(U, u2po)
  trees = lapply(pos, function(po) po2tree(po)$tree)
  
  # set recording steps
  n = nrow(U0)
  if(is.null(n_record)) {n_record = 2*n}
  
  # define likelihood function
  if (noise_model=='queue-jumping'){loglik=loglikQJ}
  if (noise_model=='mallows'){loglik=loglikM}
  
  loglkd = loglik(trees,Y,p,theta)
  
  # initiate storage lists/vectors
  U0_PATH = vector(mode= "list", length = n_itr/n_record)
  U_PATH = vector(mode="list", length = n_itr/n_record)
  RHO_PATH = vector(mode="numeric", length = n_itr/n_record)
  TAU_PATH = vector(mode= "numeric", length = n_itr/n_record)
  THETA_PATH = vector(mode= "numeric", length = (n_itr/n_record)) # for MALLOWS
  P_PATH = vector(mode= "numeric", length = (n_itr/n_record)) # for QUEUE-JUMPING
  LOGLKD_PATH = vector(mode="numeric", length = n_itr/n_record)
  
  for(itr in 1:n_itr){
    # update on U0 and U
    U0_temp = U0
    c = sample(NUM_ACTORS,1)
    U0_temp_c = rmvnorm(n=1,sigma=Sigma)
    U0_temp[c,] = U0_temp_c
    U_temp = U
    trees_temp = trees
    allvsp = TRUE
    for (j in 1:NUM_ASSESSORS){
      u = U_temp[[j]]
      u[c,] = rmvnorm(n=1,mean=tau*U0_temp_c,sigma=(1-tau^2)*Sigma)
      potree_j = po2tree(u2po(u))
      if (!potree_j$is_vsp) {
        allvsp=FALSE; break
      } else {
        trees_temp[[j]]=potree_j$tree
        U_temp[[j]] = u
      }
    }
    if (allvsp){
      loglkd_temp = loglik(trees_temp,Y,p,theta)
      log_accept_rate = loglkd_temp - loglkd
      if (log_accept_rate > log(runif(1))) {
        U0=U0_temp; U=U_temp
        trees=trees_temp; loglkd=loglkd_temp
        }
    }
    # update on rho
    # prior: rho ~ Beta(1,RHO_HYPERPAMA)
    delta = runif(1, RHO_PROPOSAL, 1/RHO_PROPOSAL)
    rho_temp = 1 - delta*(1-rho)
    Sigma_temp = matrix(rho_temp,nrow=2,ncol=2); diag(Sigma_temp)=1
    log_accept_rate = dbeta(rho_temp,1,RHO_HYPERPAMA,log=TRUE)-dbeta(rho,1,RHO_HYPERPAMA,log=TRUE)+
      pU0(U0,Sigma_temp)-pU0(U0,Sigma)+pU(U,U0,tau,Sigma_temp)-pU(U,U0,tau,Sigma)-log(delta)
    if(log_accept_rate > log(runif(1))){rho=rho_temp;Sigma=Sigma_temp}
    # update on theta - MALLOWS
    # prior: theta ~ gamma(shape=THETA_HYPERPAMA,rate=1)
    if (noise_model=='mallows'){
      theta_temp = rnorm(1,theta,0.5)
      if(theta_temp>0){
        loglkd_temp = loglik(U,theta_temp,Y)
        log_accept_rate = loglkd_temp+dgamma(theta_temp,shape=THETA_HYPERPAMA,rate=1,log=TRUE)-
          loglkd-dgamma(theta,shape=THETA_HYPERPAMA,rate=1,log=TRUE)
        if(log_accept_rate > log(runif(1))){theta=theta_temp; loglkd=loglkd_temp}
      }
    }
    # update on p - QUEUE-JUMPING - WORKS! 
    # prior: r = log(p/(1-p)) ~ Normal(0, R_HYPERPAMA)
    if (noise_model=='queue-jumping'){
      r = log(p/(1-p))
      r_temp = rnorm(1,r,1)
      p_temp = 1/(1+exp(-r_temp))
      loglkd_temp = loglik(trees,Y,p_temp,theta)
      log_accept_rate = loglkd_temp+dnorm(r_temp,0,R_HYPERPAMA,log=TRUE)-loglkd-
        dnorm(r,0,R_HYPERPAMA,log=TRUE)
      if(log_accept_rate > log(runif(1))){p = p_temp; loglkd = loglkd_temp}
    }
    # update on tau - WORKS! 
    # prior: t=log(tau/(1-tau)) ~ Normal(0, T_HYPERPAMA)
    t = log(tau/(1-tau))
    t_temp = rnorm(1,t,1)
    tau_temp = 1/(1+exp(-t_temp))
    log_accept_rate = dnorm(t_temp,0,T_HYPERPAMA,log=TRUE)-dnorm(t,0,T_HYPERPAMA,log=TRUE)+
      pU(U,U0,tau_temp,Sigma)-pU(U,U0,tau,Sigma)
    if(log_accept_rate > log(runif(1))){tau=tau_temp}
    
    if (itr %% n_record == 0){
      U0_PATH[[itr/n_record]] = U0
      U_PATH[[itr/n_record]] = U
      RHO_PATH[itr/n_record] = rho
      THETA_PATH[itr/n_record] = theta
      P_PATH[itr/n_record] = p
      TAU_PATH[itr/n_record] = tau
      LOGLKD_PATH[itr/n_record] = loglkd
      print(itr)
      output = list(U0=U0_PATH,U=U_PATH,rho=RHO_PATH,p=P_PATH,theta=THETA_PATH,
                    tau=TAU_PATH,loglkd=LOGLKD_PATH)
      save(output, file=paste0(chain_name,".RData"))
    }
    # if (itr %% (n_itr/10) == 0){
    # 
    # }
  }
  return(output)
}


# # update on U0
# U0_temp = U0
# c = sample(NUM_ACTORS,1); i = sample(2,1)
# # U0_temp[c,i] = rcmvnorm(n=1, mean=U0[c,], sigma=Sigma, dependent.ind=i, 
# #                         given.ind=(1:2)[-i], X.given=U0[c,-i])
# U0_temp[c,i]=rnorm(1,U0[c,i],1)
# is_vsp = po2tree(u2po(U0_temp))$is_vsp
# if (is_vsp){
#   log_accept_rate = pU0(U0_temp,Sigma)-pU0(U0,Sigma)+pU(U,U0=U0_temp,tau,Sigma)-
#     pU(U,U0=U0,tau,Sigma)
#   if (log_accept_rate > log(runif(1))) {U0=U0_temp}
# }
# # update on U
# for (j in 1:NUM_ASSESSORS){
#   U_temp=U
#   u_temp = U[[j]]
#   c = sample(NUM_ACTORS,1); i = sample(2,1)
#   # u_temp[c,i] = rcmvnorm(n=1, mean=u_temp[c,], sigma=(1-tau^2)*Sigma, 
#   #                        dependent.ind=i,given.ind=(1:2)[-i], X.given=u_temp[c,-i])
#   u_temp[c,i]=rnorm(1,U[[j]][c,i],1)
#   U_temp[[j]] = u_temp
#   po_j_temp = u2po(u_temp)
#   pos_temp = pos; trees_temp = trees
#   pos_temp[[j]] = po_j_temp
#   if (all(pos[[j]] == po_j_temp)){
#     loglkd_temp=loglkd
#     log_accept_rate = pU(U_temp,U0,tau,Sigma) - pU(U,U0,tau,Sigma)
#   } else {
#     potree = po2tree(po_j_temp)
#     if (potree$is_vsp){
#       trees_temp[[j]] = potree$tree
#       loglkd_temp = loglik(trees_temp,Y,p,theta)
#       log_accept_rate = loglkd_temp - loglkd + pU(U_temp,U0,tau,Sigma) - 
#         pU(U,U0,tau,Sigma)
#     }
#   }
#   if (log_accept_rate > log(runif(1))) {
#     U=U_temp; loglkd=loglkd_temp; trees=trees_temp; pos=pos_temp
#     }
# }

