# This script contains the main function for hierarchical vsp mcmc update.

vsp_hierarchy = function(initial_states,chain_name,n_itr=10000,n_record=NULL,
                         noise_model='queue-jumping',hierarchy=TRUE,clustering=TRUE){
  
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
  #     S:     list.
  #       The clustering to data list.
  #     Y:     list. 
  #       The nested data list. Each data entry contains attributes 'assessor' and 
  #         'order'. 
  # 
  # chain_name: string. 
  #   The name of the output file. 
  #
  # n_itr: int.
  #   The number of MCMC iterations. 
  #
  # n_record: int. 
  #   The step between each path record. We only record every n_record iterations. 
  #     If NULL, n_record = 2*NUM_ACTORS. 
  #
  # noise_model: string.
  #   The choice of noise model, either 'queue-jumping' (vsp) for 'mallows' 
  #     (general po). 
  # 
  # hierarchy: bool.
  #   Whether to fit the hierarchical VSP model or the independent VSP model.  
  # 
  # clustering: bool.
  #   Whether to perform clustering on data list.
  # 
  # Returns
  # -------
  #   List, with attributes:
  #     U0, U, rho, p, theta, tau, S: same as above. 
  #     loglkd: the trace of log-likelihood.
  
  # load initial states
  Y = initial_states$Y
  U0 = initial_states$U0
  U = initial_states$U
  rho = initial_states$rho
  p = initial_states$p
  theta = initial_states$theta
  if (hierarchy) {tau=initial_states$tau}else {tau=0}
  if (clustering) {
    S=initial_states$S
    if(length(S)!=length(U)){stop("Error type 1: the number of assessors does not match!")}
    if(length(S)!=max(sapply(Y,'[[','assessor'))){stop("Error type 2: the number of assessors does not match!")}
  } else {
    assessors=sapply(Y,'[[','assessor')
    S=lapply(1:length(unique(assessors)), function(t) which(assessors == t))
  }
  
  Sigma=matrix(rho,nrow=2,ncol=2); diag(Sigma)=1
  pos = lapply(U, u2po)
  trees = lapply(pos, function(po) po2tree(po)$tree)

  # set recording steps
  if(is.null(n_record)) {n_record = 2*NUM_ACTORS}
  
  # define likelihood function
  if (noise_model=='queue-jumping'){loglik=loglikQJ;loglik_y=loglikQJ_y}
  if (noise_model=='mallows')      {loglik=loglikM}
  if (noise_model=='test')         {loglik=function(trees,Y,p,theta){return(0)};loglik_y=function(tree,y,p,theta){return(0)}}
  
  loglkd = loglik(trees,Y,p,theta)
  
  # initiate storage lists/vectors
  U0_PATH = vector(mode= "list", length = n_itr/n_record)
  U_PATH = vector(mode="list", length = n_itr/n_record)
  RHO_PATH = vector(mode="numeric", length = n_itr/n_record)
  TAU_PATH = vector(mode= "numeric", length = n_itr/n_record)
  THETA_PATH = vector(mode= "numeric", length = (n_itr/n_record)) # for MALLOWS
  P_PATH = vector(mode= "numeric", length = (n_itr/n_record)) # for QUEUE-JUMPING
  S_PATH = vector(mode= "list", length = n_itr/n_record)
  LOGLKD_PATH = vector(mode="numeric", length = n_itr/n_record)
  
  for(itr in 1:n_itr){
    
    if (isTRUE(clustering)){
      # update on S & U
      e=sample(NUM_LISTS,1)
      clusterIndex=which(sapply(S, FUN=function(X) e %in% X))
      listIndex=which(S[[clusterIndex]]==e)
      
      if(length(S[[clusterIndex]])==1){
        # calculate weights
        weights = sapply(trees,loglik_y,y=Y[[e]],p=p)
        weights_constant = lengths(S)-PDP_ALPHA
        weights_constant[clusterIndex]=PDP_THETA + (length(S)-1)*PDP_ALPHA
        weights = exp(weights)*weights_constant
        # sample new cluster
        ind = sample(1:length(weights),1,prob = weights)
        if(ind!=clusterIndex){
          # dimension down
          S[[ind]] = c(S[[ind]],e)
          S = S[-clusterIndex];U = U[-clusterIndex];trees = trees[-clusterIndex]
          Y = YbyS(Y,S)
          loglkd = loglik(trees,Y,p,theta=theta)
        }
      } else {
        # generate u for new cluster (conditional on being VSP)
        u = tau*U0+rmvnorm(NUM_ACTORS,sigma=(1-tau^2)*Sigma)
        potree=po2tree(u2po(u))
        while(!potree$is_vsp){
          u = tau*U0+rmvnorm(NUM_ACTORS,sigma=(1-tau^2)*Sigma)
          potree=po2tree(u2po(u))
        }
        # calculate weights
        S_reduced = S
        S_reduced[[clusterIndex]] = S_reduced[[clusterIndex]][-listIndex]
        weights = sapply(trees,loglik_y,y=Y[[e]],p=p)
        weights = exp(weights) * (lengths(S_reduced)-PDP_ALPHA)
        weight_new = (PDP_THETA+length(S)*PDP_ALPHA)*exp(loglik_y(potree$tree,
                                                                  Y[[e]],p=p))
        weights = c(weights, weight_new)
        # sample new cluster
        ind = sample(1:length(weights),1,prob = weights)
        if(ind==length(weights)){
          # dimension up
          S = append(S_reduced,e)
          U = append(U, list(u))
          trees = append(trees, list(potree$tree))
          Y = YbyS(Y,S)
          loglkd = loglik(trees,Y,p,theta=theta)
        } else {
          if(ind!=clusterIndex){
            # dimension remain unchanged
            S_reduced[[ind]] = c(S_reduced[[ind]],e)
            S = S_reduced
            Y = YbyS(Y,S)
            loglkd = loglik(trees,Y,p,theta=theta)
          }
        }
      }
    }
    
    # update on U0
    U0_temp = U0
    c = sample(NUM_ACTORS,1); i = sample(2,1)
    U0_temp[c,i]=rnorm(1,U0[c,i],sd=1)
    log_accept_rate = pU0(U0_temp,Sigma) + pU(U,U0_temp,tau,Sigma)-
      pU0(U0,Sigma) - pU(U,U0,tau,Sigma)
    if (log_accept_rate > log(runif(1))) {U0 = U0_temp}
    
    # update on U
    for (assrIndex in 1:length(U)){
      U_temp = U
      trees_temp = trees
      u_temp = u = U_temp[[assrIndex]]
      c = sample(NUM_ACTORS,1); i = sample(2,1)
      u_temp[c,i] = rnorm(1,u[c,i],sd=1)
      po_temp = u2po(u_temp)
      trees_temp = trees
      loglkd_temp = loglkd
      if(!all(u2po(u_temp) == u2po(u))){
        potree_tmp = po2tree(po_temp)
        if (potree_tmp$is_vsp){
          trees_temp[[assrIndex]]=potree_tmp$tree
          U_temp[[assrIndex]] = u_temp
          loglkd_temp = loglik(trees_temp,Y,p,theta)
          log_accept_rate = loglkd_temp-loglkd+
            pU(U_temp,U0,tau,Sigma)-pU(U,U0,tau,Sigma)
        }
      } else {
        U_temp[[assrIndex]] = u_temp
        log_accept_rate = pU(U_temp,U0,tau,Sigma)-pU(U,U0,tau,Sigma)
      }
      if (log_accept_rate > log(runif(1))) {
        U=U_temp; trees=trees_temp; loglkd=loglkd_temp
      }
    }

    # update on U0 and U
    U0_temp = U0
    c = sample(NUM_ACTORS,1)
    U0_temp_c = rmvnorm(n=1,sigma=Sigma)
    U0_temp[c,] = U0_temp_c
    U_temp = U
    trees_temp = trees
    allvsp = TRUE
    for (assrIndex in 1:length(U)){
      u = U_temp[[assrIndex]]
      u[c,] = rmvnorm(n=1,mean=tau*U0_temp_c,sigma=(1-tau^2)*Sigma)
      potree_tmp = po2tree(u2po(u))
      if (potree_tmp$is_vsp){
        trees_temp[[assrIndex]]=potree_tmp$tree
        U_temp[[assrIndex]] = u
      } else {allvsp=FALSE; break}
    }
    if(allvsp){
      loglkd_temp = loglik(trees_temp,Y,p,theta)
      log_accept_rate = loglkd_temp-loglkd
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
    if(hierarchy){
      # prior: tau~U[0,1]
      tau_temp = runif(1)
      log_accept_rate = pU(U,U0,tau_temp,Sigma)-pU(U,U0,tau,Sigma)
      if(log_accept_rate > log(runif(1))){tau=tau_temp}
    }
    if (itr %% n_record == 0){
      U0_PATH[[itr/n_record]] = U0
      U_PATH[[itr/n_record]] = U
      S_PATH[[itr/n_record]] = S
      RHO_PATH[itr/n_record] = rho
      THETA_PATH[itr/n_record] = theta
      P_PATH[itr/n_record] = p
      TAU_PATH[itr/n_record] = tau
      LOGLKD_PATH[itr/n_record] = loglkd
      print(itr)
      output = list(U0=U0_PATH,U=U_PATH,rho=RHO_PATH,p=P_PATH,theta=THETA_PATH,
                    tau=TAU_PATH,S=S_PATH,loglkd=LOGLKD_PATH)
      save(output, file=paste0(chain_name,".RData"))
    }
  }
  return(output)
}
