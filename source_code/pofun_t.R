# POwTIE functions
## This file contains all relevent functions for 'POPLModel-Tie.Rmd'

############################## Basic Set-up Functions ##############################
CRP = function(n, alpha,list=FALSE){
  nDnr <- n  # number of diners; must be at least 2
  vDnrTbl <- rep(0, nDnr)  # table label for each diner
  
  vDnrTbl[1] <- 1
  for (dnr in 2:length(vDnrTbl)) {
    
    # compute occupation probabilities for current diner
    vOcc <- table(vDnrTbl[1:(dnr-1)])
    vProb <- c(vOcc, alpha) / (dnr - 1 + alpha)
    
    # add table label to diner
    nTbl <- as.numeric(names(vOcc)[length(vOcc)])  # avoid overhead of finding max of possibly large vector
    vDnrTbl[dnr] <- sample.int(nTbl+1, size=1, prob=vProb)
  }
  if (list == FALSE){return(vDnrTbl)} else {
    nTbl <- max(c(nTbl, vDnrTbl[dnr]))
    listTableOccupants <- lapply(1:nTbl, function(t) which(vDnrTbl == t))
    return(listTableOccupants)
  }
}
pCRP = function(list,alpha){
  A = length(list)
  term1=gamma(alpha)/(gamma(alpha/A)^A)
  term2=prod(gamma(lengths(list)+alpha/A)) / gamma(length(unlist(list))+alpha/A)
  return(term1*term2)
}
RandPar_tie = function(n,K,rho,alpha){
  # This is the only function that treat the tie scenario as 1-1, else all as incompariable
  ## IMPORTANT ##
  # treating the tie scenario as 1-1; 
  # in the rest of the code - treating tie to be equivalent as incomparible while using S to indicate tie 
  ls=CRP(n,alpha,list=TRUE)
  
  Sigma=matrix(rho, nrow=K,ncol=K)
  diag(Sigma)=1
  Z=rmvnorm(n,sigma=Sigma)
  for (i in 1:length(ls)){
    ni = length(ls[[i]])
    Z[ls[[i]],]=matrix(rep(Z[ls[[i]][1],],each=ni),nrow=ni)
  }
  if(length(ls)==1){m=matrix(1,nrow=n,ncol=n);diag(m)=0} else{
    m=matrix(NA,nrow=n,ncol=n)
    m[unlist(lapply(ls, `[[`, 1)),unlist(lapply(ls, `[[`, 1))]=order2partial(v = latent2order(unique(Z)), n = length(ls))
    for(j in 1:length(ls)){
      ind=ls[[j]]
      if(length(ind)>1){
        for(k in ind[-1]){
          m[k,-k]=m[ind[1],-k]
          m[-k,k]=m[-k,ind[1]]
        }
        m[expand.grid(ind,ind)[,1],expand.grid(ind,ind)[,2]]=1
      }
    }
    diag(m)=0
  }
  return(transitive.reduction(m))
}
# TO DO: refine showDAG function to include undirected edge - tie (add S in)

# Transformation on Z
Z_tr_tie = function(Z,S){
  # Transfer a Z-matrix into a partial order considering the ties
  PO = transitive.closure(order2partial(v = latent2order(Z), n = nrow(Z)), mat = TRUE, loop=FALSE)
  for(j in 1:length(S)){
    ind=S[[j]]
    PO[expand.grid(ind,ind)[,1],expand.grid(ind,ind)[,2]]=0}
  return(transitive.reduction(PO))
}
ZbyS = function(Z,S){
  # Refine a Z-matrix based on S
  for (i in 1:length(S)){
    ni = length(S[[i]])
    Z[S[[i]],]=matrix(rep(Z[S[[i]][1],],each=ni),nrow=ni)
  }
  return(Z)
}
S_sort = function(S){
  S=lapply(S,sort)
  S=S[order(sapply(S, function(x) x[1], simplify=TRUE), decreasing=FALSE)]
  return(S)
}
############################## MCMC #################################################
# Loglikelihood + logZ
loglikP=function(Z,S,Y,p=NULL,AlphaPL=NULL){
  tr=Z_tr_tie(Z,S)
  tc=transitive.closure(tr,mat=TRUE, loops = FALSE)
  
  if (noconflict4(Y,PO=tr)){
    if (length(unique(lengths(Y)))==1){
      N=length(Y)
      loglk=(-N*log(nle(tr)))
    } else {
      loglkp = function(y,Z,S){
        -log(nle(transitive.reduction(tc[sort(y),sort(y)])))
      }
      loglk=sum(unlist(lapply(Y,loglkp,Z=Z,S=S)))
    }
  } else {loglk = -Inf}
  return(loglk)
}
loglikQ=function(Z,S,Y,p,AlphaPL=NULL){ 
  # queue-jumping (down)
  tr=Z_tr_tie(Z,S)
  tc=transitive.closure(tr,mat=TRUE,loops=FALSE)
  
  loglik_j=function(tc,y,p,j){
    n = length(y)
    tc_temp = tc[sort(y[j:n]),sort(y[j:n])]
    tr_temp = transitive.reduction(tc_temp)
    ind = which(colnames(tr_temp) == y[j])
    if(apply(tr_temp,2,sum)[ind]==0){
      tc_sub = tc_temp[-ind,-ind]
      return(log(p/(n-j+1)+(1-p)*nle(tc_sub)/nle(tr_temp))) # nle works the same for both tc and tr
    } else {return(log(p/(n-j+1)))}
  }  
  
  loglik_y=function(y,tc,p){
    n = length(y)
    sum(unlist(lapply(1:(n-1),loglik_j,tc=tc,y=y,p=p)))
  }
  return(sum(unlist(lapply(Y,loglik_y,tc=tc,p=p))))
}

# sum(loglkd(r=p,mc=transitive.closure(Z_tr_tie(Z,S),mat=TRUE,loops=FALSE),la=Y, model='lkddown')) cross-check

# loglik = function(Z,S,Y){return(0)}
logZ=function(Z,rho,k){
  # calculates the probability of a unique Z
  Z = unique(Z)
  Sigma=matrix(rho,nrow=k,ncol=k); diag(Sigma)=1
  return(sum(apply(Z,1,dmvnorm,log=TRUE,sigma=Sigma)))
} 
# MCMC algorithm

tie_update = function(k, Z, S, rho, p=0, Y, alpha, n.itr=100000, q=0.3, w=0.50, model,chain.index){
  # there are three choices of 'model' - 'p' (error-free), 'q' (queue-jumping) and 'pl' (plackett-luce)
  if (model == "p"){loglik = loglikP}
  if (model == "q"){loglik = loglikQ}
  if (model == "pl"){loglik = loglikPL} 
  
  k_Update=function(k,Z,S,rho,p,Y,q,loglk,PO){
    
    k_temp=sample(c(k+1,k-1),1)
    
    condZ = function(Zi,rho,k){
      Sigma = matrix(rho, k, k); diag(Sigma) = 1
      return(c(Zi,rcmvnorm(1,mean=rep(0,k),sigma=Sigma,dependent.ind=k,given.ind=1:(k-1),X.given=Zi)))
    }
    
    if(k_temp>0){
      if(k_temp==k+1){
        Z_temp=t(apply(Z,1,condZ,rho=rho,k=k_temp))
        Z_temp=ZbyS(Z=Z_temp,S=S)
      } else {
        Z_temp=Z[,1:k_temp,drop=FALSE]
      }
      PO_temp=Z_tr_tie(Z=Z_temp,S=S)
      if (all(PO_temp==PO)){
        loglk_temp=loglk
        log_eta2=dgeom(x=k_temp,prob=q,log=TRUE)-dgeom(x=k,prob=q,log=TRUE)
      } else {
        loglk_temp=loglik(Z=Z_temp,S=S,Y=Y,p=p)
        log_eta2=loglk_temp-loglk+dgeom(x=k_temp,prob=q,log=TRUE)-dgeom(x=k,prob=q,log=TRUE)
      }
      if (log_eta2>log(runif(1))){Z=Z_temp;k=k_temp;loglk=loglk_temp;PO=PO_temp}
    }
    return(list(Z,k,loglk,PO))
  }
  
  S_update_MH = function(S,Z,rho,p,k,Y,alpha,PO,loglk){
    
    add_ele = function(S,e,w,Z,rho,k){
      if (w<=length(S)){
        S[[w]]=c(S[[w]],e)
        Z[e,]=Z[S[[w]][1],]
      } else {
        S=append(S,e)
        Sigma=matrix(rho,nrow=k,ncol=k); diag(Sigma)=1
        Z[e,]=rmvnorm(1,sigma=Sigma)
      }
      return(list(S,Z))
    }
    
    n=nrow(Z)
    e=sample(n,1)
    c=which(sapply(S, FUN=function(X) e %in% X))
    j=which(S[[c]]==e)
    
    S_reduced=S
    S_reduced[[c]]=S_reduced[[c]][-j]
    S_reduced=Filter(length,S_reduced)
    
    R_temp = lapply(1:(length(S_reduced)+1),add_ele,S=S_reduced,e=e,Z=Z,rho=rho,k=k)
    S_temp = lapply(R_temp, `[[`,1)
    S_temp = lapply(S_temp, function(S){S=lapply(S,sort); S=S[order(sapply(S,function(x) x[1], simplify=TRUE), decreasing=FALSE)]})
    S_temp[[length(S_temp)+1]] = S
    Z_temp = lapply(R_temp, `[[`,2)
    Z_temp[[length(Z_temp)+1]] = Z
    PO_temp = lapply(R_temp,function(X){Z_tr_tie(X[[2]],X[[1]])})
    PO_temp[[length(PO_temp)+1]] = PO
    
    set_PO=which(sapply(PO_temp, function(X,PO){all(X==PO)},PO=PO))
    set_diff=setdiff(1:(length(S_reduced)+2),set_PO)
    
    loglk_temp = rep(NA,length(S_reduced)+2)
    loglk_temp[set_PO]=loglk
    loglk_temp[set_diff] = unlist(lapply(R_temp[set_diff], function(X, Y, p){loglik(X[[2]],X[[1]],Y,p)},Y=Y,p=p))
    
    prob = c(lengths(S_reduced),rep(alpha,2))*exp(loglk_temp)
    ind=sample((length(S_reduced)+2),1,prob=prob)
    
    S = S_temp[[ind]]
    Z = Z_temp[[ind]]
    loglk = loglk_temp[ind]
    PO = PO_temp[[ind]]

    return(list(S,Z,loglk,PO))
  }
  
  ptm = proc.time()
  n = nrow(Z)
  N = length(Y)
  n.record = max(n^2, n*N)
  
  Z_res = vector(mode= "list", length = n.itr/n.record)
  rho_res = vector(mode="numeric", length = n.itr/n.record)
  p_res = vector(mode="numeric", length = n.itr/n.record)
  k_res = vector(mode="numeric", length = n.itr/n.record)
  S_res = vector(mode= "list", length = n.itr/n.record)
  PO_res = vector(mode= "list", length = (n.itr/n.record-1))
  loglk_res = vector(mode="numeric", length = (n.itr/n.record-1))
  
  Z_res[[1]]=Z
  rho_res[1]=rho
  p_res[1]=p
  k_res[1]=k
  S_res[[1]]=S
  
  PO=Z_tr_tie(Z,S)
  loglk=loglik(Z,S,Y,p) # the initial state must fit data
  
  for (i in 2:n.itr){
    
    # update on k & Z - works!
    res1=k_Update(k,Z,S,rho,p,Y,q=q,loglk,PO)
    Z=res1[[1]]
    k=res1[[2]]
    loglk=res1[[3]]
    PO=res1[[4]]
    
    ## update on S & Z - prior ok!
    # res2=S_update_MH(S,Z,rho,p,k,Y,alpha,PO,loglk)
    # S=res2[[1]]
    # Z=res2[[2]]
    # loglk=res2[[3]]
    # PO=res2[[4]]
    
    # update on Z - works!
    iZ=sample(nrow(unique(Z)),1)
    jZ=sample(ncol(Z),1)

    Z_temp=Z
    Sigma = matrix(rho, k, k);diag(Sigma) = 1
    Z_temp[S[[iZ]],jZ]=rcmvnorm(n = 1, mean = unique(Z)[iZ,], sigma = Sigma, dependent.ind = jZ, given.ind = (1:k)[-jZ], X.given = unique(Z)[iZ,-jZ])
    PO_temp=Z_tr_tie(Z_temp,S)

    if(all(PO_temp==PO)){
      loglk_temp=loglk
      log_eta1=logZ(Z=Z_temp,rho=rho,k=k)-logZ(Z=Z,rho=rho,k=k)
    } else {
      loglk_temp=loglik(Z=Z_temp,S=S,Y=Y,p=p)
      log_eta1 = loglk_temp+logZ(Z=Z_temp,rho=rho,k=k)-loglk-logZ(Z=Z,rho=rho,k=k)
    }
    if (log_eta1 > log(runif(1))) {Z = Z_temp; PO = PO_temp; loglk=loglk_temp}

    # update on rho - works!
    delta = runif(1, w, 1/w)
    rho_temp = 1 - delta*(1-rho)
    log_eta3 = logZ(Z=Z,rho=rho_temp,k=k)+dbeta(rho_temp,1,1/5,log=TRUE)-logZ(Z=Z,rho=rho,k=k)-dbeta(rho,1,1/5,log=TRUE)-log(delta)
    if(log_eta3 > log(runif(1)) & rho_temp<(1-10^(-7))){rho = rho_temp}

    # update on p
    if (model == "q"){
      deltap = runif(1, w, 1/w) # might change w potentially
      p_temp = 1 - deltap*(1-p)
      if (p_temp>0){
        loglk_temp = loglik(Z=Z,S,Y,p_temp)
        log_eta4 = loglk_temp +dbeta(p_temp,2,2,log=TRUE)-loglk-dbeta(p,2,2,log=TRUE)-log(deltap)
        if(log_eta4 > log(runif(1)) & p_temp<(1-10^(-7))){p = p_temp; loglk=loglk_temp}
      }
    }
    
    if (i %% n.record == 0){
      Z_res[[i/n.record+1]] = Z
      rho_res[i/n.record+1] = rho
      p_res[i/n.record+1] = p
      k_res[i/n.record+1] = k
      S_res[[i/n.record+1]] = S
      loglk_res[i/n.record] = loglk
      PO_res[[i/n.record]] = PO
    }
    if (i %% (n.itr/10) == 0){
    	print(i)
    	Result = list(Z=Z_res,k=k_res,S=S_res,rho=rho_res,p=p_res,loglk=loglk_res,PO=PO_res)
    	save(Result, file=paste0('Bishop',chain.index,'.RData'))
    	}
  }
  print(proc.time() - ptm)
  # return(list(Z=Z_res,k=k_res,S=S_res,rho=rho_res,p=p_res,loglk=loglk_res,PO=PO_res))
}

############################## Result Analysis ###############################################

plota=function(dag.con, display.threshold=0.5,label.threshold=0.5,fontsize=8,main, tc = TRUE){
  # the following part is adapted from the dag.conc2 function
  if(tc){dag.con2<-(dag.con>display.threshold)} else {dag.con2<-transitive.reduction(dag.con>display.threshold)}
  dag.con<-round(dag.con,2)*dag.con2
  dag.con<-new("graphAM", adjMat=dag.con, edgemode="directed", values=list(weight=1))
  dag.con<-as(dag.con,"graphNEL")
  ewv<-unlist(edgeWeights(dag.con))
  eAttrs <- list()
  n.edges=length(edgeNames(dag.con))
  e.color <- rep("red",n.edges); names(e.color)<-edgeNames(dag.con)
  e.color[ewv<0.9]<-"black"
  eAttrs$color<-e.color
  nAttrs<-list()
  n.nodes=length(nodes(dag.con))
  nAttrs$height <- nAttrs$width <- rep("0", n.nodes)
  nAttrs$style <- rep("invis",n.nodes)
  nAttrs$fontsize <- rep(fontsize, n.nodes)
  nAttrs <- lapply(nAttrs, function(x) {names(x) <- nodes(dag.con); x })
  attrs <- list(node = list(shape = "plaintext", fixedsize = "false"), edge=list(headclip="false",tailclip="true",style="dashed",arrowsize=2.0), graph=list(splines="false"))
  h = list(dag.con=dag.con,nodeAttrs=nAttrs,edgeAttrs=eAttrs,attrs=attrs)
  numEdges=length(unlist(edgeL(h$dag.con)))
  if (numEdges==0) {
    a=plot(h$dag.con,attrs=h$attrs,"dot")
  } else {
    a=plot(h$dag.con,nodeAttrs=h$nodeAttrs,edgeAttrs=h$edgeAttrs,attrs=h$attrs,"dot")
  }
  return(a)
}

str2po=function(x,n){
  matrix(as.numeric(strsplit(x," ")[[1]]),nrow=n)
}


