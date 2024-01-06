# This script contains helper functions to calculate the Mallow's likelihood 
# under the partial-order Mallow's model.

# required packages
library(rankdist)
library(combinat)

kendall <- function(l,y){
  # calculates the Kendall's tau distance
  # Parameters
  # ----------
  # l, y: vector.
  #   the ranking vectors. 
  # Return
  # ------
  #   float.
  #   The Kendall's tau distance between l and y. 
  n = length(l)
  d = 0
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      d = d + sum(which(y==l[i])>which(y==l[j]))
    }
  }
  return(d)
}

d.kendall <- function(theta,l,y,log=TRUE){
  # calculates the likelihood for the Mallow's phi model (Kendall's tau distance)
  # Parameters
  # ----------
  # theta: float in (0, inf).
  #   the dispersion parameter.
  # l, y: vector.
  #   the ranking vectors. 
  # log: Bool.
  #   whether to return log-likelihood or likelihood
  # Return
  # ------
  #   float.
  #   The (log-)likelihood for the Mallow's phi model. 
  n = length(y)
  if(length(l)!=n) stop('the rankings are not of the same lengths')
  p=exp(-theta*kendall(l,y))/prod(cumsum(exp(-c(0:(n-1))*theta)))
  if (isTRUE(log)) {return(log(p))} else {return(p)}
}

.mallows.i <- function(i,theta,y,tr,measure='kendall'){
  # calulates p(i|theta,y) for repeated selection
  # Parameters
  # ----------
  # i: int. 
  #   the actor of interest. 
  # theta: float in (0, inf).
  #   the dispersion parameter
  # y: vector. 
  #   the ranking vector
  # tr: matrix. 
  #   the transitive reduction of a partial order
  # measure: string.
  #   the distance measure in the Mallow's model; currently only support 'kendall' (Kendall tau distance). 
  # Return
  # ------
  # float. 
  #   p(i|theta,y).
  n = nrow(tr)
  if(measure=='kendall'){d = sum(which(y==i)>sapply(setdiff(sort(y),i),function(a){which(y==a)}))}
  return(exp(-theta*d)/sum(exp(-c(0:(length(y)-1))*theta)))
}

mallows <- function(tr, theta, y, p=1,measure='kendall'){
  # Mallow's log-likelihood on a partial order.
  # parameters
  # ----------
  # tr: matrix.
  #   the transitive reduction of a partial order
  # theta: float in (0, inf).
  #   the dispersion parameter in the mallows model
  # y: vector.
  #   the data list. 
  # p: float.
  #   the likelihood on the top actors. Default p=1.
  # return
  # ------
  #   float. 
  #   the Mallow's likelihood on a partial order. 
  if (length(tr)==1) {return(list(p=p,count=1,lik=0))}
  if (nrow(tr)!=length(y)) {tr = transitive.reduction(transitive.closure(tr,mat=TRUE,loops=FALSE)[sort(y),sort(y)])}
  n = nrow(tr)
  cs<-apply(tr,2,sum)
  if (sum(cs)==0) {return(list(p=p,count=factorial(n),lik=-log(factorial(n))))}
  csi<-(cs==0)
  bs<-apply(tr,1,sum)
  bsi<-(bs==0)

  tops<-which(csi)
  bots<-which(bsi)
  
  g <- 0 
  count <- 0
  for (i in tops) {
    e = sort(y)[i]
    pi = p*.mallows.i(i = e,theta,y,tr,measure=measure)
    trr<-tr[-i,-i]
    r = mallows(trr, theta, y[-which(y==e)],p=pi,measure=measure)
    g = g + r$p
    count = count + r$count
  }
  return(list(p=g,count=count,lik=log(g)-log(count)))
}

loglikM <- function(tr,theta,Y){
  # Mallow's log-likelihood on data list Y given a partial order tr.
  # Parameters
  # ----------
  # tr: matrix.
  #   the transitive reduction of a partial order
  # theta: float in (0, inf).
  #   the dispersion parameter in the mallows model
  # Y: list of vectors.
  #   the data list.  
  # Return
  # ------
  #   float. 
  #   the Mallow's log-likelihood on data list Y given a partial order tr
  sum(sapply(Y,function(y,tr,theta) {mallows(tr=tr,y=y,theta=theta)$lik}
             ,theta=theta,tr=tr))
}

### test ###

# source("~/Desktop/POwTie/code/source-code/pofun.R")
# source("~/Desktop/POwTie/code/source-code/pofun_t.R")
# 
# n = 5
# theta = 8
# 
# test=rep(NA,1000)
# 
# for (i in 1:1000){
#   tr = transitive.reduction(randDGDAG(n=n))
#   rownames(tr) = colnames(tr) = 1:n
#   tc = transitive.closure(tr,mat=TRUE,loops=FALSE)
#   y = sample(n,n)
#   l = sample(n,n)
#   
#   a= log(sum(sapply(allLE(tc),d.kendall,theta=theta,y=y,log=FALSE))/nle(tr));a
#   
#   r = mallows(tr,theta,y)
#   b=r$lik;b
#   
#   test[i]=all.equal(a,b)
# }
# 
# unique(test)

