# nem functions

library(MASS)
library(igraph)

transitive.reduction = function (g){
  if (!(class(g)[1] %in% c("matrix", "graphNEL"))) 
    stop("Input must be an adjacency matrix or graphNEL object")
  if (class(g)[1] == "graphNEL") {
    g = as(g, "matrix")
  }
  g = transitive.closure(g, mat = TRUE)
  g = g - diag(diag(g))
  type = (g > 1) * 1 - (g < 0) * 1
  for (y in 1:nrow(g)) {
    for (x in 1:nrow(g)) {
      if (g[x, y] != 0) {
        for (j in 1:nrow(g)) {
          if ((g[y, j] != 0) & sign(type[x, j]) * sign(type[x, 
                                                            y]) * sign(type[y, j]) != -1) {
            g[x, j] = 0
          }
        }
      }
    }
  }
  g
}

transitive.closure = function (g, mat = FALSE, loops = TRUE) {    
  if (!(class(g)[1] %in% c("graphNEL", "matrix")))         
    stop("Input must be either graphNEL object or adjacency matrix")    
  g <- as(g, "matrix")    
  n <- ncol(g)    
  matExpIterativ <- function(x, pow, y = x, z = x, i = 1) {        
    while (i < pow) {            
      z <- z %*% x            
      y <- y + z            
      i <- i + 1}        
    return(y)}    
  h <- matExpIterativ(g, n)    
  h <- (h > 0) * 1    
  dimnames(h) <- dimnames(g)    
  if (!loops)         
    diag(h) <- rep(0, n)    
  else diag(h) <- rep(1, n)    
  if (!mat)         
    h <- as(h, "graphNEL")    
  return(h)
  
  if (!(class(g)[1] %in% c("matrix", "graphNEL")))         
    stop("Input must be an adjacency matrix or graphNEL object")    
  if (class(g)[1] == "graphNEL") {g = as(g, "matrix")}    
  g = transitive.closure(g, mat = TRUE)    
  g = g - diag(diag(g))    
  type = (g > 1) * 1 - (g < 0) * 1    
  for (y in 1:nrow(g)) {        
    for (x in 1:nrow(g)) {            
      if (g[x, y] != 0) {                
        for (j in 1:nrow(g)) {                  
          if ((g[y, j] != 0) & sign(type[x, j]) * sign(type[x,y]) * sign(type[y, j]) != -1) {g[x, j] = 0}}}}}    
  g}

# DAG

showDAG<-function(m=NULL,...) {
  g<-graph_from_adjacency_matrix(m,mode ="directed")#as(m,"graphNEL")
  #h<-graph_from_adjacency_matrix(transitive.closure(m,mat=TRUE,loops=FALSE),mode ="directed")
  h<-g #seems like a bug in layout function puts vertices on top of one another
  plot(g,layout=layout_with_sugiyama(h)$layout,...);   	 
}

showDAGTRTC<-function(m=NULL,mr=NULL,mc=NULL,mt=NULL,t=NULL) {
  
  #draw DAG's (adjacency matrix format)
  
  #example
  #p<-0.5
  #n<-15
  #for (k in 1:10) {
  #	m<-randDGDAG(n,p,DAG=TRUE)
  #	mc <- transitive.closure(m,mat=TRUE,loops=FALSE)
  #	mr <- transitive.reduction(mc) 
  #	showDAGTRTC(m=m,mr=mr,mc=mc,t=k)
  #}
  
  pics<-0	
  if (!is.null(m)) {
    pics<-pics+1
    g<-graph_from_adjacency_matrix(m,mode ="directed")#as(m,"graphNEL")
  }
  if (!is.null(mc)) {
    pics<-pics+1
    gc<-graph_from_adjacency_matrix(mc,mode ="directed")#as(mc,"graphNEL")
  }
  if (!is.null(mr)) {
    pics<-pics+1
    gr<-graph_from_adjacency_matrix(mr,mode ="directed")#as(mr,"graphNEL")
  }
  if (!is.null(mt)) {
    pics<-pics+1
    gt<-graph_from_adjacency_matrix(mt,mode ="directed")#as(mt,"graphNEL")
  }
  if (pics==0) {
    warning('3 Null graphs passed to showDAGTRTC')
    return()
  } else {
    par(mfrow=c(1,pics)); 
    if (!is.null(m)) {
      plot(g,main=c('random DAG #',as.character(t))); 
    }
    if (!is.null(mr)) {
      plot(gr,main='reduction')
    }
    if (!is.null(mc)) {
      plot(gc,main='closure')
    }
    if (!is.null(mt)) {
      plot(gt,main='truth')
    }
    
  }
}

showDAGw<-function(m=NULL,...) {
  g<-graph_from_adjacency_matrix(m,mode ="directed", weighted = TRUE)
  E(g)$width <- E(g)$weight
  h<-g
  plot(g,layout=layout_with_sugiyama(h)$layout,...);   	 
}

randDGDAG<-function(n=10,p=0.5,DAG=TRUE) {
  
  #generate a random DAG (Brightwell dbn)
  #output is adjacency matrix
  
  #p<-0.5
  #n<-15
  #m<-randDGDAG(n,p,DAG=TRUE)
  #mc <- transitive.closure(m,mat=TRUE,loops=FALSE)
  #mr <- transitive.reduction(mc) 
  #showDAGTRTC(m=m,mr=mr,mc=mc)
  
  ft<-matrix(0+(runif(n*n)<p),n,n)
  if (DAG) {ft[lower.tri(ft,diag=TRUE)]<-0} else {diag(ft)<-0}
  ft
}

dagwidth<-function(tc) {
  #compute dag width from closure; uses library(igraph)
  #tc<-transitive.closure(tc,mat=TRUE,loops=FALSE)
  op=1-(tc+t(tc))
  opNEL<-as(op,'graphNEL') #XXX TODO this may not be right now
  opIG<-igraph.from.graphNEL(opNEL)
  clique.number(opIG)
}

dagdepth<-function(tr) {
  #compute dag depth from transitive reduction
  if (length(tr)==1) {return(1)}
  n<-dim(tr)[1]
  cs<-apply(tr,2,sum)
  if (sum(cs)==0) {return(1)}
  csi<-(cs==0)
  bs<-apply(tr,1,sum)
  bsi<-(bs==0)
  free<-which(bsi&csi)
  k<-length(free)
  if (k==n) {return(1)}
  if (k>0) { 
    tr<-tr[-free,-free]
    cs<-apply(tr,2,sum)
    csi<-(cs==0)
    bs<-apply(tr,1,sum)
    bsi<-(bs==0)
  }
  tops<-which(csi)
  bots<-which(bsi)
  if (length(bots)>length(tops)) {tops<-bots}
  return(1+dagdepth(tr[-tops,-tops]))
}

unifLE<-function(tc,le=NULL) {
  #given PO=tc sample one LE uniformly at random
  tc[tc==1 & t(tc)==1]=0
  if ((dim(tc)[1]==1) && (dim(tc)[2]==1)) {
    le<-c(le,as.numeric(rownames(tc)[1]))    
    return(le)
  } else {
    tops<-which(apply(tc,2,sum)==0)
    n.tops=length(tops)
    if (n.tops==1) {
      v=tops
    } else {
      weights=rep(NA,n.tops)
      for (k in 1:n.tops) {
        weights[k]=nle(tc[-tops[k],-tops[k]]) #do I need the transitive closure?
      }
      v<-sample(x=tops,size=1,prob=weights)
    }
    le<-c(le,as.numeric(rownames(tc)[v]))
    tc<-tc[-v,-v,drop=FALSE]
    return(unifLE(tc=tc,le=le))
  }
}

unifLE.tree <- function(tree,r=NULL){
  #given a vsp-tree, sample one LE uniformly at random
  if(is.null(r)){r=which(sapply(tree,function(x) is.root(x)))} #index
  if(tree[[r]]$type=='S'){
    child = tree[[r]]$child
    v = child[which(sapply(tree[child],'[[','order')=='+')] #index
    le=c(unifLE.tree(tree,r=v),unifLE.tree(tree,r=setdiff(child,v)))
  }
  if(tree[[r]]$type=='P'){
    child = tree[[r]]$child
    le1=unifLE.tree(tree,r=child[1])
    le2=unifLE.tree(tree,r=child[2])
    len = length(le1)+length(le2)
    le = rep(NA,len)
    le1s = sort(sample(len,length(le1)),decreasing=FALSE)
    le[le1s] = le1
    le[sort(setdiff(1:len,le1s),decreasing=FALSE)]=le2
  }
  if(tree[[r]]$type==''){le=tree[[r]]$actor}
  return(le)
}

newLE<-function(tc,le=NULL){
  #given PO=tc, list all LEs
  if ((dim(tc)[1]==1) && (dim(tc)[2]==1)) {
    le<-c(le,as.numeric(rownames(tc)[1]))    
    return(list(le))
  } else {
    tops<-as.numeric(rownames(tc))[which(apply(tc,2,sum)==0)]
    return(lapply(tops,function(tops,le){append(le,tops)},le=le))
  }
}

allLE<-function(tc) {
  #given PO=tc, list all LEs
  
  newLE<-function(tc,le=NULL){
    #given PO=tc, list all LEs
    if ((dim(tc)[1]==1) && (dim(tc)[2]==1)) {
      le<-c(le,as.numeric(rownames(tc)[1]))    
      return(list(le))
    } else {
      tops<-as.numeric(rownames(tc))[which(apply(tc,2,sum)==0)]
      return(lapply(tops,function(tops,le){append(le,tops)},le=le))
    }
  }
  
  tc[tc==1 & t(tc)==1]=0
  n = nrow(tc)
  les = newLE(tc)
  
  while(length(les[[1]])!=n){
    les = lapply(les, function(le,tc){
      tc = tc[-le,-le,drop=FALSE]
      return(newLE(tc,le))
    },tc=tc)
    les = lapply(rapply(les, enquote, how="unlist"), eval)
  }
  return(les)
}

# nle<-function(tr) {return(as.numeric(lecount(tr)))}

nle<-function(tr) {

  #count the number of linear extensions of the partial order
  #with transitive reduction tr (adjacency matrix)

  #example 1
  #system.time(a<-nle(mr)); a
  #example 2
  #p<-0.5
  #n<-15
  #reps<-20
  #a<-st<-rep(0,reps)
  #for (k in 1:reps) {
  #	m<-randDGDAG(n,p,DAG=TRUE)
  #	mc <- transitive.closure(m,mat=TRUE,loops=FALSE)
  #	mr <- transitive.reduction(mc)
  #	st[k]<-system.time(a[k]<-nle(mr))[1];
  #	showDAGTRTC(m=m,mr=mr,t=k)
  #}
  #plot(log(a),log(st))

  if (length(tr)==1) {return(1)}
  n<-dim(tr)[1]
  cs<-apply(tr,2,sum)
  if (sum(cs)==0) {return(factorial(n))}
  csi<-(cs==0)
  bs<-apply(tr,1,sum)
  bsi<-(bs==0)
  free<-which(bsi&csi)
  k<-length(free)
  if (k==n) {return(factorial(n))}
  if (k>0) {
    tr<-tr[-free,-free]
    cs<-apply(tr,2,sum)
    csi<-(cs==0)
    bs<-apply(tr,1,sum)
    bsi<-(bs==0)
    fac<-factorial(n)/factorial(n-k)
  } else {
    fac<-1
  }
  if ( (n-k)==2 ) {
    return(fac)
  }
  #if ( (n-k)==3 ) {
  #	return(fac*sum(tr[csi,]))
  #}
  tops<-which(csi)
  bots<-which(bsi)
  if (length(tops)==1 & length(bots)==1) {
    return(fac*nle(tr[-c(tops,bots),-c(tops,bots)]))
  }
  if (length(bots)<length(tops)) {tops<-bots}

  count<-0
  for (i in tops) {
    trr<-tr[-i,-i]
    count<-count+nle(trr)
  }
  return(fac*count)
}

p1<-function(tc,first.or.last='first') {
  #given PO=tc calculate the probability each node is 'first' or 'last'
  if ((dim(tc)[1]==1) && (dim(tc)[2]==1)) {
    weights=1
  } else {
    weights=0*tc[1,] #to pick up column names
    if (first.or.last=='first') {
      loc<-which(apply(tc,2,sum)==0) #top nodes
    } else {
      loc<-which(apply(tc,1,sum)==0) #bottom nodes
    }
    n.loc=length(loc)
    if (n.loc==1) {
      weights[loc]=1
    } else {
      for (k in 1:n.loc) {
        weights[loc[k]]=nle(tc[-loc[k],-loc[k]])
      }
    }
  }
  return(weights/sum(weights))
}

dag.conc2<-function(X,display.threshold=0.5,label.threshold=0.5,fontsize=8) {
  
  
  
  #compute reduction of concensus DAG from the n x n x T array
  
  #of adjacncy matrices X
  
  #include edges with support above display.threshold
  
  #label edges with support above label.threshold
  
  
  
  #example
  
  #p<-0.5
  
  #n<-8
  
  #reps<-100
  
  #X<-array(0,c(n,n,reps))
  
  #for (k in 1:reps) {
  
  #     X[,,k]<-randDGDAG(n,p,DAG=TRUE)
  
  #}
  
  #par(mfrow=c(1,2));
  
  #dag.con<-dag.conc(X,show.reduction=FALSE)
  
  #dag.con<-dag.conc(X)
  
  
  
  n<-dim(X)[1]
  
  who=colnames(X[,,1])
  
  dag.con<-apply(X,3,transitive.closure,mat=TRUE,loops=FALSE)
  
  #each sample matrix is turned into a vector - we hope columnwise
  
  dag.con<-matrix(apply(dag.con,1,mean),n,n,dimnames=list(who,who))
  
  #now written back into matrix columnwise - yuck, but seems OK
  
  dag.con2<-transitive.reduction(dag.con>display.threshold)
  
  dag.con<-round(dag.con,2)*dag.con2
  
  
  
  #if (sum(dag.con)==0) { #what on earth is this for??
  
  #  dag.con<-graph_from_adjacency_matrix(dag.con,mode ="directed")#as(dag.con,"graphNEL")
  
  #  #plot(dag.con,main=main) #how to plot in this case
  
  #  return(list(dag.con=dag.con))
  
  #} else {
  
  dag.con<-new("graphAM", adjMat=dag.con, edgemode="directed", values=list(weight=1))
  
  dag.con<-as(dag.con,"graphNEL") #graph_from_adjacency_matrix(dag.con,mode ="directed")#
  
  #}
  
  
  
  ewv<-unlist(edgeWeights(dag.con))
  
  #ew <- as.character(ewv); ew[ewv<label.threshold]<-""
  
  #ew <- ew[setdiff(seq(along = ew), removedEdges(dag.con))]
  
  #names(ew) <- edgeNames(dag.con)
  
  eAttrs <- list()
  
  #eAttrs$label <- ew
  
  #ft<-as.numeric(ew)*0+fontsize; names(ft)<-names(ew)
  
  #eAttrs$fontsize <- ft
  
  n.edges=length(edgeNames(dag.con))
  
  #e.style <- rep("dashed",n.edges); names(e.style)<-edgeNames(dag.con)
  
  #eAttrs$style<-e.style
  
  e.color <- rep("red",n.edges); names(e.color)<-edgeNames(dag.con)
  
  e.color[ewv<0.9]<-"black"
  
  eAttrs$color<-e.color
  
  
  
  #labeldistance <- rep(10,n.edges); names(labeldistance)<-edgeNames(dag.con)
  
  #eAttrs$labeldistance<-labeldistance
  
  
  
  #labelfloat<-rep("true",n.edges); names(labelfloat)<-edgeNames(dag.con)
  
  #eAttrs$labelfloat<-labelfloat
  
  
  
  nAttrs<-list()
  
  n.nodes=length(nodes(dag.con))
  
  nAttrs$height <- nAttrs$width <- rep("0", n.nodes)
  
  nAttrs$style <- rep("invis",n.nodes)
  
  nAttrs$fontsize <- rep(fontsize, n.nodes)
  
  nAttrs <- lapply(nAttrs, function(x) {names(x) <- nodes(dag.con); x })
  
  
  
  attrs <- list(node = list(shape = "plaintext", fixedsize = "false"), edge=list(headclip="false",tailclip="true",style="dashed",arrowsize=2.0), graph=list(splines="false"))
  
  #plot(dag.con,nodeAttrs=nAttrs,edgeAttrs=eAttrs,attrs=attrs,"dot") #how to plot
  
  
  
  list(dag.con=dag.con,nodeAttrs=nAttrs,edgeAttrs=eAttrs,attrs=attrs)
  
}

showDAGcon<-function(B,T,PO.con,years=1:T) {
  
  #display the concensus DAGs in PO.con arranged in a lattice and labeled by year
  
  n.years=length(years)
  
  n.row=min(4,floor(sqrt(n.years)))
  
  par(mfrow=c(n.row,ceiling(length(years)/n.row)),mai=c(0,0,0,0))
  
  for (t in years) {
    
    h=PO.con[[t]]
    
    numEdges=length(unlist(edgeL(PO.con[[t]]$dag.con)))
    
    if (numEdges==0) {
      
      a=plot(h$dag.con,attrs=h$attrs,"dot")
      
    } else {
      
      a=plot(h$dag.con,nodeAttrs=h$nodeAttrs,edgeAttrs=h$edgeAttrs,attrs=h$attrs,"dot")
      
    }
    
    text(a@boundBox@upRight@x/2,0.95*a@boundBox@upRight@y,as.character(B+t-1),col=2)
    
  }
  
}

plot.mean.PO = function(Result, burn.in = 10000, main, display.threshold=0.5,label.threshold=0.5,fontsize=8){
  # not finished yet
  n = nrow(Result$Z[[1]])
  n.record = Result$n.record
  Result.Z = Result$Z[(burn.in/n.record):length(Result$Z)]
  Result.PO = lapply(lapply(lapply(Result.Z, latent2order), order2partial, n = n), transitive.reduction)
  dag.con = Reduce('+', Result.PO)/length(Result.Z)
  # the following part is adapted from the dag.conc2 function
  dag.con2<-transitive.reduction(dag.con>display.threshold)
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

## Log-likelihoods

loglkd<-function(r=.5,mc,la,model='lkdup') {
  
  #loglikelihood r=prob for error, mc is closure of PO
  #la[[]] is a list, elements have a field which is an order sequence o
  
  ##example
  #wl<-list()
  #wl[[1]]<-list(o=c(1,2,3,4,5))
  #wl[[2]]<-list(o=c(1,3,2,4,5))
  
  ##mr is the DAG for the seq 1 and conflicts seq 2
  #mr<-seq2dag(wl[[1]]$o,5) #it is a reduction
  #mc<-transitive.closure(mr,mat=TRUE,loops=FALSE)
  #
  #loglkd(r=0.5,mc,wl,'lkdup')
  
  ##if the seq builder chooses each person at random
  #loglkd(r=1,mc,wl,'lkdup')
  #-log(factorial(5))
  
  ##if there is no order information, the seq builder
  ##places each person correctly, but this is the same
  ##as placing them at random
  #mcp<-matrix(0,5,5); colnames(mcp)<-row.names(mcp)<-c(1,2,3,4,5); 
  #loglkd(r=0.5,mcp,wl,'lkdup')
  
  if (model=='prior') {return(0)}
  
  lkddown.fac<-function(r,mc,o) {
    
    #evaluate the log-probability for the
    #placement of o[1] - the first person in the order o
    #given the PO with reduction mr
    
    #mr,o, an order of length 1, just one way to place it
    if (length(o)==1) {return(0)}
    
    #the sub-DAG for the order
    mla<-mc[o,o]
    #the sub-DAG with the first element removed
    mlb<-mc[o[-1],o[-1]]
    
    #first person may have been placed at random
    fac <- r/length(o)
    
    #if the first person is in a place that does
    #not violate the proposed PO (given by mc) then
    #they may have been placed using the distribution
    #over linear extensions
    if (sum(mla[,1])==0) {
      fac<-fac+(1-r)*nle(mlb)/nle(mla)
    }
    #return the log-likelihood for this placement
    #plus the log-likelihood for the subsequent placements
    return(log(fac)+lkddown.fac(r,mc,o[-1]))
    
  }
  
  lkdup.fac<-function(r,mc,o) {
    
    #evaluate the log-probability for the
    #placement of o[n] - the last person in the order o
    #given the PO with reduction mr
    
    #mr,o, an order of length 1, just one way to place it
    n<-length(o)
    if (n==1) {return(0)}
    
    #the sub-DAG for the order
    mla<-mc[o,o]
    #the sub-DAG with the last element removed
    mlb<-mc[o[-n],o[-n]]
    
    #last person may have been placed at random
    fac <- r/length(o)
    
    #if the last person is in a place that does
    #not violate the proposed PO (given by mc) then
    #they may have been placed using the distribution
    #over linear extensions
    if (sum(mla[n,])==0) {
      fac<-fac+(1-r)*nle(mlb)/nle(mla)
    }
    #return the log-likelihood for this placement
    #plus the log-likelihood for the subsequent placements
    return(log(fac)+lkdup.fac(r,mc,o[-n]))
    
  }
  
  lkdnat.fac<-function(r,mc,o) {
    
    #evaluate the log-probability for the
    #placement of o[1] - the first person in the order o
    #given the PO with reduction mr
    
    #mr,o, an order of length 1, just one way to place it
    n<-length(o)
    if (n==1) {return(0)}
    
    #the sub-DAG for the order
    mla<-mc[o,o]
    
    #first person may have been placed at random
    fac <- r/length(o)
    
    #if the last person is in a place that does
    #not violate the proposed PO (given by mc) then
    #they may have been placed at random from the legal placements
    if (sum(mla[,1])==0) {
      fac<-fac+(1-r)/sum(apply(mla,2,sum)==0)
    }
    #return the log-likelihood for this placement
    #plus the log-likelihood for the subsequent placements
    return(log(fac)+lkdnat.fac(r,mc,o[-1]))
    
  }
  
  n.order<-length(la)
  llkda<-matrix(0,1,n.order)
  ind<-match(1:n,as.numeric(rownames(mc))); 
  mcs<-mc[ind,ind];
  for (k in 1:n.order) {
    if (model=='lkddown') {llkda[1,k]<-lkddown.fac(r,mcs,la[[k]])} # removed la[[k]]$o - different setup
    if (model=='lkdup') {llkda[1,k]<-lkdup.fac(r,mcs,la[[k]])}
    if (model=='lkdnat') {llkda[1,k]<-lkdnat.fac(r,mcs,la[[k]])}
  }
  return(sum(llkda))
  
}

## No-conflict

noconflict<-function(mc,la) {
  #test if there are conflicts 
  stop('old noconflict had a bug - see noconflict2')
  n.order<-length(la)
  if (n.order==0) {return(TRUE)} 
  for (k in 1:n.order) {
    o<-la[[k]]$o
    ol<-la[[k]]$ll
    ind<-match(1:n,as.numeric(rownames(mc)));
    mcs<-mc[ind,ind]
    for (i in 2:ol) {
      if (mcs[o[i],o[i-1]]==1) {return(FALSE)} 
    }
  }
  TRUE
}

#the noconflict fn in pofun is dodgy - this is correct (I believe) ... see below
#old noconflict() just returns an error message saying it doesn't work
noconflict2<-function(h,y2l,cla,years=1:T) { #h must be trans.closed
  #test if there are conflicts
  for (t in years) {
    lat=cla[y2l[[t]]]
    mc=h[[t]]
    n.order<-length(lat)
    if (n.order>0) {
      for (k in 1:n.order) {
        o<-lat[[k]]$o
        ol<-lat[[k]]$ll
        ind=match(o,as.numeric(rownames(mc)))
        for (i in 1:(ol-1)) {
          for (j in (i+1):ol) {
            if (mc[ind[j],ind[i]]==1) return(list(outcome=FALSE,t=t,bl=y2l[[t]][k],b=o[c(j,i)]))
          }
        }
      }
    }
  }
  return(list(outcome=TRUE))
}

## Transfer b/w Objects

seq2dag<-function(o,n,p=1) {
  m<-matrix(0,n,n)
  for (k in 2:length(o)) {
    m[o[k-1],o[k]]<-p
  }
  m
}

latent2order<-function(z) {
  #input: z is n by k
  #output: k perms of 1:n equal rank in z columns
  if (is.matrix(z)){
    n<-dim(z)[1]
    k<-dim(z)[2]
    resm<-n-apply(z,2,rank,ties.method="random")+1
    resl<-list()
    for (i in 1:k) {
      resm[resm[,i],i]<-1:n
      resl[[i]]<-resm[,i]
    }
  } else {
    n=length(z)
    k=1
    resm<-n-rank(z,ties.method="random")+1
    resl<-list()
    resm[resm]<-1:n
    resl=list(resm)
  }
  resl
}

order2partial<-function(v,n=NULL) {
  #this
  #output is the transitive closure of 
  #the intersection of the list of complete orders v
  if (is.null(n)) {n<-max(v)}
  w<-lapply(v,seq2dag,n)
  x<-lapply(w,transitive.closure,mat=TRUE,loops=FALSE)
  z<-matrix(0,n,n); 
  colnames(z)<-rownames(z)<-1:n
  for (y in x) {z<-z+y}
  0+(z==length(v))
}

## Data Processing

makejkdata<-function(on=1122,off=1128,takebishops=TRUE,include.nodates=FALSE,include.singletons=FALSE,remove.badbishops=FALSE) {
  
  #read the data from files
  act <- read.table("act.txt", sep="\t", header = T)
  person <- read.table("person.txt", sep="\t", header = T)
  pnames <- read.table("person name.txt", sep="\t", quote = "\"", header = T, fill = T)
  
  #pull out the first name, second name and position information
  #most have just 1st and 2nd or 1st and 3rd, but see for example
  #Roger (1st), d'Abetot (2nd), sheriff of Worcester (3rd) entry 4563, id 4723
  z<-pnames[,c('name.entry.element','name.second.element','name.third.element')]
  z2<-apply(cbind(levels(z[,1])[z[,1]],levels(z[,2])[z[,2]],
                  levels(z[,3])[z[,3]]),1,paste,collapse=', ')
  z3<-sub(', ,',',', z2) 
  z4<-sub(', +$','', z3)
  z4n<-pnames$person.name.ID
  
  #a list of name-ids (not row #s!) corresponding to Bishops
  if (takebishops) {
    #select the bishops
    rw<-grep(c("^bishop"),pnames$name.third.element,ignore.case=TRUE)
    bip<-pnames$person.name.ID[rw]
    
    #standardise their dioceses
    z5<-pnames[,c('name.third.element')][rw]
    diocese<-sub("-[[:digit:]]+","",
                 sub(" [[:digit:]]+","",
                     sub(", [[:graph:]]+","",
                         sub("bishop[[:blank:]]+of ","",z5,ignore.case=TRUE))))
    diocese<-sub("St David's","St Davids",diocese)
    diocese<-sub("Chester-Coventry","Chester",diocese)
    diocese<-sub(" Stephen","",diocese)
    diocese<-sub("S[[:graph:]][[:graph:]][[:graph:]]$","Sees",diocese)
    #bip the person id's matching diocese - person bip[i] is in diocese[i]
    #diocese<-sort(unique(diocese))
    
    if (remove.badbishops) {
      
      bad.dio=c(grep("Lincoln/Chester",diocese),grep("Thetford/Norwich",diocese),
                grep("^bishop$",diocese),grep("uncertain",diocese),grep("attesting",diocese))
      diocese<-diocese[-bad.dio]
      bip<-bip[-bad.dio]
      rw<-rw[-bad.dio]
      
      #same thing - check
      #bad.pn=c(grep("Lincoln/Chester",pnames$name.third.element),grep("Thetford/Norwich",pnames$name.third.element),
      #         grep("^bishop$",pnames$name.third.element),
      #         grep("bishop, identity uncertain",pnames$name.third.element),grep("attesting",pnames$name.third.element))
      #bad.nid<-pnames[bad.pn,]$person.name.ID
      #bip<-setdiff(bip,bad.nid)
      
    }
    
    #GKN 1/1/19
    
  }
  
  #make a list of all the act id's that people appear in
  v<-unique(person$act.ID)
  la<-NULL
  count<-0
  
  for (i in v) {
    
    #find the act with id=i in the list of acts 
    a<-act[act$act.ID==i,]
    
    #if the year.2 field is empty or 0 set it equal the year.1 value
    #ie if the data is telling us it happened in a single year, make the
    #upper and lower limits equal
    if ( (is.na(a$year.2)) | (a$year.2==0) ) {a$year.2<-a$year.1}
    
    #get all the people who were witnesses in this act
    #z has columns  "person.name.ID" and "name.sequence" the latter being
    #their position in the list
    z<-person[(person$act.ID==i)&(person$name.role=='W'),c(2,4)]
    
    #Is this act of further interest? 
    #Were there any witnesses (dim...>1)
    #Are there any data on the time of the act (if year.1 is null forget it)?
    #Is the lower limit after 'on' and the upper limit before 'off'?
    if ( (dim(z)[1]>1 | include.singletons) & (!is.na(a$year.1) | include.nodates) & (is.na(a$year.1) | ((a$year.2<off) & (a$year.1>on))) ) {
      #the list order of the people in this witness list
      zs<-sort.int(z[,2],index.return=TRUE)
      #if (any(diff(zs$x)!=1)) {cat('this stinks \n'); break}
      if (takebishops) {
        #take the subset of witnesses who are bishops
        tem<-list(w=intersect(z[zs$ix,1],bip),id=i,year.1=a$year.1, year.2=a$year.2, ac=a$authenticity.code)
        tem$allnames<-unlist(lapply(z[zs$ix,1],function(x) {z4[z4n==x]}))
      } else {
        tem<-list(w=z[zs$ix,1],id=i,year.1=a$year.1, year.2=a$year.2, ac=a$authenticity.code)
      }
      tem$ll<-length(tem$w)
      tem$names<-unlist(lapply(tem$w,function(x) {z4[z4n==x]}))
      if (tem$ll>1) {
        count<-count+1
        la[[count]]<-tem
      }
    }
  }
  n.order<-count
  
  #make a visual display of the lists
  unlist(lapply(la,function(x) x$ll))
  
  #la$w has name id sorted by list location but these name id's are not packed
  #so gather together all the distinct name ids
  #that appear in our list, and give them a new, extra id
  #that just counts up to the number in our subset
  v2nid<-unique(unlist(lapply(la,function(x) x$w)))
  
  #mapping name id to vertex
  nid2v<-function(i) {
    els<-length(i)
    u<-rep(NaN,els)
    for (j in 1:els) {
      u[j]<-which(v2nid==i[j])
    }
    u
  }
  #use vertex id's in data
  for (count in 1:n.order) {la[[count]]$o<-nid2v(la[[count]]$w)}
  
  #return the selected and processed data
  if (takebishops) {
    bin=is.element(bip,v2nid); #which bishops made into the lists?
    b.name.id=bip[bin];        #thin the bishop id's
    diocese=diocese[bin]       #and matching dioceses
    rw=rw[bin]                 #the rows of pnames that ended up getting used
    return( list(la=la,
                 pnames=list(name.id=pnames[,1], name=z4),
                 qnames=list(name.id=b.name.id,name=z4[rw],node.name=nid2v(b.name.id),diocese=diocese))
            )
  } else {
    just.used=is.element(pnames[,1],v2nid)
    return( list(la=la,
                 pnames=list(name.id=pnames[,1], name=z4),
                 qnames=list(name.id=pnames[just.used,1],name=z4[just.used],node.name=nid2v(v2nid),diocese=NA))
            )
  } 
}

relabel = function(v){
  # relabel the actors to 1:n based on increasing order
  n = length(unique(unlist(v)))
  index_mat = cbind(sort(unique(unlist(v)), decreasing = FALSE), 1:n)
  v_melt = lapply(v, function(x, index_mat){
    for(i in 1:n){x[x==index_mat[i,1]]=index_mat[i,2]}
    return(x)
  }, index_mat = index_mat)
  return(list(list = v_melt, index_matrix = index_mat))
}

## MCMC - POPL

RandPar<-function(n, K, rho) {
  
  #generate a partial order based on (K, n, rho) partial order model
  #output is adjacency matrix
  
  Sigma = matrix(rep(rho, K*K), K, K)
  diag(Sigma) = 1
  Z = mvrnorm(n = n, mu = rep(0, K), Sigma)
  return(order2partial(v = latent2order(Z), n = n))
}

sample.vec <- function(x, ...) x[sample(length(x), ...)]

norAlpha = function(alpha.index, sigma, n){
  # Trivial function to help generate Alpha
  return(sort(rnorm(n, mean = 0, sd = sigma), decreasing = TRUE)[order(alpha.index)])
}

simAlpha = function(mc, sigma, N){
  
  # Simulate matrix Alpha (N x n)
  
  v = rep(list(mc), N)
  v = lapply(v, unifLE)
  
  while(!all(transitive.reduction(order2partial(v, n = nrow(mr)))==mr)){
    v = rep(list(mc), N)
    v = lapply(v, unifLE)
  }
  
  return(t(sapply(v, norAlpha, n = n, sigma = sigma)))
}

PO2Randlatent = function(mc, K){
  
  # Generates a random latent matrix that EXACTLY follows a given partial order

  n = nrow(mc)
  
  Z.index = rep(list(mc), K)
  Z.index = lapply(Z.index, unifLE)
  
  while(any(!(transitive.closure(order2partial(v = Z.index, n = n), mat = TRUE, loops = FALSE)==mc))){
    Z.index = rep(list(mc), K)
    Z.index = lapply(Z.index, unifLE)
  }
 
  Z.index = do.call(cbind, Z.index)
  Z = apply(Z.index, 2, function(v) {sort(rnorm(n), decreasing = TRUE)[order(v)]})
  return(Z)
}

seq2dag1 = function(o, n, p=1){
  m = matrix(0, n, n)
  m[t(combn(o, 2))] = p
  return(m)
}

order2partial1<-function(v,n) {
  #this
  #output is the transitive closure of 
  #the intersection of the list of incomplete orders v
  z<-Reduce('+', lapply(v,seq2dag1,n))
  z[z!=0] = 1
  z_temp = z
  z_temp[z_temp==0] = NA
  z[which((t(z_temp) - z_temp)==0, arr.ind = TRUE)] = 0
  colnames(z)<-rownames(z)<-1:n
  return(z)
}

## Conditional Posterior distributions (in log scale)

### Latent matrix Z

nle_denominator = function(data_vec, PO){
  return(nle(transitive.reduction(PO)))
}

dZ = function(Z, rho, log_nle){
  K = ncol(Z)
  n = nrow(Z)
  sigma = matrix(rep(rho, K*K), K, K)
  diag(sigma) = 1
  A = sum(apply(Z, 1, dmvnorm, sigma = sigma, log = TRUE))
  return(A-log_nle)
}

noconflict3 = function(AlphaPL, PO){
  # Tests whether the list AlphaPL corresponding to existing actors violates partial order PO
  # If TRUE -> no conflict
  ## Note: the input must be a list
  
  N = length(AlphaPL)
  test.result = rep(NA, N)
  
  for(i in 1:N){
    alpha.order = order(AlphaPL[[i]], decreasing = TRUE)
    index.matrix = t(combn(alpha.order, 2))[,c(2,1)]
    test.result[i] = (sum(transitive.closure(PO,mat=TRUE,loops=FALSE)[index.matrix]) == 0)
  }
  
  return(all(test.result))
}

noconflict4 = function(LEs, PO){
  # Test whether linear extensions violate partial order PO
  # If TRUE -> no conflict
  
  noconflict_le = function(LE, PO){
    index.matrix = t(combn(length(LE), 2, function(x) LE[x]))[,c(2,1)]
    return(sum(transitive.closure(PO, mat = TRUE, loops = FALSE)[index.matrix])==0)
  }
  return(all(unlist(lapply(LEs,noconflict_le, PO))))
}

### Actor Attributes Alpha

lglk = function(Alpha, Y){
  
  lglk_vec = function(alpha, y){
    sum(alpha[y]) - sum(log(revcumsum(exp(alpha[y]))))
  }
  
  if (is.matrix(Alpha)){
    N = nrow(Alpha)
    lglk = 0
  
    for(j in 1:N){
      lglk = lglk + lglk_vec(Alpha[j,], Y[[j]])
    }
    return(lglk)} else {
      
      return(lglk_vec(Alpha, Y))
    }
}

dAlpha = function(Alpha, sigma, Y){
  
  return(sum(dnorm(Alpha, sd = sigma, log = TRUE)) + lglk(Alpha, Y))
}

### Depth Parameter rho

drho = function(rho, Z){
  K = ncol(Z)
  Sigma = matrix(rep(rho, K*K), K, K)
  diag(Sigma) = 1
  return(sum(apply(Z, 1, dmvnorm, sigma = Sigma, log = TRUE)) + dbeta(rho, 1, 1/6, log = TRUE)) 
}

# Sampling (GS for sigma)

gibbs.sampler = function(Z_initial, rho_initial, alpha_initial, sigma_initial, sigma_prior = 5, 
                         sigma_prior_sd = 7, alpha_sd = 9, w_alpha = 0.7, data, n.sims, 
                         n.write = 100000, write.file.main, wd = "~/Desktop/POPLModel", return = FALSE){
  
  ptm = proc.time()
  
  K = ncol(Z_initial)
  n = nrow(Z_initial)
  N = length(data)
  
  n.record = max(n^2, n*N)
  
  # shape and scale for sigma
  shape_initial = (sigma_prior^4)/(sigma_prior_sd^2) + 2
  scale_initial = (sigma_prior^2)*(shape_initial-1)
  
  Z_list = vector(mode= "list", length = n.sims/n.record)
  rho_list = vector(mode="numeric", length = n.sims/n.record)
  sigma_list = vector(mode="numeric", length = n.sims/n.record)
  Alpha_list = vector(mode= "list", length = n.sims/n.record)
  loglik = vector(mode="numeric", length = (n.sims/n.record-1))
  
  # record acceptances of each parameter
  acceptance_Z = rep(0, n.sims-1)
  acceptance_rho = rep(0, n.sims-1)
  acceptance_Alpha = rep(0, n.sims-1)
  
  # record rejection messages
  rejection_Z = rep(1, n.sims-1)
  rejection_Alpha = rep(1, n.sims-1)
  
  # initial values
  Z = Z_initial
  rho = rho_initial
  sigma = sigma_initial
  Alpha = alpha_initial
  PO = order2partial(v = latent2order(Z_initial), n = nrow(Z_initial))
  log_nle = sum(log(unlist(lapply(data, nle_denominator, PO))))
  
  Z_list[[1]] = Z
  rho_list[[1]] = rho
  sigma_list[[1]] = sigma
  Alpha_list[[1]] = Alpha
  
  for (i in 2:n.sims){
    
    # Update for Z
    
    Z_temp = Z
    rZ = sample(1:n, 1)
    cZ = sample(1:K, 1)
    Sigma = matrix(rep(rho, K*K), K, K)
    diag(Sigma) = 1
    Z_temp[rZ,cZ] = rcmvnorm(n = 1, mean = Z[rZ,], sigma = Sigma, dependent.ind = cZ, given.ind = (1:K)[-cZ], X.given = Z_temp[rZ,-cZ])
    PO_temp = order2partial(v = latent2order(Z_temp), n = nrow(Z_temp))
    if(noconflict3(Alpha, PO_temp)){
      log_nle_temp = N*(log(nle(transitive.reduction(PO_temp))))
      log_eta1 = dZ(Z_temp, rho, log_nle_temp) - dZ(Z, rho, log_nle)
      if (log_eta1 > log(runif(1))){Z = Z_temp; acceptance_Z[i-1] = 1; log_nle = log_nle_temp; PO = PO_temp} else {rejection_Z[i-1] = "Likelihood"}
    } else {rejection_Z[i-1] = "Conflict_w_Alpha"}
    # print(Z) 
    
    # Update for Alpha
    
    ## Stage one - update one element in Alpha
    Alpha_temp = Alpha
    p = sample(1:N, 1)
    q = sample(1:n, 1)
    Alpha_temp[p,q] = Alpha[p,q] + rnorm(1, 0, alpha_sd)
    
    if (noconflict3(Alpha = Alpha_temp[p,], PO = PO)) { # need to modify to consider one actor
      log_eta2 = dAlpha(Alpha_temp, sigma = sigma, Y = data)-dAlpha(Alpha, sigma = sigma, Y = data)
      if (log_eta2 > log(runif(1))) {Alpha = Alpha_temp; acceptance_Alpha[i-1] = 1} else {rejection_Alpha[i-1] = "Likelihood"}
    }  else {rejection_Alpha[i-1] = "Conflict_w_PO"}
    
    ## Stage two - scale each vector in Alpha
    for (j in 1:N){
      Alpha_temp = Alpha
      w2= w_alpha
      delta2 = runif(1, w2, 1/w2)
      Alpha_mean_j = mean(Alpha_temp[j,], na.rm= TRUE)
      Alpha_temp[j,] = Alpha_mean_j + delta2*(Alpha_temp[j,] - Alpha_mean_j) 
      log_eta4 = dAlpha(Alpha_temp, sigma = sigma, Y = data)+(n-3)*log(delta2)-dAlpha(Alpha, sigma = sigma, Y = data)
      if (log_eta4 > log(runif(1))) {Alpha = Alpha_temp}
    }
    # print(Alpha)
    
    # Update for rho
    
    w = 0.50
    delta = runif(1, w, 1/w)
    rho_temp = 1 - delta*(1-rho)
    log_eta3 = drho(rho_temp, Z)-drho(rho, Z)-log(delta)
    if(log_eta3 > log(runif(1)) & rho_temp<(1-10^(-7))){rho = rho_temp; acceptance_rho[i-1] = 1} 
    # print(rho)
    
    # Update for sigma
    shape = shape_initial + length(unlist(Alpha))/2
    scale = scale_initial + .5*sum(unlist(Alpha)^2)
    sigma = sqrt(rinvgamma(n = 1, shape = shape, scale = scale))
    # print(sigma) 
    
    if (i %% n.record == 0){
      Z_list[[i/n.record+1]] = Z
      rho_list[i/n.record+1] = rho
      sigma_list[i/n.record+1] = sigma
      Alpha_list[[i/n.record+1]] = Alpha
      loglik[i/n.record] = dAlpha(Alpha, sigma, data)
    }
    
    if (return == FALSE){
      if (i %% n.write == 0){
        setwd(wd)
        Result = list(n.record = n.record, Z = Z_list[1:(i/n.record+1)], rho = rho_list[1:(i/n.record+1)], sigma = sigma_list[1:(i/n.record+1)], Alpha = Alpha_list[1:(i/n.record+1)], loglik = loglik[1:(i/n.record+1)], acceptance = list(Z = acceptance_Z[1:(i/n.record+1)], rho = acceptance_rho[1:(i/n.record+1)], Alpha = acceptance_Alpha[1:(i/n.record+1)]), rejection = list(Z = rejection_Z[1:(i/n.record+1)], Alpha = rejection_Alpha[1:(i/n.record+1)]))
        save(Result, file = paste0(write.file.main, i/n.write, ".RData"))
      }
    } 
  }
  print(proc.time() - ptm)
  if (return == TRUE){
    return(list(n.record = n.record, Z = Z_list, rho = rho_list, sigma = sigma_list, Alpha = Alpha_list, loglik = loglik, acceptance = list(Z = acceptance_Z, rho = acceptance_rho, Alpha = acceptance_Alpha), rejection = list(Z = rejection_Z, Alpha = rejection_Alpha)))
  }
}
