
#
# VSP helper functions - GKN 11-05-22
#

node<-function(actor=NA,parent=NA,child=NA,type='',order='',nc=NA) {
  list(actor=actor,parent=parent,child=child,type=type,order=order,nc=nc)
}

tree2bin<-function(tree,use.node.names=FALSE) {
  nv=length(tree)
  h=matrix(0,nv,nv)
  nm=rep(NA,nv)
  for (i in 1:nv) {
    if (is.anst(tree[[i]])) {
      child=tree[[i]]$child
      h[i,child]<-1
      if (use.node.names) {
        nm[i]=i
      } else {
        nm[i]=paste(tree[[i]]$type, tree[[i]]$order, sep='') #substitute(list(a[b]),list(a=tree[[i]]$type,b=tree[[i]]$order))  #
      }
    } 
    if (is.leaf(tree[[i]])) {
      if (use.node.names) {
        nm[i]=i
      } else {
        nm[i]=paste(tree[[i]]$actor, tree[[i]]$order, sep='') #substitute(list(a[b]),list(a=tree[[i]]$actor,b=tree[[i]]$order))
      }
    }
  }
  
  row.names(h)<-colnames(h)<-nm

  #when deleting leaves we create zombie nodes - dont include these
  del=which(sapply(tree, function(x) is.zombie(x)))
  if (length(del)>0) h<-h[-del,-del]

  return(h)
}

showTREE<-function(tree=NULL,use.node.names=FALSE,...) {
  m=tree2bin(tree,use.node.names)
  g<-graph_from_adjacency_matrix(m,mode ="directed")#as(m,"graphNEL")
  #h<-graph_from_adjacency_matrix(my.transitive.closure(m),mode ="directed") #ed 5-4-22 for mnem
  h<-g #seems like a bug in layout function puts vertices on top of one another
  plot(g,layout=layout_as_tree(h),...);     
}

is.leaf<-function(node) { (!is.zombie(node) && is.na(node$child[1])) } 
is.anst<-function(node) { (!is.zombie(node) && !is.na(node$child[1])) }
is.root<-function(node) { (!is.zombie(node) && is.na(node$parent))}
is.zombie<-function(node) {identical(node,list())}
nodes.under<-function(tree,v){
  # returns the indices of the nodes under v (including v)
  if(all(is.na(tree[[v]]$child))) return(v)
  nodes = c(v,nodes.under(tree,tree[[v]]$child[1]),nodes.under(tree,tree[[v]]$child[2]))
  return(nodes)
}

is.vsp<-function(tc) {
  # function to detect if a transitive closure is a series-parallel partial order
  tcg<-graph_from_adjacency_matrix(tc,mode ="directed")
  b<-matrix(c(0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0),4,4,byrow=TRUE)
  bg<-graph_from_adjacency_matrix(b,mode ="directed")
  return(!subgraph_isomorphic(bg,tcg,induced=TRUE))
}

po2tree<-function(h,tree=vector('list',0)) {
  
  # generates a tree (vsp) from a partial order
  # 
  # parameter
  # ---------
  # h: matrix. 
  #   The transitive closure of the partial order. 
  # tree: tree (vsp), default tree = vector('list',0).
  #   The working tree (vsp). 
  #
  # return
  # ------
  #   list. 
  #
  #   Attribute
  #   ---------
  #     h: matrix.  
  #       The unstandardised partial order (transitive closure).
  #     tree: tree (vsp). 
  #       The tree (vsp) corresponding to the input partial order. 
  # 
  # example
  # -------
  # N=6; actors=1:N; tree=rvsp(actors,0.5); h=tree2po(tree)
  # par(mfrow=c(1,2)); showDAG(transitive.reduction(h),edge.arrow.size=3/N); showTREE(tree,edge.arrow.size=3/N)

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

is.cherry<-function(h,i,j) {
  a=all(h[i,-c(i,j),drop=FALSE]==h[j,-c(i,j),drop=FALSE])
  b=all(h[-c(i,j),i,drop=FALSE]==h[-c(i,j),j,drop=FALSE])
  return(a&b)
}

tree2po<-function(tree,r=NA) {
  
  # generates the corresponding transitive closure from a vsp tree
  # 
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The vsp tree of interest with N actors. 
  # r: int., default r = NA. 
  #   The node index of the root. 
  # 
  # return
  # ------
  #   matrix (N x N). 
  #   The transitive closure of the corresponding partial order. 

  if (is.na(r)) r=which(sapply(tree,function(x) is.root(x)))

  if (is.leaf(tree[[r]])) {
    h=matrix(0,1,1); row.names(h)<-colnames(h)<-tree[[r]]$actor
    return(h)
  }

  c1=tree[[r]]$child[1]
  c2=tree[[r]]$child[2]
  if(tree[[c1]]$order=='+'|tree[[c1]]$order==''){h1=tree2po(tree,c1)} else {h1=tree2po(tree,c2)}
  if(tree[[c2]]$order=='-'|tree[[c2]]$order==''){h2=tree2po(tree,c2)} else {h2=tree2po(tree,c1)}
  
  n1=dim(h1)[1]
  n2=dim(h2)[1]
  is.S=(tree[[r]]$type=='S')+0
  ub=matrix(is.S,n1,n2)
  lb=t(0*ub)
  h.ub=cbind(h1,ub)
  h.lb=cbind(lb,h2)
  h=rbind(h.ub,h.lb); 
  row.names(h)<-colnames(h)<-c(row.names(h1),row.names(h2))

  return(standardise(h))
}

rvsp.tree<-function(actors,prob.snode=0.5) {
  
  # creates a random vsp tree
  #
  # parameter
  # ---------
  # actors: vector, int. 
  #   A vector with integer values. The actors of interest, e.g. actors = 1:N. 
  # prob.snode: float. in [0,1]. 
  #   The probability of getting a type = "S" (Serial) internal node. 
  #
  # return
  # ------
  #   Tree (vsp), list of length 2N-1.
  #
  #   Attribute
  #   ---------
  #   actor: int. 
  #     The actor represented by node i (i in 1:N). If internal node, actor = NA. 
  #   parent: int. 
  #     The node index of parent node. If root, parent = NA. 
  #   child: int. 
  #     The node index of child node. If leaf node, child = NA. 
  #   type: character. 
  #     The type of operation on internal node, type = "S" (Serial) or "P" (Parallel). If leaf node, type = "". 
  #   order: character. 
  #     The order for serial operation, order = "+" or "-". Order "+" nodes beat the order "-" nodes. If root or parented by "P" nodes, order = "". 
  #   nc: int.
  #     The number of leaves covered by the node. A leaf node covers one leaf. 
  
  n=length(actors)
  tree=vector('list',2*n-1)
  tree[[1]]=node(actor=actors[1],parent=3)
  tree[[2]]=node(actor=actors[2],parent=3)
  tree[[3]]=node(child=c(1,2))
  v=3; root=3
  if (n>2) {
    for (k in 3:n) {
      a=sample(1:v,1)
      new.leaf=v+1;
      new.anst=v+2
      tree[[new.leaf]]=node(actor=actors[k],parent=new.anst)
      tree[[new.anst]]=node(child=c(a,new.leaf))
      if (a==root) {
        root=new.anst
      } else {
        ap=tree[[a]]$parent
        ca=which(tree[[ap]]$child==a)
        tree[[ap]]$child[ca]=new.anst
        tree[[new.anst]]$parent=ap
      }
      tree[[a]]$parent=new.anst
      v=v+2
    }
  }

  for (i in 1:(2*n-1)) {
    if (is.anst(tree[[i]])) {
      node.SP=sample(c('S','P'),1,replace=FALSE,prob=c(prob.snode,1-prob.snode))
      tree[[i]]$type=node.SP
      if (node.SP=='S') {
        o=sample(tree[[i]]$child);
        tree[[i]]$child=o
        a=o[1]; b=o[2]; tree[[a]]$order="+"; tree[[b]]$order="-"
      }
    }
  }
  # add children counts
  v=which(sapply(tree,function(x) is.root(x)))
  cc=get.child.count(tree,v)
  n=dim(cc)[2]
  for (i in 1:n) tree[[cc[1,i]]]$nc=cc[2,i]
  return(tree)
}

rvsp.old<-function(actors,p=0.5) {
  #samples ordered histories UAR not topologies
  n=length(actors)
  tree=vector('list',2*n-1)
  for (i in 1:n) {tree[[i]]=node(actor=actors[i])}

  L=1:n; k=length(L); v=n
  while (k>1) {
    v=v+1
    m=sample(L,2); a=m[1]; b=m[2]
    L=c(setdiff(L,m),v)
    tree[[a]]$parent=tree[[b]]$parent=v
    node.SP=sample(c('S','P'),1,replace=FALSE,prob=c(p,1-p))
    tree[[v]]=node(child=c(a,b),type=node.SP)
    if (node.SP=='S') {
      o=sample(c(a,b)); tree[[v]]$child=o
      a=o[1]; b=o[2]; tree[[a]]$order="+"; tree[[b]]$order="-"
      #if (o[1]>o[2]) tree<-swap.nodes(tree,o)
    }
    k=length(L)
  }
  return(tree)
}

formPclusters<-function(r=NA) {

  if (is.na(r)) r=which(sapply(G_tree,function(x) is.root(x)))

  if (is.leaf(G_tree[[r]])) {G_tree[[r]]$Pclust<<-list(); return()}

  c1=G_tree[[r]]$child[1]
  c2=G_tree[[r]]$child[2]
  formPclusters(c1)
  formPclusters(c2)
  cl1=G_tree[[c1]]$Pclust
  cl2=G_tree[[c2]]$Pclust

  if (G_tree[[r]]$type=='P') {
    c1.merge=c2.merge=FALSE
    if (length(cl1)>0) {
      tc1=sapply(cl1,function(x) any(x==c1))
      if (any(tc1)) {c1.merge=TRUE; li1=which(tc1)}
    }
    if (length(cl2)>0) {
      tc2=sapply(cl2,function(x) any(x==c2))
      if (any(tc2)) {c2.merge=TRUE; li2=which(tc2)}
    }
    if (!c1.merge & !c2.merge) {G_tree[[r]]$Pclust<<-c(cl1,cl2,list(r))}
    if (c1.merge & !c2.merge) {cl1[[li1]]=c(cl1[[li1]],r); G_tree[[r]]$Pclust<<-c(cl1,cl2)}
    if (!c1.merge & c2.merge) {cl2[[li2]]=c(cl2[[li2]],r); G_tree[[r]]$Pclust<<-c(cl1,cl2)}
    if (c1.merge & c2.merge) {a=cl1[[li1]]; b=cl2[[li2]]; cl1<-cl1[-li1]; cl2<-cl2[-li2]; nc=c(a,b,r); G_tree[[r]]$Pclust<<-c(cl1,cl2,list(nc))}
  } else {G_tree[[r]]$Pclust<<-c(cl1,cl2)}
  return()   
} 

form.type.clusters<-function(r=NA,type) {

  if (is.na(r)) r=which(sapply(G_tree,function(x) is.root(x)))

  if (is.leaf(G_tree[[r]])) {G_tree[[r]]$clust<<-list(); return()}

  c1=G_tree[[r]]$child[1]
  c2=G_tree[[r]]$child[2]
  form.type.clusters(c1,type)
  form.type.clusters(c2,type)
  cl1=G_tree[[c1]]$clust
  cl2=G_tree[[c2]]$clust

  if (G_tree[[r]]$type==type) {
    c1.merge=c2.merge=FALSE
    if (length(cl1)>0) {
      tc1=sapply(cl1,function(x) any(x==c1))
      if (any(tc1)) {c1.merge=TRUE; li1=which(tc1)}
    }
    if (length(cl2)>0) {
      tc2=sapply(cl2,function(x) any(x==c2))
      if (any(tc2)) {c2.merge=TRUE; li2=which(tc2)}
    }
    if (!c1.merge & !c2.merge) {G_tree[[r]]$clust<<-c(cl1,cl2,list(r))}
    if (c1.merge & !c2.merge) {cl1[[li1]]=c(cl1[[li1]],r); G_tree[[r]]$clust<<-c(cl1,cl2)}
    if (!c1.merge & c2.merge) {cl2[[li2]]=c(cl2[[li2]],r); G_tree[[r]]$clust<<-c(cl1,cl2)}
    if (c1.merge & c2.merge) {a=cl1[[li1]]; b=cl2[[li2]]; cl1<-cl1[-li1]; cl2<-cl2[-li2]; nc=c(a,b,r); G_tree[[r]]$clust<<-c(cl1,cl2,list(nc))}
  } else {G_tree[[r]]$clust<<-c(cl1,cl2)}
  return()   
} 

form.type.clusters.new<-function(tree,r=NA,type) {

  if (is.na(r)) r=which(sapply(tree,function(x) is.root(x)))

  if (is.leaf(tree[[r]])) {return(list())}

  c1=tree[[r]]$child[1]
  c2=tree[[r]]$child[2]
  cl1=form.type.clusters.new(tree,c1,type)
  cl2=form.type.clusters.new(tree,c2,type)
  
  if (tree[[r]]$type==type) {
    c1.merge=c2.merge=FALSE
    if (length(cl1)>0) {
      tc1=sapply(cl1,function(x) any(x==c1))
      if (any(tc1)) {c1.merge=TRUE; li1=which(tc1)}
    }
    if (length(cl2)>0) {
      tc2=sapply(cl2,function(x) any(x==c2))
      if (any(tc2)) {c2.merge=TRUE; li2=which(tc2)}
    }
    if (!c1.merge & !c2.merge) {clust<-c(cl1,cl2,list(r))}
    if (c1.merge & !c2.merge) {cl1[[li1]]=c(cl1[[li1]],r); clust<-c(cl1,cl2)}
    if (!c1.merge & c2.merge) {cl2[[li2]]=c(cl2[[li2]],r); clust<-c(cl1,cl2)}
    if (c1.merge & c2.merge) {a=cl1[[li1]]; b=cl2[[li2]]; cl1<-cl1[-li1]; cl2<-cl2[-li2]; nc=c(a,b,r); clust<-c(cl1,cl2,list(nc))}
  } else {clust<-c(cl1,cl2)}
  return(clust)   
} 

POvsp.prior<-function(PO,prob.snode=0.5) {
  tree=po2tree(PO)$tree
  return(TreeWeight(tree,prob.snode))
}

Ntrees<-function(x) {if (x<3) {return(1)} else { return((2*x-3)*Ntrees(x-1))}} 

TreeWeight<-function(tree,prob.snode=0.5) {
  root=which(sapply(tree,function(x) is.root(x)))

  Pclust=form.type.clusters.new(tree,root,'P')
  Sclust=form.type.clusters.new(tree,root,'S')
  
  Pcount=sum(sapply(tree,function(x) x$type=='P'))
  Scount=sum(sapply(tree,function(x) x$type=='S'))
  
  if (length(Pclust)==0) {
    Pf=1
  } else {
    NP=1+sapply(Pclust,length)
    Pf=prod(sapply(NP,Ntrees))
    if (!(Pcount=sum(sapply(Pclust,length)))) stop('error Pcount not matching')
  }
  if (length(Sclust)==0) {
    Sf=1
  } else {
    NS=sapply(Sclust,length)
    Sf=prod(sapply(NS,Catalan))  #prod(2^(NS-1)) #Pell(Scount+1)
    if (!(Scount=sum(sapply(Sclust,length)))) stop('error Scount not matching')
  }
  So=1/2^Scount
  Prob=Pf*Sf*So*prob.snode^Scount*(1-prob.snode)^Pcount/Ntrees( Pcount+Scount+1 )
  return(Prob)
}

Pell<-function(s) {
  if (s==0) return(0)
  if (s==1) return(1)
  return(2*Pell(s-1)+Pell(s-2))
}

Catalan<-function(s) {
  return(choose(2*s,s)/(s+1))
}

standardise<-function(PO) {
  
  # standardise a partial order adjacency matrix
  #
  # parameter
  # ---------
  # PO: acjacency matrix. 
  #   The partial order adjacency matrix. 
  # 
  # return
  # ------
  #   matrix. 
  #   The standardised adjacency matrix with row and columns sorted.
  
  o=order(as.numeric(row.names(PO)))  
  return(PO[o,o])
}

get.child.count<-function(tree,i=NA) {
  #get a 2 x length(tree) matrix (not counting zombie nodes)
  # - top row is the node index in tree[] 
  # - bot row is the number of leaves that node covers.
  # - a leaf covers one leaf

  if (is.na(i)) i=which(sapply(tree,function(y) is.root(y)))

  if (is.anst(tree[[i]])) {
    gcc1<-get.child.count(tree,tree[[i]]$child[1])
    gcc2<-get.child.count(tree,tree[[i]]$child[2])
    gcci<-matrix(c(i,gcc1[2,1]+gcc2[2,1]),2,1)
    gcc<-cbind(gcci,gcc1,gcc2)
  } else {
    gcc=matrix(c(i,1),2,1)
  }
  return(gcc)
}

#not used
insert.leaf.cover<-function(x,i=NA) {
  #WARNING! This is a side effect function that writes into each node
  #of tree in the calling environment the labels of the leaves it covers
  #a leaf covers one leaf
  
  if (is.na(i)) i=which(sapply(x,function(y) is.root(y)))

  if (is.anst(x[[i]])) {
    tree[[i]]$leaf.cover<<-c(insert.leaf.cover(x,x[[i]]$child[1]),insert.leaf.cover(x,x[[i]]$child[2]))
  } else {
    tree[[i]]$leaf.cover<<-i; 
  }
  return(tree[[i]]$leaf.cover)
}

nle.tree <- function(tree,v=NA){
  # count the number of le's of the PO expressed by tree 
  
  if (is.na(v)) {
    v=which(sapply(tree,function(x) is.root(x)))
    # cc=get.child.count(tree,v)
    # n=dim(cc)[2]
    # for (i in 1:n) tree[[cc[1,i]]]$nc=cc[2,i]
  }
  
  if (is.leaf(tree[[v]])) {
    return(1)
  }
  
  c1=tree[[v]]$child[1]
  c2=tree[[v]]$child[2]
  nle1=nle.tree(tree,c1)
  nle2=nle.tree(tree,c2)
  tot=nle1*nle2
  
  if (tree[[v]]$type=='P') {    
    nc1=tree[[c1]]$nc
    nc2=tree[[c2]]$nc
    tot=tot*choose(nc1+nc2,nc1)
  }
  
  return(tot)
}

nle.tree.old <- function(tree,v=NA){
  # count the number of le's of the PO expressed by tree 
  
  if (is.na(v)) {
    v=which(sapply(tree,function(x) is.root(x)))
    cc=get.child.count(tree,v)
    n=dim(cc)[2]
    for (i in 1:n) tree[[cc[1,i]]]$nc=cc[2,i]
  }

  if (is.leaf(tree[[v]])) {
    return(1)
  }

  c1=tree[[v]]$child[1]
  c2=tree[[v]]$child[2]
  nle1=nle.tree(tree,c1)
  nle2=nle.tree(tree,c2)
  tot=nle1*nle2

  if (tree[[v]]$type=='P') {    
    nc1=tree[[c1]]$nc
    nc2=tree[[c2]]$nc
    tot=tot*choose(nc1+nc2,nc1)
  }

  return(tot)
}

QP.LEProb.vsp<-function(tree,le,p,q) {

  #if tree is a vsp and le is one list then calculate the likeligood for the list
  #given the vsp-PO tree. Here p (err prob) and q (prob choose top down at a given le entry insertion)
  #call this on the suborder so le and actors.in.tree have same content
  
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

delete<-function(tree,i) {
  
  # deletes leaf node i from tree
  # reconnect tree - notice that the deleted leaf and parent node stay in tree-list but are unconnected zombies
  # 
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The vsp tree of interest. 
  # i: int. 
  #   The index of a leaf (!) node. 
  #
  # return
  # ------
  #   tree (vsp). 
  #   The vsp tree with leaf node i removed - i and its parent node stays as zombie nodes (list(NULL)). 
  
  pd=tree[[i]]$parent
  pp=tree[[pd]]$parent
  cd.i=which(tree[[pd]]$child!=i)
  cd=tree[[pd]]$child[cd.i]
  
  tree[[cd]]$parent=pp
  tree[[cd]]$order=tree[[pd]]$order
  if (!is.na(pp)) {
    pp.ci=which(tree[[pp]]$child==pd)
    tree[[pp]]$child[pp.ci]=cd
    parents.i = find.parents(tree,pp)
    for(parent in parents.i){
      tree[[parent]]$nc = tree[[parent]]$nc - 1
    }
  }
  tree[[pd]]=tree[[i]]=list() #make the deleted nodes zombies
  return(tree)
}

is.top.bot<-function(tree,i,extreme='top') {
  if (extreme=='top') {od='-'} else {od='+'} 
  #test if a leaf node of tree is a top or bottom node in the PO 
  if (!is.leaf(tree[[i]])) stop('is.top.bot input node i is not a leaf')
  if (is.root(tree[[i]])) return(TRUE) #tree has just one leaf which is also root
  j=i
  if (tree[[j]]$order==od) return(FALSE)
  while (!is.na(tree[[j]]$parent)) {
    j=tree[[j]]$parent
    if (tree[[j]]$order==od) return(FALSE)    
  }
  return(TRUE)
}

find.parents <- function(tree,i){
  # returns the index of all nodes above node index i (including i)
  if (is.na(i)) stop('The index i is NA.')
  parents.i = i
  while(!is.na(tree[[i]]$parent)){
    parents.i = c(parents.i,tree[[i]]$parent)
    i = tree[[i]]$parent
  }
  return(parents.i)
}

find.actor <- function(tree,actor){
  # returns the index of an actor name 
  return(which(sapply(tree,function(node){if(is.zombie(node)){return(NA)} else {return(node$actor)}})==actor))
}

sub.tree <- function(tree,y){
  
  # creates a subtree given data list y
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The vsp tree of interest. 
  # y: int, vector. 
  #   A data list representing the order relation between actors. 
  #
  # return
  # ------
  #   tree (vsp). 
  #   The subtree with unrelated nodes removed (zombie nodes). 
  
  n = (length(tree)+1)/2
  unused = setdiff(1:n,y)
  
  for (actor in unused){
    i = find.actor(tree,actor)
    tree = delete(tree,i)
  }
  return(tree)
}

noconflict.tree <- function(tree,y){
  
  # tests if a data list y violates the vsp tree
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The vsp tree of interest. 
  # y: int, vector. 
  #   A data list representing the order relation between actors. 
  # 
  # return
  # ------
  #   bool. 
  #   If the data list y violates the vsp tree, TRUE for no violation and FALSE for violation. 
  
  tree = sub.tree(tree,y)
  i = find.actor(tree,y[1])
  test = is.top.bot(tree,i)
  
  if (!test) return(FALSE)
  
  while(length(y)>1){
    tree = delete(tree,i)
    y = y[-1]
    i = find.actor(tree,y[1])
    test = is.top.bot(tree,i)
    if (!test) return(FALSE)
  }
  return(TRUE)
}


