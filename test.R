# An example script to run the code.

source("~/Desktop/POwTie/code/source-code/pofun.R")
source("~/Desktop/vsp/Code/vspfun.R")
# source("~/Desktop/POwTie/code/mallows.R")
source("~/Desktop/vsp-hierarchy/code/vsp_hierarchy_func.R")
source("~/Desktop/vsp-hierarchy/code/vsp_hierarchy_mcmc.R")

dir='~/Desktop/vsp-hierarchy/output'
setwd(dir)

RHO_HYPERPAMA = 1/4
THETA_HYPERPAMA = 5    # MALLOWS
R_HYPERPAMA = 1.5  # QUEUE-JUMPING
RHO_PROPOSAL = 1/2
NUM_ACTORS = 5
NUM_ASSESSORS = 3
NUM_ORDERS = c(10,12,8)

# Simulation

# po0=matrix(c(0,1,0,1,rep(0,3),1,rep(0,6),1,rep(0,4),1,rep(0,5)),nrow=5,byrow=TRUE)
# rownames(po0)=colnames(po0)=1:5
# showDAG(po0)
# U0=PO2Randlatent(transitive.closure(po0,mat=TRUE,loops=FALSE),K=2)
# showDAG(transitive.reduction(u2po(U0)))
# tau=0.5
# rho=0.9
# Sigma=matrix(rho,nrow=2,ncol=2); diag(Sigma)=1
# U_error = lapply(rep(NUM_ACTORS,NUM_ASSESSORS),function(n) rmvnorm(n,sigma=(1-tau^2)*Sigma))
# U = lapply(U_error, function(ue) tau*U0+ue)
# while (!all(sapply(lapply(U,u2po),function(po){po2tree(po)$is_vsp}))){
#   U0=PO2Randlatent(transitive.closure(po0,mat=TRUE,loops=FALSE),K=2)
#   U = lapply(U_error, function(ue) tau*U0+ue)
# }
# 
# tcs = lapply(U, u2po)
# trs = lapply(tcs, transitive.reduction)
# showDAG(trs[[1]])
# showDAG(trs[[2]])
# showDAG(trs[[3]])
# 
# simulation = data_simulation(NUM_ACTORS, NUM_ASSESSORS, NUM_ORDERS, U0=U0, U=U)
# Y = simulation$Y
# save(simulation,file='simulation1.RData')

load('simulation1.RData')
initial_states = initialisation(NUM_ACTORS,NUM_ASSESSORS)

output = vsp_hierarchy(initial_states,Y=simulation$Y,chain_name='test',n_itr=500000,n_record=100,noise_model='queue-jumping')
