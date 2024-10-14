source("~/source-code/pofun.R")
source("~/source-code/pofun_t.R")
source("~/source-code/vspfun.R")
source("~/vsp_hierarchy_func.R")
source("~/vsp_hierarchy_mcmc.R")

dir='~/output'
setwd(dir)

clustering = TRUE

# Hyperparameters

RHO_HYPERPAMA = 1/3
THETA_HYPERPAMA = 5    # MALLOWS
R_HYPERPAMA = 1.5  # QUEUE-JUMPING
PDP_THETA = .5 # PDP_THETA > -PDP_ALPHA
PDP_ALPHA = 0.15
RHO_PROPOSAL = 1/2
TAU_PROPOSAL = 1/2

NUM_ACTORS = 10
NUM_ASSESSORS = 3
NUM_ORDERS = rep(5,NUM_ASSESSORS)
NUM_LISTS = sum(NUM_ORDERS)


# Simulation

# po0 = transitive.reduction(randDGDAG(n=NUM_ACTORS))
po0=matrix(c(0,1,0,1,rep(0,3),1,rep(0,6),1,rep(0,4),1,rep(0,5)),nrow=5,byrow=TRUE)
rownames(po0)=colnames(po0)=1:NUM_ACTORS
showDAG(po0)
U0=PO2Randlatent(transitive.closure(po0,mat=TRUE,loops=FALSE),K=2)
showDAG(transitive.reduction(u2po(U0)))
# po2tree(transitive.closure(po0,mat=TRUE,loops=FALSE))
tau=0
rho=0.9
Sigma=matrix(rho,nrow=2,ncol=2); diag(Sigma)=1
U_error = lapply(rep(NUM_ACTORS,NUM_ASSESSORS),function(n) rmvnorm(n,sigma=(1-tau^2)*Sigma))
U = lapply(U_error, function(ue) tau*U0+ue)
while (!all(sapply(lapply(U,u2po),function(po){po2tree(po)$is_vsp}))){
  U0=PO2Randlatent(transitive.closure(po0,mat=TRUE,loops=FALSE),K=2)
  U = lapply(U_error, function(ue) tau*U0+ue)
}

tcs = lapply(U, u2po)
trs = lapply(tcs, transitive.reduction)
showDAG(trs[[1]])
showDAG(trs[[2]])
showDAG(trs[[3]])
showDAG(trs[[4]])
showDAG(trs[[5]])

simulation = data_simulation(NUM_ACTORS, NUM_ASSESSORS, NUM_ORDERS, U0=U0, U=U)
Y = simulation$Y
save(simulation,file='simulation1.RData')

load('simulation1.RData')
initial_states = initialisation(NUM_ACTORS,Y=simulation$Y,clustering=TRUE)

output = vsp_hierarchy(initial_states,chain_name='test',n_itr=1e5,n_record=10,
                       noise_model='queue-jumping',hierarchy=TRUE,clustering=clustering)
