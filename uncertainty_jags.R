# Fit multistate model with uncertainty in the state assignement to Sooty shearwater data (titis)
# Bayesian approach
# see Gimenez et al. (2012), Pradel (2005)

# Read in the data: 
mydata = read.table('titis2.txt')
head(mydata)
dim(mydata)
N = dim(mydata)[1]
K = dim(mydata)[2]

# Compute the date of first capture for each individual:
e <- NULL
for (i in 1:N){
temp <- 1:K
e <- c(e,min(temp[mydata[i,]>=1]))}

# Let's define the model. To do so, some notation first:

# OBSERVATIONS (+1)
# 0 = non-detected
# 1 = seen and ascertained as non-breeder
# 2 = seen and ascertained as breeder
# 3 = not ascertained

# STATES
# 1 = alive non-breeder
# 2 = alive breeder
# 3 = dead

# PARAMETERS
# phiNB  survival prob. of non-breeders
# phiB  survival prob. of breeders
# pNB  detection prob. of non-breeders
# pB  detection prob. of breeders
# psiNBB transition prob. from non-breeder to breeder
# psiBNB transition prob. from breeder to non-breeder
# piNB prob. of being in initial state non-breeder
# deltaNB prob to ascertain the breeding status of an individual encountered as non-breeder
# deltaB prob to ascertain the breeding status of an individual encountered as breeder

# Now the model:
M <- function() {
  # DEFINE PARAMETERS	
  # probabilities for each INITIAL STATES
  px0[1] <- piNB # prob. of being in initial state NB
  px0[2] <- 1 - piNB # prob. of being in initial state B
  px0[3] <- 0 # prob. of being in initial state dead
  
  # OBSERVATION PROCESS: probabilities of observations (columns) at a given occasion given states (rows) at this occasion
  # step 1: detection
  po1[1,1] <- 1 - pNB
  po1[1,2] <- pNB
  po1[1,3] <- 0
  po1[2,1] <- 1 - pB
  po1[2,2] <- 0
  po1[2,3] <- pB
  po1[3,1] <- 1
  po1[3,2] <- 0
  po1[3,3] <- 0
  
  po1.init[1,1] <- 0
  po1.init[1,2] <- 1
  po1.init[1,3] <- 0
  po1.init[2,1] <- 0
  po1.init[2,2] <- 0
  po1.init[2,3] <- 1
  po1.init[3,1] <- 1
  po1.init[3,2] <- 0
  po1.init[3,3] <- 0
  # step 2: assignement
  po2[1,1] <- 1
  po2[1,2] <- 0
  po2[1,3] <- 0
  po2[1,4] <- 0
  po2[2,1] <- 0
  po2[2,2] <- deltaNB
  po2[2,3] <- 0
  po2[2,4] <- 1 - deltaNB
  po2[3,1] <- 0
  po2[3,2] <- 0
  po2[3,3] <- deltaB
  po2[3,4] <- 1 - deltaB
  # form the matrix product
  po <- po1 %*% po2
  po.init <- po1.init %*% po2
  
  # STATE PROCESS: probabilities of states at t+1 (columns) given states at t (rows)
  # step 1: survival
  px1[1,1] <- phiNB
  px1[1,2] <- 0
  px1[1,3] <- 1 - phiNB
  px1[2,1] <- 0
  px1[2,2] <- phiB
  px1[2,3] <- 1 - phiB
  px1[3,1] <- 0
  px1[3,2] <- 0
  px1[3,3] <- 1
  # step 2: transition
  px2[1,1] <- 1 - psiNBB
  px2[1,2] <- psiNBB
  px2[1,3] <- 0
  px2[2,1] <- psiBNB
  px2[2,2] <- 1 - psiBNB
  px2[2,3] <- 0
  px2[3,1] <- 0
  px2[3,2] <- 0
  px2[3,3] <- 1
  # form the matrix product
  px <- px1 %*% px2
  
  for (i in 1:N)  # for each ind
  {
    
    # estimated probabilities of initial states are the proportions in each state at first capture occasion
    alive[i,First[i]] ~ dcat(px0[1:3])
    mydata[i,First[i]] ~ dcat(po.init[alive[i,First[i]],1:4])
    
    for (j in (First[i]+1):Years)  # loop over time
      
    {
      
      ## STATE EQUATIONS ##
      # draw states at j given states at j-1
      alive[i,j] ~ dcat(px[alive[i,j-1],1:3])
      
      ## OBSERVATION EQUATIONS ##
      # draw observations at j given states at j
      mydata[i,j] ~ dcat(po[alive[i,j],1:4])
      
    }
    
  }
  
  # PRIORS 
  phiNB ~ dunif(0, 1)
  phiB ~ dunif(0, 1)
  pNB ~ dunif(0, 1)
  pB ~ dunif(0, 1)
  psiNBB ~ dunif(0, 1)
  psiBNB ~ dunif(0, 1)
  piNB ~ dunif(0, 1)
  deltaNB ~ dunif(0, 1)
  deltaB ~ dunif(0, 1)
}

# Form the list of data
mydatax = list(N=N,Years=K,mydata=as.matrix(mydata+1),First=e)

# Generate inits for the latent states:
alive = mydata
# assign 2s to 3s
alive [alive==3] = 2
for (i in 1:N) {
  for (j in 1:K) {
    if (j > e[i] & mydata[i,j]==0) {alive[i,j] = which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j < e[i]) {alive[i,j] = NA}
  }
}
alive1 <- as.matrix(alive)
alive = mydata
# assign 1s to 3s
alive [alive==3] = 1
for (i in 1:N) {
  for (j in 1:K) {
    if (j > e[i] & mydata[i,j]==0) {alive[i,j] = which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j < e[i]) {alive[i,j] = NA}
  }
}
alive2 <- as.matrix(alive)

# Now form the list of initial values:
init1 <- list(pB=0.5,phiNB=0.3,alive=alive1)
# second list of inits
init2 <- list(pB=0.5,phiNB=0.6,alive=alive2)
# concatenate list of initial values
inits <- list(init1,init2)

# Specify the parameters to be monitored
parameters <- c("phiB","phiNB","psiNBB","psiBNB","piNB","pB","pNB","deltaNB","deltaB")

# Tadaaaaaaan, fit the model:
library(R2jags)
out <- jags(mydatax,inits,parameters, M,n.chains=2,n.iter=2000,n.burnin=500)

# Check convergence:
traceplot(out,ask=F)

# Nice plots:
library(lattice)
jagsfit.mcmc <- as.mcmc(out)
densityplot(jagsfit.mcmc)

# Print results
print(out)
# These results are to be compared to the results obtained using E-SURGE 
# (Table 1 in Gimenez et al. 2012):
# deltaB | 0.737576452 0.616347581 0.831002071 0.055236538 
# deltaNB | 0.187797018 0.162175744 0.216420470 0.013831853 
# pB | 0.597718828 0.532890645 0.659302061 0.032413822 
# pNB | 0.564643587 0.508403711 0.619267746 0.028396321 
# phiB | 0.837489708 0.796415345 0.871612997 0.019139418 
# phiNB | 0.814037056 0.779702384 0.844090508 0.016414429 
# piNB | 0.704217686 0.646248176 0.756270393 0.028149150 
# psiBNB | 0.226471935 0.144866984 0.335984429 0.048899140 
# psiNBB | 0.219402142 0.173907703 0.272866821 0.025255272 
