# Fit Cormack-Jolly-Seber model to the Dipper data 
# Bayesian approach
# see Gimenez et al. (2007), Royle (2008)
# Jags implementation

# Read in the data: 
mydata <- read.table('dipper.txt')
head(mydata)
dim(mydata)
# remove counts
mydata <- mydata[,-8]
N <- dim(mydata)[1]
K <- dim(mydata)[2]

# Compute the date of first capture for each individual:
e <- NULL
for (i in 1:N){
temp <- 1:K
e <- c(e,min(temp[mydata[i,]>=1]))}

# Let's define the model. To do so, some notation first:

# OBSERVATIONS (+1)
# 0 = non-detected
# 1 = detected

# STATES
# 1 = alive
# 2 = dead

# PARAMETERS
# phi  survival
# p detection

# Now the model:
model <- function() {
  # DEFINE PARAMETERS	
  # probabilities for each INITIAL STATES
  px0[1] <- 1 # prob. of being in initial state alive
  px0[2] <- 0 # prob. of being in initial state dead
  
  # OBSERVATION PROCESS: probabilities of observations (columns) at a given occasion given states (rows) at this occasion
  # step 1: detection
  po[1,1] <- 1 - p
  po[1,2] <- p
  po[2,1] <- 1
  po[2,2] <- 0
  po.init[1,1] <- 0
  po.init[1,2] <- 1
  po.init[2,1] <- 1
  po.init[2,2] <- 0
  
  # STATE PROCESS: probabilities of states at t+1 (columns) given states at t (rows)
  # step 1: survival
  px[1,1] <- phi
  px[1,2] <- 1 - phi
  px[2,1] <- 0
  px[2,2] <- 1
  
  for (i in 1:N)  # for each indiv
  {
    
    # estimated probabilities of initial states are the proportions in each state at first capture occasion
    alive[i,First[i]] ~ dcat(px0[1:2])
    mydata[i,First[i]] ~ dcat(po.init[alive[i,First[i]],1:2])
    
    for (j in (First[i]+1):Years)  # loop over time
      
    {
      
      ## STATE EQUATIONS ##
      # draw states at j given states at j-1
      alive[i,j] ~ dcat(px[alive[i,j-1],1:2])
      
      ## OBSERVATION EQUATIONS ##
      # draw observations at j given states at j
      mydata[i,j] ~ dcat(po[alive[i,j],1:2])
      
    }
    
  }
  
  # PRIORS 
  phi ~ dunif(0, 1)
  p ~ dunif(0, 1)
}

# Form the list of data
mydatax <- list(N=N,Years=K,mydata=as.matrix(mydata+1),First=e)

# Generate inits for the latent states
x.init <- mydata
for (i in 1:N){
	if (e[i] == 1) next
	if (e[i] > 1) x.init[i,1:(e[i]-1)] <- NA
}
x.init[x.init==0] <- 1
z <- as.matrix(x.init)

# Now form the list of initial values:
init1 <- list(p=0.1,phi=0.1,alive=z)
# second list of inits
init2 <- list(p=0.8,phi=0.8,alive=z)
# concatenate list of initial values
inits <- list(init1,init2)

# Specify the parameters to be monitored
parameters <- c("phi","p")

# Tadaaaaaaan, fit the model:
library(R2jags)
out <- jags(mydatax,inits,parameters, model,n.chains=2,n.iter=2000,n.burnin=500)

# Check convergence:
traceplot(out,ask=T)

# Nice plots:
library(lattice)
jagsfit.mcmc <- as.mcmc(out)
densityplot(jagsfit.mcmc)

# Print results
print(out)
# These results might be compared to the results obtained using E-SURGE 
# 0.560243, 0.9025836
