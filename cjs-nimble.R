# Fit Cormack-Jolly-Seber model to the Dipper data 
# Bayesian approach 
# see Gimenez et al. (2007), Royle (2008)
# Nimble implementation

# packages
library(nimble)
library(basicMCMCplots)
library(coda)

# Read in the data: 
mydata <- as.matrix(read.table('dipper.txt'))
head(mydata)
dim(mydata)

# remove counts
mydata <- mydata[,-8]

# get nb of individuals and capture occasions
N <- dim(mydata)[1]
K <- dim(mydata)[2]

# Compute the date of first capture for each individual:
e <- NULL
for (i in 1:N){
temp <- 1:K
e <- c(e,min(temp[mydata[i,]>=1]))}

# Let's define the model. To do so, some notation first:

code <- nimbleCode({

# notation used
# STATES
# A for alive
# D for dead

# OBSERVATIONS
# 0 = non-observed (coded 1)
# 1 = observed (coded 2)

# prior on survival 
phi ~ dunif(0,1)

# prior on detection
p ~ dunif(0,1)

# probabilities for each initial state
px0[1] <- 1 # prob. of being in initial state A
px0[2] <- 0 # prob. of being in initial state D

# define probabilities of observations at t given states at t
po[1,1] <- 1-p
po[1,2] <- p
po[2,1] <- 1
po[2,2] <- 0

po.init[1,1] <- 0
po.init[1,2] <- 1
po.init[2,1] <- 1
po.init[2,2] <- 0

# define probabilities of states at t given states at t-1
px[1,1] <- phi
px[1,2] <- 1-phi
px[2,1] <- 0
px[2,2] <- 1

for (i in 1:nind){  # for each indiv
  # initial states
  z[i,first[i]] ~ dcat(px0[1:2])
  y[i,first[i]] ~ dcat(po.init[z[i,first[i]],1:2])
  for (j in (first[i]+1):nyear){  # loop over time

    # state equations
    z[i,j] ~ dcat(px[z[i,j-1],1:2])

	  # observation equations
	  y[i,j] ~ dcat(po[z[i,j],1:2])
    }
}
  
  
  })

# constant values 
constants <- list(nyear = K, nind = N, first = e)

# data
data <- list(y = mydata + 1)

# initial values
z_init <- mydata + 1
z_init[z_init == 2] <- 1
inits <- list(phi = 0.5, p = 0.5, z = z_init)

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## -1175.578

conf <- configureMCMC(Rmodel)

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, 2000, nburnin = 1000)

samplesPlot(samples, scale = TRUE, var = c("phi","p"))

summary(samples[,c("phi","p")])

# These results might be compared to the results obtained using E-SURGE 
# 0.560243, 0.9025836

