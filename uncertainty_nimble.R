# Fit multistate model with uncertainty in the state assignement to Sooty shearwater data (titis)
# Bayesian approach
# see Gimenez et al. (2012), Pradel (2005)

# packages
library(nimble)
library(basicMCMCplots)
library(coda)


# Read in the data: 
mydata <- read.table('titis2.txt')
head(mydata)
dim(mydata)
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
code <- nimbleCode({
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
  po[1:3,1:4] <- po1[1:3,1:3] %*% po2[1:3,1:4]
  po.init[1:3,1:4] <- po1.init[1:3,1:3] %*% po2[1:3,1:4]
  
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
  px[1:3,1:3] <- px1[1:3,1:3] %*% px2[1:3,1:3]
  
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
})

# Form the list of data and constants
constants <- list(N = N, 
                  Years = K,
                  First = e)
data = list(mydata=as.matrix(mydata+1))

# Generate inits for the latent states:
alive <- mydata
alive[alive==3] = 2
for (i in 1:N) {
  for (j in 1:K) {
    if (j > e[i] & mydata[i,j]==0) {alive[i,j] <- which(rmultinom(1, 1, c(1/2,1/2))==1)}
    if (j < e[i]) {alive[i,j] <- 0}
  }
}
alive1 <- as.matrix(alive)

# Now form the list of initial values:
inits <- list(
phiB = 0.5,
phiNB = 0.5,
psiNBB = 0.2,
psiBNB = 0.2,
piNB = 0.2,
pB = 0.2,
pNB = 0.5,
pB = 0.5, 
deltaNB = 0.1,
deltaB = 0.1,
alive = alive1)


Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## -317.776

conf <- configureMCMC(Rmodel, monitors = c("phiB","phiNB","psiNBB","psiBNB","piNB","pB","pNB","deltaNB","deltaB","alive"))

conf$printMonitors()
conf$printSamplers(byType = TRUE)
conf$printSamplers()

Rmcmc <- buildMCMC(conf, enableWAIC = TRUE)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

## Posterior inference

samples <- runMCMC(Cmcmc, 2500, nburnin = 1000, WAIC = TRUE)

samples$WAIC

samplesPlot(samples, scale = TRUE, var = c("phiB","phiNB","psiNBB","psiBNB","piNB","pB","pNB","deltaNB","deltaB"))

# Print results
summary(samples$samples[,c("phiB","phiNB","psiNBB","psiBNB","piNB","pB","pNB","deltaNB","deltaB")])


#      phiB            phiNB            psiNBB           psiBNB            piNB       
# Min.   :0.7611   Min.   :0.7653   Min.   :0.1505   Min.   :0.1238   Min.   :0.6027  
# 1st Qu.:0.8227   1st Qu.:0.8002   1st Qu.:0.2027   1st Qu.:0.1812   1st Qu.:0.6731  
# Median :0.8369   Median :0.8134   Median :0.2233   Median :0.2054   Median :0.6895  
# Mean   :0.8358   Mean   :0.8121   Mean   :0.2234   Mean   :0.2065   Mean   :0.6884  
# 3rd Qu.:0.8503   3rd Qu.:0.8226   3rd Qu.:0.2409   3rd Qu.:0.2311   3rd Qu.:0.7061  
# Max.   :0.8937   Max.   :0.8554   Max.   :0.3310   Max.   :0.3369   Max.   :0.7541  
#       pB              pNB            deltaNB           deltaB      
# Min.   :0.5228   Min.   :0.4848   Min.   :0.1493   Min.   :0.5889  
# 1st Qu.:0.5715   1st Qu.:0.5493   1st Qu.:0.1853   1st Qu.:0.6801  
# Median :0.5921   Median :0.5675   Median :0.1933   Median :0.7044  
# Mean   :0.5943   Mean   :0.5668   Mean   :0.1941   Mean   :0.7065  
# 3rd Qu.:0.6151   3rd Qu.:0.5825   3rd Qu.:0.2019   3rd Qu.:0.7309  
# Max.   :0.6878   Max.   :0.6512   Max.   :0.2356   Max.   :0.8209  

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

