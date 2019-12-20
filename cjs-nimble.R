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

# OBSERVATIONS (+1)
# 0 = non-detected
# 1 = detected

# STATES
# 1 = alive
# 2 = dead

# PARAMETERS
# phi  survival
# p detection

# input code
code <- nimbleCode({
phi ~ dunif(0,1)
p ~ dunif(0,1)
for(i in 1:nind) {
	x[i, first[i]] <- 1
	y[i, first[i]] <- 1
	for(t in (first[i]+1):k) {
		x[i,t] ~ dbern(phi * x[i,t-1])
		y[i,t] ~ dbern(p * x[i,t])
							} # t loop
				} # i loop
})

# constant values 
constants <- list(k = K, nind = N, first = e)

# data
data <- list(y = mydata)

# initial values
x_init <- mydata
x_init[x_init ==0] <- 1
inits <- list(phi = 0.5, p = 0.5, x = x_init)

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## -317.776

conf <- configureMCMC(Rmodel)

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, 2000, nburnin = 1000)

samplesPlot(samples, scale = TRUE, var = c("phi","p"))

summary(samples[,c("phi","p")])

# These results might be compared to the results obtained using E-SURGE 
# 0.560243, 0.9025836

