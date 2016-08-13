# Fit Cormack-Jolly-Seber model to the Dipper data 
# Maximum-lik approach
# see Pradel (2005)

# -log(lik) 
devCJS <- function(b,data,eff,e,garb,nh,km1){
# data encounter histories, eff counts# e vector of dates of first captures# garb vector of initial states # km1 nb of recapture occasions (nb of capture occ - 1)
# nh nb ind

# logit linkphi <- 1/(1+exp(-b[1]))p <- 1/(1+exp(-b[2]))

# prob of obs (rows) cond on states (col)
B <- array(0,c(2,2,km1)) # allows for time-dep in the parameters
for (t in 1:km1){
	B[,,t] <- matrix(c(1-p,1,p,0),nrow=2,ncol=2,byrow=T)
}

# first encounter
BE <- matrix(c(0,1,1,0),nrow=2,ncol=2,byrow=T)    # prob of states at t+1 given states at tA <- array(0,c(2,2,km1))
	for (t in 1:km1){
		A[,,t] <- matrix(c(phi, 1 - phi,0,1),nrow=2,ncol=2,byrow=T)
}

# init statesPI <- array(0,c(1,2,km1+1))
for (t in 1:(km1+1)){
		PI[,,t] <- c(1,0)
}

# likelihood
   l <- 0   for (i in 1:nh) # loop on ind
   {      ei <- e[i] # date of first det      oe <- garb[i] + 1 # init obs      evennt <- data[,i] + 1 # non-det become 1s, det 2s      ALPHA <- PI[,,ei]*BE[oe,]      for (j in (ei+1):(km1+1)) # cond on first capture
      {
      	if ((ei+1)>(km1+1)) {break} # sous MATLAB la commande >> 8:7 rend >> null, alors que sous R, ça rend le vecteur c(8,7)!        ALPHA <- (ALPHA %*% A[,,j-1])*B[evennt[j],,j-1]      }      l <- l + logprot(sum(ALPHA))*eff[i]   }    l <- -l
    l
}

# avoid explosion of log(v) for small values of v
logprot <- function(v){
eps <- 2.2204e-016
u <- log(eps) * (1+vector(length=length(v)))
index <- (v>eps)
u[index] <- log(v[index])
u
}

# read in data
data <- read.table('dipper.txt')

# various quantities to be defined
nkp1 <- dim(data)[2]
nh <- dim(data)[1]
k <- nkp1-1
km1 <- k-1

# counts
eff <- data[,nkp1]

# select ecounter histories only - remove counts
data <- data[,1:k]

# compute the date of first capture fc, and state at initial capture init.state
fc <- NULL
init.state <- NULL
for (i in 1:nh){
temp <- 1:k
fc <- c(fc,min(temp[data[i,]==1]))
init.state <- c(init.state,data[i,fc[i]])
}

# init valuesbinit <- rep(0,2)

# transpose data
data <- t(data)

# fit model
deb=Sys.time()
tmpmin <- optim(binit,devCJS,NULL,hessian=FALSE,data,eff,fc,init.state,nh,km1,method="BFGS",control=list(trace=1, REPORT=1))
fin=Sys.time()
fin-deb

# get estimates and back-transform
x <- tmpmin$parphi <- exp(x[1])/(1+exp(x[1]))p <- exp(x[2])/(1+exp(x[2]))
phi
p
