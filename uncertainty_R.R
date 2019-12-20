# Fit multistate model with uncertainty in the state assignement to Sooty shearwater data (titis)
# Maximum-likelihood approach
# see Pradel (2005), Gimenez et al. (2012)

# -log(lik) 
devMULTIEVENT <- function(b,data,eff,e,garb,nh,km1){
	
# data encounter histories, eff counts
# e vector of dates of first captures
# garb vector of initial states 
# km1 nb of recapture occasions (nb of capture occ - 1)
# nh nb ind

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

# logit link for all parameters
# note: below, we decompose the state and obs process in two steps composed of binomial events, 
# which makes the use of the logit link appealing; 
# if not, a multinomial (aka generalised) logit link should be used
piNB <- 1/(1+exp(-b[1]))
phiNB <- 1/(1+exp(-b[2]))
phiB <- 1/(1+exp(-b[3]))
psiNBB <- 1/(1+exp(-b[4]))
psiBNB <- 1/(1+exp(-b[5]))
pNB <- 1/(1+exp(-b[6]))
pB <- 1/(1+exp(-b[7]))
deltaNB <- 1/(1+exp(-b[8]))
deltaB <- 1/(1+exp(-b[9]))

# prob of obs (rows) cond on states (col)
B1 <- matrix(c(1-pNB,pNB,0,1-pB,0,pB,1,0,0),nrow=3,ncol=3,byrow=T)
B2 <- matrix(c(1,0,0,0,0,deltaNB,0,1-deltaNB,0,0,deltaB,1-deltaB),nrow=3,ncol=4,byrow=T)
B <- t(B1 %*% B2)
#B = t(matrix(c(1-pNB,pNB*deltaNB,0,pNB*(1-deltaNB),1-pB,0,pB*deltaB,pB*(1-deltaB),1,0,0,0),nrow=3,ncol=4,byrow=T))

# first encounter
BE1 <- matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3,byrow=T)BE2 <- matrix(c(1,0,0,0,0,deltaNB,0,1-deltaNB,0,0,deltaB,1-deltaB),nrow=3,ncol=4,byrow=T)
BE <- t(BE1 %*% BE2) 
#BE = t(matrix(c(0,deltaNB,0,(1-deltaNB),0,0,deltaB,(1-deltaB),1,0,0,0),nrow=3,ncol=4,byrow=T))
# prob of states at t+1 given states at t
A1 <- matrix(c(phiNB,0,1-phiNB,0,phiB,1-phiB,0,0,1),nrow=3,ncol=3,byrow=T)
A2 <- matrix(c(1-psiNBB,psiNBB,0,psiBNB,1-psiBNB,0,0,0,1),nrow=3,ncol=3,byrow=T)
A <- A1 %*% A2
#A <- matrix(c(phiNB*(1-psiNBB),phiNB*psiNBB,1-phiNB,phiB*psiBNB,phiB*(1-psiBNB),1-phiB,0,0,1),nrow=3,ncol=3,byrow=T)

# init states
PI <- c(piNB,1-piNB,0)
# likelihood
   l <- 0   for (i in 1:nh) # loop on ind
   {      ei <- e[i] # date of first det      oe <- garb[i] + 1 # init obs      evennt <- data[,i] + 1 # add 1 to obs to avoid 0s in indexing      ALPHA <- PI*BE[oe,]      for (j in (ei+1):(km1+1)) # cond on first capture
      {
      	if ((ei+1)>(km1+1)) {break} # sous MATLAB la commande >> 8:7 rend >> null, alors que sous R, ça rend le vecteur c(8,7)!        ALPHA <- (ALPHA %*% A)*B[evennt[j],]      }      l <- l + logprot(sum(ALPHA))*eff[i]   }    l <- -l
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
data <- read.table('titis2.txt')

# define various quantities
nh <- dim(data)[1]
k <- dim(data)[2]
km1 <- k-1

# counts
eff <- rep(1,nh)

# compute the date of first capture fc, and state at initial capture init.state
fc <- NULL
init.state <- NULL
for (i in 1:nh){
temp <- 1:k
fc <- c(fc,min(which(data[i,]!=0)))
init.state <- c(init.state,data[i,fc[i]])
}

# init valuesbinit <- runif(9)

# transpose data
data <- t(data)

# fit model
deb <- Sys.time()
tmpmin <- optim(binit,devMULTIEVENT,NULL,hessian=FALSE,data,eff,fc,init.state,nh,km1,method="BFGS",control=list(trace=1, REPORT=1))
fin <- Sys.time()
fin-deb 

# get estimates and back-transform
x <- tmpmin$parpiNB <- 1/(1+exp(-x[1]))
phiNB <- 1/(1+exp(-x[2]))
phiB <- 1/(1+exp(-x[3]))
psiNBB <- 1/(1+exp(-x[4]))
psiBNB <- 1/(1+exp(-x[5]))
pNB <- 1/(1+exp(-x[6]))
pB <- 1/(1+exp(-x[7]))
deltaNB <- 1/(1+exp(-x[8]))
deltaB <- 1/(1+exp(-x[9]))

piNB 
phiNB 
phiB 
psiBNB 
psiNBB 
pNB
pB 
deltaNB 
deltaB 

# These estimates should be compared to those we get in E-SURGE:
# Par#     1#    IS(  1,  1)(  1,  1)(  1   1) | 0.704217686 0.646248176 0.756270393 0.028149150 
# Par#    15#    S(  1,  1)(  1,  1)(  1   1) | 0.814037056 0.779702384 0.844090508 0.016414429 
# Par#    16#    S(  2,  2)(  1,  1)(  1   1) | 0.837489708 0.796415345 0.871612997 0.019139418 
# Par#    46#    T(  2,  1)(  1,  1)(  1   2) | 0.226471935 0.144866984 0.335984429 0.048899140 
# Par#    47#    T(  1,  2)(  1,  1)(  1   2) | 0.219402142 0.173907703 0.272866821 0.025255272 
# Par#   113#    E(  1,  2)(  2,  2)(  1   1) | 0.564643587 0.508403711 0.619267746 0.028396321 
# Par#   114#    E(  2,  3)(  2,  2)(  1   1) | 0.597718828 0.532890645 0.659302061 0.032413822 
# Par#   141#    D(  2,  2)(  1,  1)(  1   2) | 0.187797018 0.162175744 0.216420470 0.013831853 
# Par#   142#    D(  3,  3)(  1,  1)(  1   2) | 0.737576452 0.616347581 0.831002071 0.055236538 
