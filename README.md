# Fitting multievent models in R, Nimble and JAGS
### by Olivier Gimenez
### August 12, 2016 and December 19, 2019.

## What it does

Here I provide R codes to fit multievent capture-recapture models (Pradel et al. 2005). Multievent models are hidden Markov models that are helpful in lots of situations to analyse capture-recapture data (see this [list of applications](https://multievent.sciencesconf.org/resource/page/id/9) for example). 

I show how to obtain maximum-likelihood estimates using R and Bayesian estimates using Nimble and JAGS. Two examples are considered. First a simple Cormack-Jolly-Seber model is illustrated with the classical Dipper dataset (Pradel 2005; Gimenez et al. 2007). Second, a multistate model with uncertainty in the state assignement is illustrated with a dataset on Sooty shearwaters (Pradel 2005; Gimenez et al. 2012).

## What it contains

### Cormack-Jolly-Seber example using the Dipper dataset
* `cjs_nimble.R`: Bayesian fitting using R and Nimble
* `cjs_jags.R`: Bayesian fitting using R and JAGS
* `cjs_R.R`: maximum-likelihood fitting using R
* `dipper.txt`: the Dipper dataset

### Multistate with uncertain state example using the Sooty shearwater dataset
* `uncertainty_nimble.R`: Bayesian fitting using R and Nimble
* `uncertainty_jags.R`: Bayesian fitting using R and JAGS
* `uncertainty_R.R`: maximum-likelihood fitting using R
* `titis2.txt`: the Sooty shearwater dataset

## References

Gimenez, O., Lebreton, J.-D., Gaillard, J.-M., Choquet, R. and R. Pradel (2012). [Estimating demographic parameters using hidden process dynamic models](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Gimenezetal2012TPB.pdf). Theoretical Population Biology 82: 307-316.

Gimenez, O., V. Rossi, R. Choquet, C. Dehais, B. Doris, H. Varella, J.-P. Vila and R. Pradel (2007). [State-space modelling of data on marked individuals](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Gimenezetal2007EcologicalModelling.pdf). Ecological Modelling 206: 431-438.

Pradel,R. (2005). [Multievent: an extension of multistate capture–recapture models to uncertain states](http://200.46.218.171/bds-cbc/sites/default/files/Pradel%20Biometrics%202005.pdf). Biometrics 61: 442–447.
