source(here::here("code", "score-functions-helpers.R"))
source(here::here("code", "SHM-helper-functions.R"))
source(here::here("code", "DHM-TV-general-data-generation.R"))
source(here::here("code", "PDA-2-pw-general.R"))

# Generate DHM with time-varying damping: ---------------------------------

grid_range <- c(0, 4 * pi)
n_grid <- 200
grid_points <- seq(from = grid_range[1],
                     to = grid_range[2],
                     length.out = n_grid)
n_i <- 500
num_iter <- 10
set.seed(1996)

damping_fun_test <- function(t) {
  - 0.01 *  (t - 2 * pi) ^2
}


dataset_i <- generate_DHM_TV_dataset_general(n = n_i,
                                             grid_points = grid_points,
                                             sigma = 0.5,
                                             mu_init = c(0, 1),
                                             intensity = 1,
                                             damping_fun = damping_fun_test,
                                             sigma_init = diag(rep(0.05, 2)))


# Translating it into Ramsay dialect. 

bbasis = create.bspline.basis(range(grid_points),nbasis = 51)

betaC = cbind(rep(-1,bbasis$nbasis),rep(0,bbasis$nbasis))
betafd = fd(betaC,bbasis)

Xfd = Data2fd(grid_points,dataset_i$x,bbasis)


tfine = seq(min(grid_points),max(grid_points),len=201)
nfine = length(tfine)
# Sig = (0.75 ^ 2) * dnorm(x = 2 * outer(tfine,tfine,'-'))


# We'll use this to try solving

# Note that this takes _forever_ ~90 minutes; I use a much
# smaller basis for beta below. 

#library(nleqslv)

#hmm = nleqslv(c(rep(-1,bbasis$nbasis),rep(0,bbasis$nbasis)),ScoreFn,
#           Xfd=Xfd,bbasis=bbasis,tfine=tfine,Sig=Sig)

#betaest = matrix(hmm$x,nbasis,2)

#plot(fd(betaest,bbasis))

# Now lets' check what happens with pda.fd

betaPar = fdPar(bbasis,int2Lfd(2),0)

pdaest = fRegress(deriv.fd(Xfd,2),list(Xfd,deriv.fd(Xfd,1)),list(betaPar,betaPar))

betaest.pda = cbind(pdaest$betaestlist[[1]]$fd$coef,pdaest$betaestlist[[2]]$fd$coef)

plot(fd(betaest.pda,bbasis))

matplot(dataset_i$x, type = "l")
# Now let's just try this with a few less basis functions

betabasis = create.bspline.basis(range(tfine),nbasis=20)
#betabasis = create.constant.basis(range(tfine))

#betabasis = bbasis

betaC = c(rep(-1,betabasis$nbasis),rep(0,betabasis$nbasis))

betafd= fd( matrix(betaC,betabasis$nbasis,2),betabasis)

#hmm.small = nleqslv(betaC,ScoreFn,Xfd=Xfd,bbasis=betabasis,tfine=tfine,Sig=Sig)

#betaest.s = matrix(hmm.small$x,betabasis$nbasis,2)


betaPar = fdPar(betabasis,int2Lfd(2),0)

pdaest.small = fRegress(deriv.fd(Xfd,2),list(Xfd,deriv.fd(Xfd,1)),list(betaPar,betaPar))

betaest.pda = cbind(pdaest.small$betaestlist[[1]]$fd$coef,pdaest.small$betaestlist[[2]]$fd$coef)


# X11()
plot(fd(betaest.pda,betabasis))
abline(h=-1)


# Let's try our own Gauss-Newton code

best = betaC


best = c(betaest.pda[,1],betaest.pda[,2])

for(i in 1:10){
  #	Sest = ScoreFndc(best,Xfd,betabasis,tfine,Sig)
  Sest = ScoreFndcSig(best,Xfd,betabasis,tfine)
  best = best - 0.5*solve(Sest$scoresdc,Sest$scores)
  print(c(i, sum(Sest$scores^2)) )
}

betaest.s = best
iGest = solve(Sest$scoresdc)


# While that runs we'll look at the error process

bvals = eval.basis(betabasis,tfine)
betafine = bvals%*%matrix(betaest.s,betabasis$nbasis,2,byrow=FALSE)


Xfine = eval.fd(Xfd,tfine)
DXfine = eval.fd(Xfd,tfine,1)
D2Xfine = eval.fd(Xfd,tfine,2)


errs = D2Xfine - betafine[,1]*Xfine - betafine[,2]*DXfine

# now we can look at the variance of the projections of Xeta 
# and DXeta onto the basis space

svals = rbind(t(bvals)%*%(Xfine*errs)*diff(betabasis$range)/nfine,
              t(bvals)%*%(DXfine*errs)*diff(betabasis$range)/nfine)

scov = cov(t(svals))		

# Should give us a variance

betavar = cbind(bvals,bvals)%*%iGest%*%scov%*%t(iGest)%*%t(cbind(bvals,bvals))*
  (nfine/diff(betabasis$range))^2/ncol(svals)

betaest.sd = matrix(sqrt(diag(betavar)),nfine,2,byrow=FALSE)

par(mar=c(5,5,1,1))
matplot(tfine,betafine,type='l',ylim=c(-1.5,0.5),lwd=2,lty=1,col=c('darkblue','darkred'),
        xlab='t',ylab=expression(beta(t)),cex.lab=1.5,cex.axis=1.5)
matplot(tfine,betafine+2*betaest.sd,type='l',lty=2,add=TRUE,col=c('darkblue','darkred'))
matplot(tfine,betafine-2*betaest.sd,type='l',lty=2,add=TRUE,col=c('darkblue','darkred'))
abline(h=c(-1),col=c('black'), lwd = 3)
lines(fd(betaest.pda,betabasis),lty=4,lwd=3,col=c('blue','red'))
lines(tfine, damping_fun_test(tfine), col = "black", lwd = 3)
legend('bottomleft',c(expression(beta[0](t)),expression(beta[1](t)),'Score Est','PDA', "PDA (bc)", "Truth(s)"),
       lty=c(1,1,1,4,3, 1),lwd=3,col=c('darkblue','darkred','darkblue','red', "red", "black"),cex=1.5,ncol=3)


x_i_test <- dataset_i$x
Dx_i_test <- dataset_i$dx
eps_i_test <- dataset_i$noise
D2x_i_test <- - x_i_test + damping_fun_test(grid_points) * Dx_i_test + eps_i_test
try(expr = {pda_test_bc <- do_pda_iterative_br_parallel(grid_points = grid_points,
                                                     x = x_i_test,
                                                     Dx = Dx_i_test,
                                                     D2x = D2x_i_test,
                                                     n = n_i, 
                                                     num_iter = num_iter,
                                                     verbose = TRUE, 
                                                     n_cores = 8)})


lines(grid_points, pda_test_bc$beta[, 2,11], col = "blue", lty = 3, lwd = 2)
lines(grid_points, pda_test_bc$beta[, 3,11], col = "red", lty = 3, lwd = 2)





