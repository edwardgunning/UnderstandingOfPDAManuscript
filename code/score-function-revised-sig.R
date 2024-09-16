library(fda)
library(deSolve)


Ord2ODE = function(t,x,betafd){
	if(t < betafd$basis$range[1] ){ t = betafd$basis$range[1] }
	if(t > betafd$basis$range[2] ){ t = betafd$basis$range[2] }
	bvals = as.vector(eval.fd(betafd,t))
	dx = c(x[2], bvals[1]*x[1] + bvals[2]*x[2])
	return(list(dx))
}


# Here is the ODE that also calculates the 
# derivative with respect to coefficients

Ord2ODEdc = function(t,x,betafd){
	if(t < betafd$basis$range[1] ){ t = betafd$basis$range[1] }
	if(t > betafd$basis$range[2] ){ t = betafd$basis$range[2] }
	
	bvals = eval.basis(betafd$basis,t)
	betavals = bvals%*%betafd$coef
	
	J = matrix(c(0,1,betavals[1],betavals[2]),2,2,byrow=TRUE)
	bvals = diag(2)%x%t(bvals)%x%matrix(c(0,1),2,1)
	 
	bigJ = rbind(cbind(J,matrix(0,2,betafd$basis$nbasis*4)),
				 cbind(bvals, diag(2*betafd$basis$nbasis)%x%J) )
				 
	dx = bigJ%*%x
	return(list(dx))
}


Ord2Trans = function(tfine,betafd){

	nfine = length(tfine)
	
	Phi = array(0,c(nfine,nfine,2,2))
	
	for(i in 1:(nfine-1)){
		val = lsoda(c(1,0),tfine[i:nfine],Ord2ODE,betafd)
		Phi[i:nfine,i,,1] = val[,2:3]
		val = lsoda(c(0,1),tfine[i:nfine],Ord2ODE,betafd)
		Phi[i:nfine,i,,2] = val[,2:3]		
	}
	Phi[nfine,nfine,,1] = c(1,0)
	Phi[nfine,nfine,,2] = c(0,1)
	
	return(Phi)
}

# Here we can calculate both the transition matrix
# and its derivative with respect to coefficients

Ord2Transdc = function(tfine,betafd){

	nfine = length(tfine)
	nbasis = betafd$basis$nbasis
	
	Phi = array(0,c(nfine,nfine,2,2))
	dPhi = array(0,c(nfine,nfine,2,2,2*nbasis))	
		
	for(i in 1:(nfine-1)){
		val = lsoda(c(1,0,rep(0,4*nbasis)),tfine[i:nfine],Ord2ODEdc,betafd)
		Phi[i:nfine,i,,1] = val[,2:3]
		dPhi[i:nfine,i,,1,] = matrix(val[,3+(1:(4*nbasis))],2,2*nbasis,byrow=TRUE)
		val = lsoda(c(0,1,rep(0,4*nbasis)),tfine[i:nfine],Ord2ODEdc,betafd)
		Phi[i:nfine,i,,2] = val[,2:3]	
		dPhi[i:nfine,i,,2,] = matrix(val[,3+(1:(4*nbasis))],2,2*nbasis,byrow=TRUE)
	}
	Phi[nfine,nfine,,1] = c(1,0)
	Phi[nfine,nfine,,2] = c(0,1)
	
	return(list(Phi=Phi,dPhi=dPhi))
}


ScoreFn = function(beta,Xfd,bbasis,tfine,Sig){

	beta = matrix(beta,bbasis$nbasis,2,byrow=FALSE)
	totm = diff(bbasis$range)  #total area to integrate

	nfine = length(tfine)

	Xfine = eval.fd(Xfd,tfine)
	DXfine = eval.fd(Xfd,tfine,1)
	D2Xfine = eval.fd(Xfd,tfine,2)

	n = ncol(Xfine)
	
	bvals = eval.basis(bbasis,tfine)

	betafd = fd(beta,bbasis)
	bfine = bvals%*%beta
	
	errs = D2Xfine - bfine[,2]*DXfine - bfine[,1]*Xfine
	Xerrs = apply(Xfine*errs,1,mean)
	DXerrs = apply(DXfine*errs,1,mean)


	Phi = Ord2Trans(tfine,betafd)
	sBias1 = (Phi[,,1,2]*Sig)%*%rep(totm/nfine,nfine)
	sBias2 = (Phi[,,2,2]*Sig)%*%rep(totm/nfine,nfine)
	
	
	return( c(t(bvals)%*%(Xerrs - sBias1), t(bvals)%*%(DXerrs - sBias2) ))
	
}



ScoreFndc = function(beta,Xfd,bbasis,tfine,Sig){

	beta = matrix(beta,bbasis$nbasis,2,byrow=FALSE)
	totm = diff(bbasis$range)  #total area to integrate

	nfine = length(tfine)

	Xfine = eval.fd(Xfd,tfine)
	DXfine = eval.fd(Xfd,tfine,1)
	D2Xfine = eval.fd(Xfd,tfine,2)

	n = ncol(Xfine)
	
	bvals = eval.basis(bbasis,tfine)

	betafd = fd(beta,bbasis)
	bfine = bvals%*%beta
	
	errs = D2Xfine - bfine[,2]*DXfine - bfine[,1]*Xfine
	Xerrs = apply(Xfine*errs,1,mean)
	DXerrs = apply(DXfine*errs,1,mean)

	Phi = Ord2Transdc(tfine,betafd)
	sBias1 = (Phi$Phi[,,1,2]*Sig)%*%rep(totm/nfine,nfine)
	sBias2 = (Phi$Phi[,,2,2]*Sig)%*%rep(totm/nfine,nfine)
		
	scores =  c(t(bvals)%*%(Xerrs - sBias1), t(bvals)%*%(DXerrs - sBias2) )
	
	## Now let's work on the derivative of the scores
	
	regdc = rbind(cbind(t(bvals)%*%diag(apply(Xfine^2,1,mean))%*%bvals,
				    t(bvals)%*%diag(apply(Xfine*DXfine,1,mean))%*%bvals),
				  cbind(t(bvals)%*%diag(apply(Xfine*DXfine,1,mean))%*%bvals,
				    t(bvals)%*%diag(apply(DXfine^2,1,mean))%*%bvals))
	
	sBias1dc = matrix(NA,bbasis$nbasis,2*bbasis$nbasis)
	sBias2dc = matrix(NA,bbasis$nbasis,2*bbasis$nbasis)
	
	for(i in 1:(2*bbasis$nbasis)){
		sBias1dc[,i] = t(bvals)%*%(Phi$dPhi[,,1,2,i]*Sig)%*%rep(totm/nfine,nfine)
		sBias2dc[,i] = t(bvals)%*%(Phi$dPhi[,,2,2,i]*Sig)%*%rep(totm/nfine,nfine)
	}
	
	scoresdc = -regdc  - rbind(sBias1dc,sBias2dc)
	
	# Return both of these
	
	return(list(scores = scores,scoresdc=scoresdc))
}



# Now including the estmation of Sigma in the score

ScoreFndcSig = function(beta,Xfd,bbasis,tfine){

	beta = matrix(beta,bbasis$nbasis,2,byrow=FALSE)
	totm = diff(bbasis$range)  #total area to integrate

	nfine = length(tfine)

	Xfine = eval.fd(Xfd,tfine)
	DXfine = eval.fd(Xfd,tfine,1)
	D2Xfine = eval.fd(Xfd,tfine,2)

	n = ncol(Xfine)
	
	bvals = eval.basis(bbasis,tfine)

	betafd = fd(beta,bbasis)
	bfine = bvals%*%beta
	
	errs = D2Xfine - bfine[,2]*DXfine - bfine[,1]*Xfine
	Xerrs = apply(Xfine*errs,1,mean)
	DXerrs = apply(DXfine*errs,1,mean)

	Sig = (1/n)*errs%*%t(errs)


	Phi = Ord2Transdc(tfine,betafd)
	sBias1 = (Phi$Phi[,,1,2]*Sig)%*%rep(totm/nfine,nfine)
	sBias2 = (Phi$Phi[,,2,2]*Sig)%*%rep(totm/nfine,nfine)
		
	scores =  c(t(bvals)%*%(Xerrs - sBias1), t(bvals)%*%(DXerrs - sBias2) )
	
	## Now let's work on the derivative of the scores
	
	regdc = rbind(cbind(t(bvals)%*%diag(apply(Xfine^2,1,mean))%*%bvals,
				    t(bvals)%*%diag(apply(Xfine*DXfine,1,mean))%*%bvals),
				  cbind(t(bvals)%*%diag(apply(Xfine*DXfine,1,mean))%*%bvals,
				    t(bvals)%*%diag(apply(DXfine^2,1,mean))%*%bvals))
	
	sBias1dc = matrix(NA,bbasis$nbasis,2*bbasis$nbasis)
	sBias2dc = matrix(NA,bbasis$nbasis,2*bbasis$nbasis)
	
	for(i in 1:2){
     if(i == 1){ XX = Xfine }else{ XX = DXfine }
	 for(j in 1:bbasis$nbasis){
		sBias1dc[,(i-1)*bbasis$nbasis+j] = t(bvals)%*%(Phi$dPhi[,,1,2,(i-1)*bbasis$nbasis+j]*Sig)%*%rep(totm/nfine,nfine)
						- t(bvals)%*%(Phi$Phi[,,1,2]*(errs%*%t(bvals[,j]*XX)))%*%rep(totm/(n*nfine),nfine)
						- t(bvals)%*%(Phi$Phi[,,1,2]*(bvals[,j]*XX)%*%t(errs))%*%rep(totm/(n*nfine),nfine)
				
		sBias2dc[,(i-1)*bbasis$nbasis+j] = t(bvals)%*%(Phi$dPhi[,,2,2,(i-1)*bbasis$nbasis+j]*Sig)%*%rep(totm/nfine,nfine)
						- t(bvals)%*%(Phi$Phi[,,1,2]*(errs%*%t(bvals[,j]*XX)))%*%rep(totm/(n*nfine),nfine)
						- t(bvals)%*%(Phi$Phi[,,1,2]*(bvals[,j]*XX)%*%t(errs))%*%rep(totm/(n*nfine),nfine)
	 }
	}
	
	scoresdc = -regdc  - rbind(sBias1dc,sBias2dc)
	
	# Return both of these
	
	return(list(scores = scores,scoresdc=scoresdc))
}




# From Ed's code to generate the SMH


source(here::here("code/SHM-helper-functions.R"))

grid_range <- c(0, 2 * pi)
n_grid <- 500
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)

dat = generate_SHM_dataset(n = 100,
                        grid_points = grid_points,
                        sigma = 0.5,
                        mu_init = c(1, 0),
                        intensity = 1,
                        sigma_init = diag(rep(0.05, 2)))


# Translating it into Ramsay dialect. 
			
bbasis = create.bspline.basis(range(grid_points),nbasis = 51)

betaC = cbind(rep(-1,bbasis$nbasis),rep(0,bbasis$nbasis))
betafd = fd(betaC,bbasis)

Xfd = Data2fd(grid_points,dat$x,bbasis)




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


# Now let's just try this with a few less basis functions

betabasis = create.bspline.basis(range(tfine),nbasis=10)
#betabasis = create.constant.basis(range(tfine))

#betabasis = bbasis

betaC = c(rep(-1,betabasis$nbasis),rep(0,betabasis$nbasis))

betafd= fd( matrix(betaC,betabasis$nbasis,2),betabasis)

#hmm.small = nleqslv(betaC,ScoreFn,Xfd=Xfd,bbasis=betabasis,tfine=tfine,Sig=Sig)

#betaest.s = matrix(hmm.small$x,betabasis$nbasis,2)


betaPar = fdPar(betabasis,int2Lfd(2),0)

pdaest.small = fRegress(deriv.fd(Xfd,2),list(Xfd,deriv.fd(Xfd,1)),list(betaPar,betaPar))

betaest.pda = cbind(pdaest.small$betaestlist[[1]]$fd$coef,pdaest.small$betaestlist[[2]]$fd$coef)


X11()
plot(fd(betaest.pda,betabasis))
abline(h=-1)


# OK, let's get a gradient by finite differencing. This is pretty
# sticky so we may want to come up with something analytical. The
# gradient increases markedly as you decrease the size of the finite 
# difference. 
#xest = betaC

#Sest = ScoreFndc(xest,Xfd,betabasis,tfine,Sig)
#Gest = matrix(NA,2*betabasis$nbasis,2*betabasis$nbasis)

#for(i in 1:length(xest)){
#	txest = xest
#	txest[i] = txest[i] + 1e-2
	
#	Gest[,i] = (ScoreFndc(txest,Xfd,betabasis,tfine,Sig)$scores-Sest$scores)/1e-2
#}



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
abline(h=c(0,-1),col=c('red','blue'))
lines(fd(betaest.pda,betabasis),lty=4,lwd=3,col=c('blue','red'))
legend('bottomleft',c(expression(beta[0](t)),expression(beta[1](t)),'Score Est','PDA'),
	lty=c(1,1,1,4),lwd=3,col=c('darkblue','darkred','darkblue','red'),cex=1.5,ncol=2)
