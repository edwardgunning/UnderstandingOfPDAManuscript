require(fda)
require(deSolve)


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
