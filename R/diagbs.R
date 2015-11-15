#' Diagnostic Analysis - Local Influnce
#' @description Diagnostics for the RBS model
#'
#'@usage envelope(model,k=19,alpha=0.05,res="deviance", precision = "fixed")
#'
#' @param model object of class \code{gamlss} holding the fitted model.
#' @param k number of replications for envelope construction. Default is 19.
#' @param alpha value giving the size of the envelope. Default is 0.05 which is equivalent to a 95\% band.
#' @param res type of residuals to be extracted. Default is deviance.
#' @param precision If \code{precision = "fixed"} a model with precision fixed is used;
#' else a model with precision variable is used.
#'
#'@return A simulated envelope of the class RBS.
#'
#' @author
#'Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}, F.J.A. Cysneiros \email{cysneiros@de.ufpe.br}, Victor Leiva \email{victorleivasanchez@gmail.com} and Michelli Barros \email{michelli.karinne@gmail.com}
#'
#'@references
#'Atkinson, A. C. (1985) Plots, transformations and regression : an introduction to graphical methods of diagnostic regression analysis. Oxford Science Publications, Oxford.
#'@export

diag.bs=function(model,npoints=NULL,plot.diag = "FALSE",perturbation="PC")
{

  X= model.matrix(model)# model matrix
  p=ncol(X) #number parameters
  n=nrow(X) #number observations
  vp=c(model$mu.coefficients,model$sigma.coefficients)   #estimate parameters
  Y=model$y  # response variable
  H = -solve(vcov(model)) # Hessian Matrix !!!!!!Estou Aqui!!!!!
  L1 = vcov(model) # Hessian Matrix inverse
  link = model$mu.link
  linkstr = link
  linkobj = make.link(linkstr)
  linkfun = linkobj$linkfun
  linkinv = linkobj$linkinv
  mu.eta = linkobj$mu.eta


  ################################################## Curvatures #############################################################

  DeltaPC=function(vp)
  {
    betab = vp[1:p]
    delta = vp[-(1:p)]
    eta = as.vector(X%*%betab)
    mu = linkinv(eta)
    vt = Y
    vu=mu
    vd=delta
    Da = diag(mu.eta(eta))
    ve=(-1/(2*vu)) + vd /((vt*vd) + vt + (vd*vu)) +  ((vd+1)*vt)/(4*(vu^2)) - (vd^2)/(4*vt*(vd+1)) #ok!
    vb= 1/2 - (1/2)*(vd+1)^(-1) + (vt+vu)/((vt*vd) + vt + (vd*vu)) - (1/4)*(vt/vu) - (vd*(vd+2)*vu)/(4*vt*((vd+1)^2)) #ok!
    De=diag(ve)
    Deltab=t(X)%*%Da%*%De
    Deltad=t(as.matrix(vb))
    Cpd=rbind(Deltab,Deltad)
    return(Cpd)
  }


  DeltaPR=function(vp)
  {
    betab = vp[1:p]
    delta = vp[-(1:p)]
    eta = as.vector(X%*%betab)
    mu = linkinv(eta)
    vt = Y
    vu=mu
    vd=delta
    Da = diag(mu.eta(eta))
    phi=((2*vd)+5)/((vd+1)^2)
    vk=sqrt((vu^2)*phi)
    Dk=diag(vk)
    vpsi=-(vd*(vd+1))/( ((vd*vt)  + vt  + (vd*vu) )^2) + (vd+1)/(4*(vu^2)) +  ((vd^2)/(4*(vd+1)*(vt^2))) #ok!!
    Dpsi=diag(vpsi)
    vro=(-vu/( ((vd*vt) + vt + (vd*vu))^2))  -  1/(4*vu) + (vd*(vd+2)*vu)/(4*(vt^2)*((vd+1)^2))  #ok!!
    Deltab=t(X)%*%Da%*%Dpsi%*%Dk
    Deltad=t(as.matrix(vro))%*%Dk
    Delta= rbind(Deltab,Deltad)
    return(Delta)
  }


  DeltaPP=function(vp)
  {
    betab = vp[1:p]
    delta = vp[-(1:p)]
    eta = as.vector(X%*%betab)
    mu = linkinv(eta)
    vt = Y
    vu=mu
    vd=delta
    Da = diag(mu.eta(eta))
    qsi=-(vd*vt)/(((vd*vt) + vt +(vd*vu))^2) - (vd*vt)/(4*(vu^2)) + ((vd^2)*(vd+2))/(4*vt*((vd+1)^2))
    vpsi=-(1/2) + 1/(2*((vd+1)^2)) - (vt*(vt+vu))/(((vd*vt)+vt+(vd*vu))^2) + vt/(4*vu)+ ((vd^2)*vu*(vd+3))/(4*vt*((vd+1)^3)) + (vd*vu)/(vt*((vd+1)^3))
    Dqsi=diag(qsi)
    Deltab=t(X)%*%Da%*%Dqsi
    Deltad=t(as.matrix(vpsi))
    Delta= rbind(Deltab,Deltad)
    return(Delta)
  }

  ################################################ Alavancagem generalizadas  ############################################

  Lbetat=function(vp)
  {
    betab = vp[1:p]
    delta = vp[-(1:p)]
    eta = as.vector(X%*%betab)
    mu = linkinv(eta)
    vt = Y
    vu=mu
    vd=delta
    Da = diag(mu.eta(eta))
    phi=(2*((vd+1)^2))/((2*vd) +5)
    vpsi= -(vd*(vd+1))/( ((vd*vt)  + vt  + (vd*vu) )^2) + (vd+1)/(4*(vu^2)) +  ((vd^2)/(4*(vd+1)))*(1/(vt^2))
    Dpsi=diag(vpsi)
    vro= (-vu/( ((vd*vt) + vt + (vd*vu))^2))  +  1/(4*vu^2) + (vd*(vd+2)*vu)/(4*(vt^2)*((vd+1)^2))
    Deltab= t(X)%*%Da%*%Dpsi
    Deltad= t(as.matrix(vro))
    Delta= rbind(Deltab,Deltad)
    return(Deltab)
  }

  B=function(Delta,I,M)
  {
    #I ? o inverso da matrix hessiana
    #M determina sobre qual sobre quais cojuntos de par?metros deseja-se estudar o influ?ncia
    #Delta ? a curvatura que determina qual o tipo de perturba??o deseja-se estudar
    B=(t(Delta)%*%(I-M)%*%Delta)
    return(B)
  }

  caseweights=DeltaPC(vp)  #Case-weights perturbation
  responsepert=DeltaPR(vp) #Response perturbation
  precisionpert=DeltaPP(vp) #Pertubation on precision parameter

  #################################################Matrizes Auxiliares#############################
  Ldelta= H[p+1,p+1]
  Lbeta=H[1:p,1:p]
  b11=cbind(matrix(0, p, p), matrix(0, p, 1))
  b12=cbind(matrix(0, 1, p), -Ldelta^(-1))
  B1= rbind(b11, b12)  #parameter beta
  b211 =cbind(-solve(Lbeta), matrix(0, p, 1))
  b212= cbind(matrix(0, 1, p), matrix(0, 1, 1))
  B2=rbind(b211,b212)  # parameter delta

  b311 =cbind(matrix(0, p, p), matrix(0, p, 1))
  b312= cbind(matrix(0, 1, p), matrix(0, 1, 1))
  B3=rbind(b311,b312)  # parameter theta


  ##################################################Calculation of dmax,Ci,hii#########################
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Pertubation on precision parameter$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  FPP=B(precisionpert,L1,B3)           #theta
  autovmaxthetaPP=eigen(FPP,symmetric=TRUE)$val[1]
  vetorpcthetaPP=eigen(FPP,symmetric=TRUE)$vec[,1]

  FPP1=B(precisionpert,L1,B1)       #Beta
  autovmaxbPP=eigen(FPP1,symmetric=TRUE)$val[1]
  vetorpcbPP=eigen(FPP1,symmetric=TRUE)$vec[,1]

  FPP2=B(precisionpert,L1,B2)      #delta
  autovmaxdPP=eigen(FPP2,symmetric=TRUE)$val[1]
  vetorpcdPP=eigen(FPP2,symmetric=TRUE)$vec[,1]


  vCithetaPP=2*abs(diag(FPP))
  vCiPP=2*abs(diag(FPP1))
  vCidPP=2*abs(diag(FPP2))
  ###################################################################################################################
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Response perturbation$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  FPR=B(responsepert,L1,B3)           #Theta
  autovmaxthetaPR=eigen(FPR,symmetric=TRUE)$val[1]
  vetorpcthetaPR=eigen(FPR,symmetric=TRUE)$vec[,1]

  FPR1=B(responsepert,L1,B1)       #beta
  autovmaxbPR=eigen(FPR1,symmetric=TRUE)$val[1]
  vetorpcbPR=eigen(FPR1,symmetric=TRUE)$vec[,1]

  FPR2=B(responsepert,L1,B2)     #delta
  autovmaxdPR=eigen(FPR2,symmetric=TRUE)$val[1]
  vetorpcdPR=eigen(FPR2,symmetric=TRUE)$vec[,1]

  vCithetaPR=2*abs(diag(FPR))
  vCiPR=2*abs(diag(FPR1))
  vCidPR=2*abs(diag(FPR2))
  ############################################################################################################################
  ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Pondera??o de casos$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  FPC=B(caseweights,L1,B3)            #Theta
  autovmaxthetaPC=eigen(FPC)$val[1]
  vetorpcthetaPC=eigen(FPC)$vec[,1]

  FPC1=B(caseweights,L1,B1)      #beta
  autovmaxbPC=eigen(FPC1)$val[1]
  vetorpcbPC=eigen(FPC1)$vec[,1]

  FPC2=B(caseweights,L1,B2)     #delta
  autovmaxdPC=eigen(FPC2)$val[1]
  vetorpcdPC=eigen(FPC2)$vec[,1]

  vCithetaPC=2*abs(diag(FPC))
  vCiPC=2*abs(diag(FPC1))
  vCidPC=2*abs(diag(FPC2))

  ################################################### #############################################################
  yest=muest=model$mu.fv #fitted.values
  Lby=Lbetat(vp)
  betas=model$mu.coefficients
  vO=rep(0,n)
  eta = as.vector(X%*%betas)
  Da = diag(mu.eta(eta))
  Do=(Da%*%X)
  GL=diag(Do%*%(L1[1:p,1:p])%*%Lby )    #Alavancagem Generalizada

  ####################################################Graphics#####################################################

  res1=residuals(model, res="pearson")  #residuals 1
  res2=residuals(model, res="score")  #residuals 2
  res3=residuals(model, res="deviance")  #residuals 3

  ############################################### number points ############################
  npointsr1=npoints[1]
  npointsr2=npoints[2]
  npointsr3=npoints[3]
  npointsdt=npoints[4]
  npointsdb=npoints[5]
  npointsdd=npoints[6]
  npointsCbt=npoints[7]
  npointsCbb=npoints[8]
  npointsCbd=npoints[9]
  npointsGL=npoints[10]

  if(perturbation == "PC")
  {

    #############################################dmax e Ci###########################################################
    #dmax
    infl=vetorpcthetaPC/sqrt(autovmaxthetaPC)
    dmaxG=abs(infl)

    infl1=vetorpcbPC/sqrt(autovmaxbPC)
    dmaxG1=abs(infl1)

    infl2=vetorpcdPC/sqrt(autovmaxdPC)
    dmaxG2=abs(infl2)

    #Ci
    Cb=vCithetaPC
    Cb=Cb/sum(Cb)

    Cb1=vCiPC
    Cb1=Cb1/sum(Cb1)

    Cb2=vCidPC
    Cb2=Cb2/sum(Cb2)
    ################################################################################################################

    if(plot.diag == "TRUE")
    {
      ########################################################Residuals ###############################################
      x11()
      plot(yest,res1,pch=8,lwd=2,xlab="fv",ylab="r1",sub="",cex=0.7)
      if(!is.null(npoints) && npointsr1 != 0) identify(yest,res1,n=npointsr1)
      x11()
      plot(yest,res2,pch=8,lwd=2,xlab="fv",ylab="r2",sub="",cex=0.7)
      if(!is.null(npoints) && npointsr2 !=0) identify(yest,res2,n=npointsr2)
      x11()
      plot(yest,res3,pch=8,lwd=2,xlab="fv",ylab="r3",sub="",cex=0.7)
      if(!is.null(npoints) && npointsr3 !=0 )  identify(yest,res3,n=npointsr3)

      ############################################## Alavancagem Generalizada ####################################################3
      x11()
      plot(yest,GL,pch=8,xlab="fv",ylab="al",cex=0.7)
      if(!is.null(npoints) && npointsGL !=0) identify(yest,GL,n=npointsGL)
      ####################################################dmax#######################################################
      x11()
      plot(1:length(dmaxG),dmaxG,pch=8,xlab="id",ylim=c(0,1),ylab="dmaxpctheta",sub="",cex=0.7)
      if(!is.null(npoints)&& npointsdt !=0) identify(1:length(dmaxG),dmaxG,n=npointsdt)

      x11()
      plot(1:length(dmaxG1),dmaxG1,pch=8,xlab="id",ylim=c(0,1),ylab="dmaxpcbeta",sub="",cex=0.7)
      if(!is.null(npoints) &&  npointsdb !=0) identify(1:length(dmaxG1),dmaxG1,n=npointsdb)

      x11()
      plot(1:length(dmaxG2),dmaxG2,pch=8,xlab="id",ylim=c(0,1),ylab="dmaxpcdelta",sub="",cex=0.7)
      if(!is.null(npoints) && npointsdd !=0) identify(1:length(dmaxG2),dmaxG2,n=npointsdd)

      ########################################################Ci#########################################################
      x11()
      limi=2*mean(Cb)
      plot(1:length(Cb),Cb,pch=8,xlab="id",ylab="Ctheta",sub="",ylim=c(0,1),cex=0.7)
      abline(h=limi,lty=2,lwd=2)
      if(!is.null(npoints) && npointsCbt !=0) identify(1:length(Cb),Cb,n=npointsCbt)

      x11()
      limi=2*mean(Cb1)
      plot(1:length(Cb1),Cb1,pch=8,xlab="id",ylab="Cbeta",sub="",ylim=c(0,1),cex=0.7)
      abline(h=limi,lty=2,lwd=2)
      if(!is.null(npoints) && npointsCbb !=0) identify(1:length(Cb1),Cb1,n=npointsCbb)

      x11()
      limi=2*mean(Cb2)
      plot(1:length(Cb2),Cb2,pch=8,xlab="id",ylab="Cdelta",sub="",ylim=c(0,1),cex=0.7)
      abline(h=limi,lty=2,lwd=2)
      if(!is.null(npoints) && npointsCbd !=0) identify(1:length(Cb2),Cb2,n=npointsCbd)

    }

    else{
      diagbs=list(dmax = list(theta = vetorpcthetaPC , beta = vetorpcbPC , delta = vetorpcdPC), Ci = list(theta = Cb,beta=Cb1 ,delta=Cb2),GL=GL)
      return(diagbs)
    }

  }

  else{

    #dmax
    infl=vetorpcthetaPR/sqrt(autovmaxthetaPR)
    dmaxG=abs(infl)

    infl1=vetorpcbPR/sqrt(autovmaxthetaPR)
    dmaxG1=abs(infl1)

    infl2=vetorpcdPR/sqrt(autovmaxthetaPR)
    dmaxG2=abs(infl2)

    #Ci
    Cb=vCithetaPR
    Cb=Cb/sum(Cb)

    Cb1=vCiPR
    Cb1=Cb1/sum(Cb1)

    Cb2=vCidPR
    Cb2=Cb2/sum(Cb2)

    if(plot.diag == "TRUE")
    {
      ########################################################Residuals ###############################################
      x11()
      plot(yest,res1,pch=8,lwd=2,xlab="fv",ylab="r1",sub="",cex=0.7)
      if(!is.null(npoints) && npointsr1 != 0) identify(yest,res1,n=npointsr1)
      x11()
      plot(yest,res2,pch=8,lwd=2,xlab="fv",ylab="r2",sub="",cex=0.7)
      if(!is.null(npoints) && npointsr2 !=0) identify(yest,res2,n=npointsr2)
      x11()
      plot(yest,res3,pch=8,lwd=2,xlab="fv",ylab="r3",sub="",cex=0.7)
      if(!is.null(npoints) && npointsr3 !=0 )  identify(yest,res3,n=npointsr3)

      ############################################## Alavancagem Generalizada ####################################################3
      x11()
      plot(yest,GL,pch=8,xlab="fv",ylab="al",cex=0.7)
      if(!is.null(npoints) && npointsGL !=0) identify(yest,GL,n=npointsGL)


      ####################################################dmax#######################################################
      x11()
      plot(1:length(dmaxG),dmaxG,pch=8,xlab="id",ylim=c(0,1),ylab="dmaxpctheta",sub="",cex=0.7)
      if(!is.null(npoints)&& npointsdt !=0) identify(1:length(dmaxG),dmaxG,n=npointsdt)

      x11()
      plot(1:length(dmaxG1),dmaxG1,pch=8,xlab="id",ylim=c(0,1),ylab="dmaxpcbeta",sub="",cex=0.7)
      if(!is.null(npoints) &&  npointsdb !=0) identify(1:length(dmaxG1),dmaxG1,n=npointsdb)

      x11()
      plot(1:length(dmaxG2),dmaxG2,pch=8,xlab="id",ylim=c(0,1),ylab="dmaxpcdelta",sub="",cex=0.7)
      if(!is.null(npoints) && npointsdd !=0) identify(1:length(dmaxG2),dmaxG2,n=npointsdd)

      ########################################################Ci#########################################################
      x11()
      limi=2*mean(Cb)
      plot(1:length(Cb),Cb,pch=8,xlab="id",ylab="Ctheta",sub="",ylim=c(0,1),cex=0.7)
      abline(h=limi,lty=2,lwd=2)
      if(!is.null(npoints) && npointsCbt !=0) identify(1:length(Cb),Cb,n=npointsCbt)

      x11()
      limi=2*mean(Cb1)
      plot(1:length(Cb1),Cb1,pch=8,xlab="id",ylab="Cbeta",sub="",ylim=c(0,1),cex=0.7)
      abline(h=limi,lty=2,lwd=2)
      if(!is.null(npoints) && npointsCbb !=0) identify(1:length(Cb1),Cb1,n=npointsCbb)

      x11()
      limi=2*mean(Cb2)
      plot(1:length(Cb2),Cb2,pch=8,xlab="id",ylab="Cdelta",sub="",ylim=c(0,1),cex=0.7)
      abline(h=limi,lty=2,lwd=2)
      if(!is.null(npoints) && npointsCbd !=0) identify(1:length(Cb2),Cb2,n=npointsCbd)
    }

    else{
      diagbs=list(dmax = list(theta = dmaxG , beta = dmaxG1 , delta = dmaxG2), Ci = list(theta = Cb, beta= Cb1 ,delta=Cb2),GL=GL)
      return(diagbs)
    }

  }

}
