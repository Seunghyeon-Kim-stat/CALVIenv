
#' Generates the initial values
#' U: centered version of the response matrix
#' X: centered version of the predictor variable X
#' u: dimension of the envelope subspace. u can not be greater than total number of responses.
# @export
startingValue_genU<-function(U,X,u){
  fit<-lm(U~X-1)

  var_U=var(U)
  r=dim(U)[2]
  var_r=var(fit$residuals)
  gamma_start=FindGammaStart(var_U,var_r,u,t(fit$coef))
  beta_start=coef(fit)
  if(u==0){


    omega_start=NULL
    Q_gamma=diag(replicate(expr=1,n=r))
    gamma0_start=eigen(Q_gamma%*% var_U%*% Q_gamma)$vectors[,1:(r-u)]
    omega0_start=eigen(Q_gamma%*% var_U%*% Q_gamma)$values[1:(r-u)]
    O_start=gamma0_start

  }

  if((0<u)*(u<r)){
    omega_start=eigen(t(gamma_start)%*% var_r %*% gamma_start)$values
    Q_gamma=diag(replicate(expr=1,n=r))-gamma_start%*% t(gamma_start)
    gamma0_start=eigen(Q_gamma%*% var_U%*% Q_gamma)$vectors[,1:(r-u)]
    omega0_start=eigen(Q_gamma%*% var_U%*% Q_gamma)$values[1:(r-u)]
    O_start=cbind(gamma_start,gamma0_start)



  }

  if(u==r){
    omega_start=eigen(t(gamma_start)%*% var_r %*% gamma_start)$values

    gamma0_start=NULL
    omega0_start=NULL
    O_start=gamma_start



  }

  startpoint=list()
  startpoint[["omega"]]=omega_start
  startpoint[["omega0"]]=omega0_start
  startpoint$O=O_start
  if(   !is.null(gamma0_start)  ){
    startpoint$gamma0 <- gamma0_start;
  }
  if( !is.null(gamma_start)  ){
    startpoint$gamma=   as.matrix(gamma_start)

  }
  if(dim(X)[2]==1){startpoint$beta_start=t(as.matrix(beta_start))}
  if(dim(X)[2]>1){startpoint$beta_start=t(as.matrix(beta_start))}


  return(startpoint)
}

# See  cook s paper page 21
#Envelope Models for Parsimonious and E???cient Multivariate Linear Regression

fValue<-function(x,S1,S2){

  return(  log(det(t(x)%*%( S1) %*% x ))+log(det(t(x)%*%( S2) %*% x )))

}


fValue_alt<-function(x,S1,S2){

  return(  -log(.0001+det(t(x)%*%( S1) %*% x ))+log(.0001+det(t(x)%*%( S2) %*% x )))

}

##########################################################################
# @export
FindGammaStart<-function(var_U,var_r,u,beta_hat){
  if(u==0){return(NULL)}
  if(u>0){
    var_U_inv=try(solve(var_U),silent = TRUE)

    #print(var_U_inv)
    if(  ( is.numeric(var_U_inv))*(det(var_U)>0 )  ) {
      E1=eigen(var_U)$vectors
      E2=eigen(var_r)$vectors
      r=dim(E1)[1]
      allPair=combn((1:r),u)
      f1=f2=0
      for(i in 1:dim(allPair)[2]){
        f1[i]=fValue(E1[,allPair[,i]],var_U_inv,var_r)
        f2[i]=fValue(E2[,allPair[,i]],var_U_inv,var_r)

      }

    if(min(f1)<=min(f2)) {   Gamma=E1[,allPair[,which(f1==min(f1))]]     }
    if(min(f2)<min(f1)){   Gamma=E2[,allPair[,which(f2==min(f2))]]     }
    }


#     if(!is.numeric(var_U_inv)){
#
#      E1=eigen(var_U)$vectors
#     E2=eigen(var_r)$vectors
#     r=dim(E1)[1]
#     allPair=combn((1:r),u)
#     f1=f2=0
#      #print(f1,f2)
#     for(i in 1:dim(allPair)[2]){
#       f1[i]=fValue_alt(E1[,allPair[,i]],var_U,var_r)
#       f2[i]=fValue_alt(E2[,allPair[,i]],var_U,var_r)
#
#     }
#     #print(cbind(f1,f2))
#     if(min(f1)<=min(f2)) {   Gamma=E1[,allPair[,which(f1==min(f1))]]     }
#     if(min(f2)<min(f1)){   Gamma=E2[,allPair[,which(f2==min(f2))]]     }
#      #Gamma=E1[,1:u]
#   }


    if(! ( is.numeric(var_U_inv))* (det(var_U)>0)){

      E1=eigen(var_U)$vectors
      E2=eigen(var_r)$vectors
      r=dim(E1)[1]
      allPair=combn((1:r),u)
      f1=f2=0
      #print(f1,f2)
      for(i in 1:dim(allPair)[2]){
        f1[i]=norm(t(E1[,allPair[,i]])%*%beta_hat)
        f2[i]= norm(t(E2[,allPair[,i]])%*%beta_hat)

      }
      #print(cbind(f1,f2))
      if(max(f1) >= max(f2)) {   Gamma=E1[,allPair[,which(f1==max(f1))]]     }
      if(max(f2)> max(f1)){   Gamma=E2[,allPair[,which(f2==max(f2))]]     }
      #Gamma=E1[,1:u]
    }

    return(Gamma)
  }

}




