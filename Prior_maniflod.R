#' Prior constructor
#' @export
Benvlp_prior<-function(number_of_responses=NULL, number_of_covariates=NULL, response_envelope_dim=NULL,
                     eta_hyper_e= NULL,
                         eta_hyper_C=NULL,
                         O_hyper_D=NULL,
                         O_hyper_G=NULL,
                         omega_hyper_shape=NULL,
                         omega_hyper_rate=NULL,
                         omega0_hyper_shape=NULL,
                         omega0_hyper_rate=NULL, model="Manifold"){

# Checking for  number_of_responses and response_envelope_dim conditions
  error_flag=FALSE
  if(is.null(number_of_responses)){ stop("number_of_responses: Number of the response variables must be specified to construct the prior.") }
  if(is.null(response_envelope_dim)){ stop(" envelope_dim: Dimension of the response envelope space must be specified to construct the prior.") }
  if(is.null(number_of_covariates)){stop( "number_of_covariates: Number of covariates must specified to construct the prior." )}


  if(!is.numeric(number_of_responses)){stop("number_of_responses: must be an integer.")}
  if(!is.numeric(response_envelope_dim)){stop("response_envelope_dim: must be an integer between 0 and 'number_of_responses' ")}
  if(!is.numeric(number_of_covariates)){stop( "number_of_covariates: Number of covariates must Integer." )}

  number_of_responses=as.integer(number_of_responses); number_of_covariates=as.integer(number_of_covariates); response_envelope_dim=as.integer(response_envelope_dim)
  if( number_of_responses<=0 ){stop("number_of_responses: NUmber of responses must be a nonzero integer.")}
  if( response_envelope_dim<0 ){stop("response_envelope_dim: dimension of the envelope space must be a nonnegative integer.")}
  if( response_envelope_dim>number_of_responses ){stop(" response_envelope_dim<=number_of_responses: dimension of the envelope space must be a nonnegative integer less than equal to the number of responses.")}

  if( number_of_covariates <=0 ){ stop(" number_of_covariates: must be a positive integer. ")}


  ######### ######### ######### #########
  ######### ######### ######### #########
# Checking conditions on the hyper priors of the Omega and Omega0:
  ######### ######### ######### #########
  ######### ######### ######### #########
  if(!is.numeric(omega_hyper_shape)){stop('The variable omega_hyper_shape: must be numeric.')}
  if(!is.numeric(omega_hyper_rate)){stop('The variable omega_hyper_rate: must be numeric.')}
  if(!is.numeric(    omega0_hyper_shape)){stop('The variable     omega0_hyper_shape: must be numeric.')}
  if(!is.numeric(   omega0_hyper_rate)){stop('The variable    omega0_hyper_rate: must be numeric.')}

#
  omega_hyper_shape=as.vector(omega_hyper_shape);omega_hyper_rate=as.vector(omega_hyper_rate)
  omega0_hyper_shape=as.vector(omega0_hyper_shape); omega0_hyper_rate=as.vector(omega0_hyper_rate)
#browser()
  if(response_envelope_dim>=0){
    if(!(length(omega_hyper_shape)==response_envelope_dim)){  stop("length of the numeric vector omega_hyper_shape must be equal to the  response_envelope_dim. omega_hyper_shape must be set to NULL in case response_envelope_dim=0") }
    if(!(length(omega_hyper_rate)==response_envelope_dim)){  stop("length of the numeric vector omega_hyper_rate must be equal to the  response_envelope_dim. omega_hyper_rate must be set to NULL in case response_envelope_dim.=0 ") }

    if( any(omega_hyper_shape < (-1-number_of_covariates/2 )) ){ stop("Each component of the omega_hyper_shape parameters are required to be greater than equl to -1-p/2 where p=number_of_covariates."); error_flag=TRUE  }
    if( any(omega_hyper_rate <0)  ){ stop("Each component of the omega_hyper_rate parameters are required to be non-negative.")  }

    }

  if(response_envelope_dim<=number_of_responses){
    dim_orthogonal_envelope_space=number_of_responses-response_envelope_dim;
    if(!(length(omega0_hyper_shape)==dim_orthogonal_envelope_space)){  stop("length of the numeric vector omega0_hyper_shape must be equal to the  number_of_responses-response_envelope_dim. omega0_hyper_shape must be set to NULL in case number_of_responses=response_envelope_dim ") }
    if(!(length(omega0_hyper_rate)==dim_orthogonal_envelope_space)){  stop("length of the numeric vector omega0_hyper_rate must be equal to the  number_of_responses-response_envelope_dim. omega0_hyper_rate must be set to NULL in case number_of_responses=response_envelope_dim ") }

    if( any(omega0_hyper_shape <  (-1))  ){ stop("Each component of the omega0_hyper_shape parameters are required to be greater than equals to -1."); error_flag=TRUE }
    if( prod(omega0_hyper_rate >=0)==0  ){ stop("Each component of the omega0_hyper_rate parameters are required to be non-negative.") }

    }
  ######### ######### ######### #########
  ######### ######### ######### #########
  ######### Hyper parameter on eta
  ######### ######### ######### #########
  ######### ######### ######### #########
   #browser()
 # if(response_envelope_dim==0){
  #  eta_hyper_C=(diag(0));
  #}
  if(response_envelope_dim>0){
  if(!is.null(eta_hyper_C)){
    eta_hyper_C=as.matrix(eta_hyper_C)
    if(!is.numeric(eta_hyper_C)){stop(" eta_hyper_C: must be numeric u times u matrix where u is equal to the dimension of the response envelope space. ")}
    if( prod(dim(eta_hyper_C)==c(number_of_covariates,number_of_covariates))!=1 ){stop("Dimension of the matrix eta_hyper_C has to be u times u where u=response_envelope_dim ")}
    if(any((eigen(eta_hyper_C)$value)<0)){ stop("eta_hyper_C must be a nonnegative definite matrix.")}
  }
  }

  if(is.null(eta_hyper_C)){
    print("Unspecified parameter eta_hyper_C, inverse of a variance covariance matrix for the prior distribution of the parameter eta. Uniform improper prior specification will be considered as a default stratigy.")
    eta_hyper_C=0*diag(response_envelope_dim)
    #eta_hyper_e= matrix(0, nrow=number_of_responses, ncol=number_of_covariates)
  }

  ######### ######### ######### ######### ######### ######### ######### #########
  ## eta_hyper_e
  ######### ######### ######### ######### ######### ######### ######### #########
  if(all(eta_hyper_C==0)){
    eta_hyper_e= matrix(0, nrow=number_of_responses, ncol=number_of_covariates)
  }
  if(!all(eta_hyper_C==0)){
    if(!is.null(eta_hyper_e)){
             if(!is.numeric(eta_hyper_e)){stop(" eta_hyper_e must be numeric matrix. non numeric values are allowed.")}
              eta_hyper_e=as.matrix(eta_hyper_e)
              if(any(dim(eta_hyper_e)!=c(number_of_responses,number_of_covariates))){ stop("dimension of the numeric matrix eta_hyper_e mtst be number_of_responses \times number_of_covariates ")}

    }
    if(is.null(eta_hyper_e)){    eta_hyper_e= matrix(0, nrow=number_of_responses, ncol=number_of_covariates) }
  }


  ######## ######### ######### ######### ##########
  ######### ######### ######### ######### #########
  ###Hyper parameter on O ##
  ######### ######### ######### ######### #########
  ######### ######### ######### ######### #########
  ######## ######### ######### ######### ##########


  if(!is.null(O_hyper_G)){
    if(!is.numeric(O_hyper_G)){ stop(" O_hyper_G must be numeric. ")}
    O_hyper_G=as.matrix(O_hyper_G)
    if( any(dim(O_hyper_G) != c(number_of_responses, number_of_responses)) ) {stop("dimension of the input O_hyper_G must be r\times r where r=number_of_responses. ") }
    if(is.complex(eigen(O_hyper_G)$values)){stop(" O_hyper_G must be a nonnegetive definite matrix. The current input of O_hyper_Gleads to   complex eigen values.i.e. eigen(O_hyper_G) returns complex eigrnvalues, that is not appropriate. ")}
    if(any(  (eigen(O_hyper_G)$values)<0) ){ stop("At least one of the eigenvalues of the input argument O_hyper_G apears to be negative. The intut values for  O_hyper_G must be a nonnegative definite matrix.")}

  }


  if(is.null(O_hyper_G)){
    O_hyper_G=matrix(0, ncol=number_of_responses, nrow=number_of_responses)
  }
  ######### ######### ######### ######### #########
  ######## ######### ######### ######### ##########


  if(all(O_hyper_G==0)){
    O_hyper_D=diag(number_of_responses);
  }


  if(!is.null(O_hyper_D)){
      if(!isDiagonal(O_hyper_D)){stop("O_hyper_D must be a diagonal matrix")}
      if(any(dim(O_hyper_D)!=c(number_of_responses,number_of_responses))){stop(" O_hyper_D must be a diagonal matrix of dimension r times r where r=number_of_responses.")}
      if(any(diag(O_hyper_D<=0))){stop("O_hyper_D must be a diagonal matrix with positive diagonal  elements.")}
  }

  if(is.null(O_hyper_D)){
    O_hyper_D=diag(number_of_responses);
  }


  ######### ######### ######### ######### #########
  ######### ######### ######### ######### #########
  ######## ######### ######### ######### ##########
  ######### ######### ######### ######### #########
  ######### ######### ######### ######### #########
  ######## ######### ######### ######### ##########

  #if( length(omega_hyper_shape)!=omega_hyper_rate){print(" Number of shape parameter and the rate parameter has to be same for the parameter omega. ");error_flag=TRUE}
  #if( length(omega0_hyper_shape)!=omega_hyper_rate){print(" Number of shape parameter and the rate parameter has to be same for omga0. "); error_flag=TRUE}

  #if( (length(omega_hyper_shape)  + length(omega0_hyper_shape) )!= dim(O_hyper_D) ){ print(" Dimension of the O_hyper_D and total number of the omega and omega0 has to match");error_flag=TRUE}

  if(error_flag){ obj=NULL   }
  if(!error_flag){
   obj=list(  eta_hyper_e= eta_hyper_e,
              eta_hyper_C=eta_hyper_C,
              O_hyper_D=O_hyper_D,
              O_hyper_G=O_hyper_G,
              omega_hyper_shape=omega_hyper_shape,
              omega_hyper_rate=omega_hyper_rate,
              omega0_hyper_shape=omega0_hyper_shape,
              omega0_hyper_rate=omega0_hyper_rate , model=model)

     class(obj)="Benvlp_prior"
    }
    return(obj)
}


##################################################################
####################################################################
####################################################################
#' @export
Benvlp_prior_uniform<-function( number_of_responses=NULL, number_of_covariates=NULL, response_envelope_dim=NULL ){

  if(is.null(number_of_responses)){ stop("number_of_responses: Number of the response variables must be specified to construct the prior.") }
  if(is.null(response_envelope_dim)){ stop(" envelope_dim: Dimension of the response envelope space must be specified to construct the prior.") }
  if(is.null(number_of_covariates)){stop( "number_of_covariates: Number of covariates must specified to construct the prior." )}


  if(!is.numeric(number_of_responses)){stop("number_of_responses: must be an integer.")}
  if(!is.numeric(response_envelope_dim)){stop("response_envelope_dim: must be an integer between 0 and 'number_of_responses' ")}
  if(!is.numeric(number_of_covariates)){stop( "number_of_covariates: Number of covariates must Integer." )}

  number_of_responses=as.integer(number_of_responses); number_of_covariates=as.integer(number_of_covariates); response_envelope_dim=as.integer(response_envelope_dim)
  if( number_of_responses<=0 ){stop("number_of_responses: NUmber of responses must be a nonzero integer.")}
  if( response_envelope_dim<0 ){stop("response_envelope_dim: dimension of the envelope space must be a nonnegative integer.")}
  if( response_envelope_dim>number_of_responses ){stop(" response_envelope_dim<=number_of_responses: dimension of the envelope space must be a nonnegative integer less than equal to the number of responses.")}

  if( number_of_covariates <=0 ){ stop(" number_of_covariates: must be a positive integer. ")}

  ########################################################
  ########  Hyper parameter on Omega and omega0 ##########
  ########################################################
  if(response_envelope_dim==0){omega_hyper_shape=NULL;omega_hyper_rate=NULL}
  if(response_envelope_dim>=0){ omega_hyper_shape=rep(x =  (-1-number_of_covariates/2), length.out=response_envelope_dim); omega_hyper_rate=rep(x = 0, length.out=response_envelope_dim) }

  #browser()
  if(response_envelope_dim==number_of_responses){omega0_hyper_shape=NULL;omega0_hyper_rate=NULL}
  if(response_envelope_dim<=number_of_responses){ omega0_hyper_shape=rep(x = -1, length.out=(number_of_responses-response_envelope_dim)); omega0_hyper_rate=rep(x = 0, length.out=(number_of_responses-response_envelope_dim)) }


  ########################################################
  #########   Hyper parameter on Eta ####################
  #######################################################
  eta_hyper_C=0*diag(number_of_covariates)
  eta_hyper_e= matrix(0, nrow=number_of_responses, ncol=number_of_covariates)

  ######## ######### ######### ######### ##########
  ######### ######### ######### ######### #########
  ###Hyper parameter on O ##
  ######### ######### ######### ######### #########
  ######### ######### ######### ######### #########
  ######## ######### ######### ######### ##########
  O_hyper_G=matrix(0, ncol=number_of_responses, nrow=number_of_responses)
  O_hyper_D=diag(number_of_responses)


  prior_object<- Benvlp_prior(number_of_responses=number_of_responses, number_of_covariates=number_of_covariates, response_envelope_dim=response_envelope_dim,    eta_hyper_e= eta_hyper_e,
                            eta_hyper_C=eta_hyper_C,
                            O_hyper_D=O_hyper_D,
                            O_hyper_G=O_hyper_G,
                            omega_hyper_shape=omega_hyper_shape,
                            omega_hyper_rate=omega_hyper_rate,
                            omega0_hyper_shape=omega0_hyper_shape,
                            omega0_hyper_rate=omega0_hyper_rate, model="Manifold")

  return(prior_object)
}

##################################################################################################################
##################################################################################################################
##################################################################################################################
# print function for the prior
#' @export
Benvlp_prior_empirical<-function(Y, X, response_envelope_dim,comfidence_level="LOW"){

  ### Need to add check on Y and X

  ###############################
  Y=as.matrix(Y); X=as.matrix(X)
  n=dim(Y)[1]; number_of_responses=dim(Y)[2]; number_of_covariates=dim(X)[2]

  if(dim(Y)[1]!=dim(X)[1]){stop("Y and X are incompatible for Envelope Regression. number of rows for Y and number of rows of X does not match.")}
  ##### Che k on response_envelope_dim ########
  if(is.null(response_envelope_dim)){ stop(" envelope_dim: Dimension of the response envelope space must be specified to construct the prior.") }
  if(!is.numeric(response_envelope_dim)){stop("response_envelope_dim: must be an integer between 0 and 'number_of_responses' ")}
  response_envelope_dim=as.integer(response_envelope_dim)
  if( response_envelope_dim<0 ){stop("response_envelope_dim: dimension of the envelope space must be a nonnegative integer.")}
  if( response_envelope_dim>number_of_responses ){stop(" response_envelope_dim<=number_of_responses: dimension of the envelope space must be a nonnegative integer less than equal to the number of responses.")}
  #########################################


  ############################################
  ############################################
  ## Extracton of the information from data ##
  ###############***************##############
  lst=startingValue_genU(Y,X,response_envelope_dim)
  ###############***************#############
  gamma=lst[["gamma"]]
  gamma0=lst[["gamma0"]]
  omega=lst[["omega"]]
  omega0=lst[["omega0"]]
  O=cbind(gamma,gamma0)
  ###########################################
  ###########################################

  r=number_of_responses
  u=response_envelope_dim
  p=number_of_covariates

  ########################################################
  #########   Hyper parameter on Eta #####################
  ########################################################
  eta_hyper_C=0*diag(number_of_covariates)
  if(r>n){C_prior=diag(number_of_covariates)}
  eta_hyper_e= lst$beta_start

  ######## ######### ######### ######### ################
  ######### ######### ######### ######### ###############
  ### Hyper parameter on O ##
  ######### ######### ######### ######### ###############
  ######### ######### ######### ######### ###############
 # O_hyper_G=G_prior_EB_selection(O,omega,omega0)
  O_hyper_G=G_prior_EB_selection_new(O = O,omega = omega,omega0 = omega0,sample_size = n,comfidence_level=comfidence_level)

  if(length(c(omega,omega0))>1){O_hyper_D=diag(c(omega,omega0))}
  if(length(c(omega,omega0))==1){O_hyper_D=(c(omega,omega0))}
  ########################################################
  ########  Hyper parameter on Omega and omega0 ##########
  ########################################################



  if(u==0){
    omega_hyper_rate=numeric(0);omega_hyper_shape=numeric(0)
    e_prior=matrix(0,nrow=r,ncol=p) ######
    omega0_hyper_rate_common=mean(1/omega0)/(  mean(1/omega0^2)-(mean(1/omega0))^2 )
    omega0_hyper_shape_common=mean(1/omega0)*omega0_hyper_rate_common

    #####
    omega0_hyper_rate=rep(x =omega0_hyper_rate_common , length.out=(r-u))
    omega0_hyper_shape=rep(x =omega0_hyper_shape_common , length.out=(r-u))
  }

  if(u==1){
    omega0_hyper_rate_common=mean(1/omega0)/(  mean(1/omega0^2)-(mean(1/omega0))^2 )
    omega0_hyper_shape_common=mean(1/omega0)*omega0_hyper_rate_common
    #####
    omega_hyper_shape=1/omega
    omega_hyper_rate= 1
    omega0_hyper_rate=rep(x =omega0_hyper_rate_common , length.out=(r-u))
    omega0_hyper_shape=rep(x =omega0_hyper_shape_common , length.out=(r-u))
  }
  if(u==(r-1)){
    omega_hyper_rate_common=mean(1/omega)/ (  mean(1/omega^2)-(mean(1/omega))^2 )
    omega_hyper_shape_common=mean(1/omega)*omega_hyper_rate_common
    ########
    omega_hyper_rate=rep(x =omega_hyper_rate_common , length.out=(u))
    omega_hyper_shape=rep(x =omega_hyper_shape_common , length.out=(u))
    omega0_hyper_shape=1/omega0
    omega0_hyper_rate= 1

  }
  if(   (u<(r-1))*(u>1)   ){
    omega_hyper_rate_common=mean(1/omega)/ (  mean(1/omega^2)-(mean(1/omega))^2 )
    omega0_hyper_rate_common=mean(1/omega0)/(  mean(1/omega0^2)-(mean(1/omega0))^2 )
    omega_hyper_shape_common=mean(1/omega)*omega_hyper_rate_common
    omega0_hyper_shape_common=mean(1/omega0)*omega0_hyper_rate_common
    #print(c(omega_hyper_shape,omega0_hyper_shape))

    ############
    omega_hyper_rate=rep(x =omega_hyper_rate_common , length.out=(u))
    omega_hyper_shape=rep(x =omega_hyper_shape_common , length.out=(u))
    omega0_hyper_rate=rep(x =omega0_hyper_rate_common , length.out=(r-u))
    omega0_hyper_shape=rep(x =omega0_hyper_shape_common , length.out=(r-u))
  }

  if(u==r){
    omega_hyper_rate_common=(mean(1/omega)/ (  mean(1/omega^2)-(mean(1/omega))^2 ))
    omega_hyper_shape_common=abs(mean(1/omega)*omega_hyper_rate_common)+r-p+4
    omega_hyper_rate_common= abs(omega_hyper_rate_common)+5
    ########
    omega_hyper_rate=rep(x =omega_hyper_rate_common , length.out=(u))
    omega_hyper_shape=rep(x =omega_hyper_shape_common , length.out=(u))
    omega0_hyper_rate=numeric(0);omega0_hyper_shape=numeric(0)

  }
  #browser()

  prior_object <- Benvlp_prior(number_of_responses=number_of_responses, number_of_covariates=number_of_covariates, response_envelope_dim=response_envelope_dim,    eta_hyper_e= eta_hyper_e,
                             eta_hyper_C=eta_hyper_C,
                             O_hyper_D=O_hyper_D,
                             O_hyper_G=O_hyper_G,
                             omega_hyper_shape=omega_hyper_shape,
                             omega_hyper_rate=omega_hyper_rate,
                             omega0_hyper_shape=omega0_hyper_shape,
                             omega0_hyper_rate=omega0_hyper_rate, model="Manifold")

  return(prior_object)
}

#u=dim(gamma)[2]
############################################################
############################################################
############################################################
############################################################

validation_of_Benvlp_prior_object<-function(prior_obj){

  if(class(prior_obj)!="Benvlp_prior"){print("Input is not a valid Benvlp_prior object."); return(NULL)}
  if(class(prior_obj)=="Benvlp_prior"){
    Benvlp_prior_attributes<-c("eta_hyper_e", "eta_hyper_C" ,       "O_hyper_D" ,         "O_hyper_G"  ,        "omega_hyper_shape",  "omega_hyper_rate", "omega0_hyper_shape", "omega0_hyper_rate",  "model")
    extra_attributes<-setdiff( names(prior_obj),Benvlp_prior_attributes)
    missing_attributes<-setdiff( Benvlp_prior_attributes,names(prior_obj))

    if(length(extra_attributes)>0){ warning(paste0("The attributes  ",extra_attributes ," present in the input object is not required. The attributes will be ignored while creating the output."))}
    if(length(missing_attributes)>0){ warning(paste0("The required attributes  ",missing_attributes ," are missing in the input object. The required but the missing attributes will be assumed NULL to create the output prior object."))}


    if(is.null(prior_obj$eta_hyper_e)){stop(" Input object is not a valid Benvlp_prior object. prior_obj$eta_hyper_e can not be NULL. ")}
    r=dim(prior_obj$eta_hyper_e)[1];p=dim(prior_obj$eta_hyper_e)[2]
    if(  (length(prior_obj$omega_hyper_shape)+length(prior_obj$omega0_hyper_shape))==0){print("prior_obj is not valid.")

      if(  (length(prior_obj$omega_hyper_rate)+length(prior_obj$omega0_hyper_rate))==0){
        stop("prior_obj is not valid.")
      }
    }

    u_0=max(length(prior_obj$omega_hyper_rate),length(prior_obj$omega_hyper_shape))
    r_minus_u=max(length(prior_obj$omega0_hyper_rate),length(prior_obj$omega0_hyper_shape))
    u=max(u_0, r-r_minus_u)


    prior_object<- Benvlp_prior(number_of_responses=r, number_of_covariates=p, response_envelope_dim=u,    eta_hyper_e= prior_obj$eta_hyper_e,
                              eta_hyper_C=prior_obj$eta_hyper_C,
                              O_hyper_D=prior_obj$O_hyper_D,
                              O_hyper_G=prior_obj$O_hyper_G,
                              omega_hyper_shape=prior_obj$omega_hyper_shape,
                              omega_hyper_rate=prior_obj$omega_hyper_rate,
                              omega0_hyper_shape=prior_obj$omega0_hyper_shape,
                              omega0_hyper_rate=prior_obj$omega0_hyper_rate, model="Manifold")

    return(prior_object)





  }
}




