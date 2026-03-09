#' Generates 2 cross 2 orthogonal matrix from a Generalized Bingham distribution.
#' @param A : First matrix parameter 2 cross 2 symmetric positive definite matrix
#' @param B : second matrix parameter . A 2 cross 2 symmetric positive definite matrix
#' @param discrete: lgical . If TRUE then discrete sampler is used to generate approximate sample. If set discrete= FALSE then initially  a rejection sampler is used to generate exact sample and if there is too many rejection then only the discrete sampler is used.
#' @param Max_rejection: This function tries to use a Rejection sampler approach to get a sample from generalized Bingham distribution
#' In the case of too many rejection the function uses a discrete approximation to generate approximate
#'  sample from the generalized Bingham distribution. Max_rejection is the maximum number of consecutive rejection
#'  before using the discrete sampler.
#'  @param partiton: In the case when the discrete approximation is used partition is the number of partition
#'  considered for the interval 0 to 2*Pi to get approximate disctere distribution to the Generalized bingham distribution.
#' @export
rGbing.O2<-function(A,B,discrete=F,MAX_rejection=10000,partition=400000){

  if(!discrete){ return(Grbing.O2(A,B,MAX_rejection,partition))  }
  if(discrete){ return( discreteGrbing.O2(A,B,N=partition) )  }

}


rvm <- function(mean, kappa) {

}




Grbing.O2<-function(A,B,MAX_rejection=10000,
                    partition=400000, method = 1)  {
   if (method == 1) {

    a <- -(A[1, 1] + B[2, 2])
    b <- -(A[2, 2] + B[1, 1])
    c <- (B[1, 2] + B[2, 1] - A[1, 2] - A[2, 1])

    r <- sqrt((a-b)^2 + c^2)
    alpha <- atan2(y = c, x = a-b)

    theta <- tryCatch(
      {
        th <- BAMBI::minuspi_to_pi(BAMBI::rvm(n = 1, mu = alpha, kappa = r/2)) / 2
        as.numeric(th)
      },
      error = function(e) NA_real_
    )

    # ---- 핵심: theta가 이상하면 method=2로 우회 ----
    if (!is.finite(theta) || length(theta) != 1L) {
      return(Grbing.O2(A, B, MAX_rejection=MAX_rejection, partition=partition, method=2))
    }

    s <- ifelse(runif(2) < 0.5, 1, -1)
    s1 <- s[1]
    s2 <- s[2]

    Z <- s2 * matrix(
      c(cos(theta), s1*sin(theta),
        sin(theta), -s1*cos(theta)),
      ncol = 2, byrow = TRUE
    )

    return(Z)
  }


  else if (method == 2) {
    accept=FALSE
    counter=0
    a=-  (A[1,1]+B[2,2]-A[2,2]-B[1,1])
    b=   B[1,2]+B[2,1]-A[1,2]-A[2,1]

    if(a>0)   {
      counter=1
      #C=A;#A=B;#B=C;

      a=-a;b=-b; }
    reject_count=0

    beta=.573
    gamma=.223
    if(b<0){
      i=1;
      k1=.5*beta(.5-gamma-gamma/2,.5-gamma/2)*beta^2/(a*b)^gamma
      k2=.5*beta(.5-gamma,.5)*beta*exp(-.5*b)/(-a)^gamma
      accept=FALSE  # extra line , just for being sure
      while((!accept)*(reject_count<MAX_rejection)){
        bin=ifelse(k1==Inf,1, rbinom(1,1,k1/(k1+k2))   )

        if(    bin==1  ){
          x=sqrt(rbeta(1,.5-gamma-gamma/2,.5-gamma/2))
          lr=(  a*x^2+b*x*sqrt(1-x^2)  )  -   2*log(beta)+gamma*log(-a*x^2)+gamma*log(-b*x*sqrt(1-x^2))
        }
        if(   bin==0     ){
          x=sqrt(rbeta(1,.5-gamma,.5))
          lr=(   a*x^2-b*x*sqrt(1-x^2)    )  -  log(beta)+.5*b+gamma*log(-a*x^2)
          x=-x # choosing negative square root .
        }


        u=runif(1)
        accept=( log(u) < lr  )
        if(accept){w<-x}
        reject_count=reject_count+1
      } # (end of while !accept)
    } # (end of if b<0)





    if(b>0){
      #i=1
      k1=.5*beta(.5-gamma,.5)*beta*exp(.5*b)/(-a)^gamma
      k2=.5*beta(.5-gamma-gamma/2,.5-gamma/2)*beta^2/(-a*b)^gamma

      accept=FALSE
      while((!accept)*(reject_count<MAX_rejection)){

        bin= ifelse(k1==Inf,1,rbinom(1,1,k1/(k1+k2)))



        if(    bin ==1    ){
          x=sqrt(rbeta(1,.5-gamma,.5))
          lr=(   a*x^2+b*x*sqrt(1-x^2)    )  -  log(beta)-.5*b+gamma*log(-a*x^2)
        }

        if(    bin==0  ){
          x=sqrt(rbeta(1,.5-gamma-gamma/2,.5-gamma/2))
          lr=(  a*x^2-b*x*sqrt(1-x^2)  )  -   2*log(beta)+gamma*log(-a*x^2)+gamma*log(b*x*sqrt(1-x^2))
          x=-x  # choosing negative square root .
        }
        u=runif(1)

        accept=( log(u) < lr  )
        if(accept){w<-x}
        reject_count=reject_count+1
      } # (end of while !accept)


    } # (end of if b>0)

    if(!accept){
      #N=ifelse(N>10000, N, 10000)
      #print("Too many rejection. Approximate sample is generated using Discrete Sampler with default number of partition=400000, if not specified")
      Z=discreteGrbing.O2(A,B,partition)
    }

    if(accept){
      Z=getMfromW(w)
      if(counter==1){Z= cbind(Z[,2],Z[,1])}
    }


    return(Z)
  }
}


getMfromW<-function(w){

  if(w>0){
    x1 <- c((w), sqrt(1 - w^2))  *  (-1)^rbinom(1, 1, 0.5)
    x2 <- (x1[2:1]  *  c(1, -1)  *  (-1)^rbinom(1, 1, 0.5))
  }

  if(w<0){
    x1 <- c(abs(w), -sqrt(1 - w^2))  *  (-1)^rbinom(1, 1, 0.5)  # two elements have to of opposite sign
    x2 <- (x1[2:1]  *  c(1, -1)   *  (-1)^rbinom(1, 1, 0.5))
  }

  return(cbind(x1, x2))
}



########################  discrete approximation of the code #############################
discreteGrbing.O2<-function(A,B,N=400000){
  ######N= number of partition for this portion of code######
  #print("Warning:  Discrete approximation is being used to generate approximate sample.")
  a=-(A[1,1]+B[2,2]-A[2,2]-B[1,1])
  b=B[1,2]+B[2,1]-A[1,2]-A[2,1]
  S=2*pi*(  seq(1:(N))/N-.5/N  )
  #P=a*  (cos(S))^2   + b* cos(S)*sin(S)
  Lprob=(    a*(cos(S))^2   +   b* cos(S)*sin(S)  -  (a)    ) ##
  prob=exp(  Lprob-max(Lprob)  )

  S1=S[which(prob>0)]
  prob1=prob[which(prob>0)]
  x=sample(S1, size=1, replace = FALSE, prob = prob1)

  x1=c(cos(x),sin(x))
  x2=c(sin(x),-cos(x))*(-1)^rbinom(1,1,.5)
  return(cbind(x1,x2))
}




