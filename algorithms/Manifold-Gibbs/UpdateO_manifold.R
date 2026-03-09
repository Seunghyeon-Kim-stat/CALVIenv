
# UpdateO function undates the Orthogonal matrix and gets to the next step of the Markov chian .
# @param O: Current value of the Orthogonal matrix.
# @param G_tilde
# @param G_prior
# @param D_prior
# @param U: matrix containing the response variables
# @param omega: current value of the omega matrix
# @param omega0: current value of the omega0 matrix
# @param allStep :  takes values in c("ALLPOSSIBLEPAIR","RANDOMPAIR","ALLCOLUMN1"). Default is "ALLPOSSIBLEPAIR"
# @param discrete: logical . If set to True discrete approximation is used without using the rejection sampler. If set to FALSE
# then a rejection sampler is used to update 1 pair of column and only in case of too many rejection  the discrete approximation is used.
# @param Max_rejection: This function tries to use a Rejection sampler approach to get a sample from generalized Bingham distribution
# In the case of too many rejection the function uses a discrete approximation to generate approximate
#  sample from the generalized Bingham distribution. Max_rejection is the maximum number of consecutive rejection
#  before using the discrete sampler.
#  @param partiton: In the case when the discrete approximation is used partition is the number of partition
#  considered for the interval 0 to 2*Pi to get approximate disctere distribution to the Generalized bingham distribution.
#  @export
UpdateO <- function(O, G_tilde, G_prior, D_prior, U, omega, omega0, allstep = "ALLPOSSIBLEPAIR", discrete = F, MAX_rejection = 1000, partition = 400000) {
  r <- dim(O)[2]
  allstep <- validator_central_space_sampling_type(central_space_sampling_type = allstep, number_of_response = r)
  if (!is.numeric(allstep)) {
    # if( !is.element(toupper(allstep) ,c("ALLPOSSIBLEPAIR","RANDOMPAIR","ALLCOLUMN1"))){allstep="ALLPOSSIBLEPAIR";print(" The argument for allstep is ignored as it is invalid. The default argument allstep='ALLPOSSIBLEPAIR' is used to run the sampler.")}


    if (toupper(allstep) == "RANDOMPAIR") {
      return(Update2Col(O, G_tilde = G_tilde, G_prior = G_prior, D_prior = D_prior, U, omega, omega0, discrete, MAX_rejection, partition))
    }
    # browser()
    if (toupper(allstep) == "ALLPOSSIBLEPAIR") {
      allPair <- combn((1:r), 2)
      for (i in 1:dim(allPair)[2]) {
        O <- Update2Col(O, G_tilde = G_tilde, G_prior = G_prior, D_prior = D_prior, U, omega, omega0, discrete, MAX_rejection, partition, sel = allPair[, i])
      }
      return(O)
    }

    ################# 3

    if (toupper(allstep) == "ALLCOLUMN1") {
      # r=dim(O)[2]
      # randCol=sample(r)
      for (kk in 1:1) {
        basecol <- r

        randCol <- 1:r
        allPair <- rbind(seq(basecol, basecol, length.out = r - 1), randCol[-basecol])
        # allPair=rbind((r-1):1,seq(basecol,basecol,length.out=r-1))
        # print(allPair)
        for (i in 1:dim(allPair)[2]) {
          O <- Update2Col(O, G_tilde = G_tilde, G_prior = G_prior, D_prior = D_prior, U, omega, omega0, discrete, MAX_rejection, partition, sel = allPair[, i])
        }
      }
      return(O)
    }
    #############



    if (toupper(allstep) == "ALLCOLUMN") {
      # r=dim(O)[2]
      # randCol=sample(r)
      randCol <- 1:r
      allPair <- rbind(randCol[1:(r - 1)], randCol[2:r])
      for (i in 1:dim(allPair)[2]) {
        O <- Update2Col(O, G_tilde = G_tilde, G_prior = G_prior, D_prior = D_prior, U, omega, omega0, discrete, MAX_rejection, partition, sel = allPair[, i])
      }
      for (kk in 1:2) {
        r <- dim(O)[2]
        randCol <- sample(r)
        # randCol=1:r
        allPair <- rbind(randCol[1:(r - 1)], randCol[2:r])
        for (i in 1:dim(allPair)[2]) {
          O <- Update2Col(O, G_tilde = G_tilde, G_prior = G_prior, D_prior = D_prior, U, omega, omega0, discrete, MAX_rejection, partition, sel = allPair[, i])
        }
      }

      return(O)
    }
  }

  if (is.numeric(allstep)) {
    # allstep=as.integer(allstep)
    # r=dim(O)[2]
    allPair <- as.matrix(combn((1:r), 2)[, sample(r * (r - 1) / 2, size = allstep, replace = F)])
    for (i in 1:dim(allPair)[2]) {
      O <- Update2Col(O, G_tilde = G_tilde, G_prior = G_prior, D_prior = D_prior, U, omega, omega0, discrete, MAX_rejection, partition, sel = allPair[, i])
    }
    return(O)
  }
}




#########################################################################
signUpdate <- function(M) {
  return(apply(M, 2, "COLmaxPositive"))
}

COLmaxPositive <- function(a) {
  a <- a * sign(a[which(abs(a) == max(abs(a)))])
  return(a)
}

##########################################################################


Update2Col <- function(O, G_tilde, G_prior, D_prior, U, omega, omega0, discrete = F, MAX_rejection = 1000, partition = 400000, sel = NULL) {
  r <- dim(O)[2]
  u <- length(omega)

  if (is.null(sel)) {
    sel <- sort(sample(1:r, 2, replace = F))
  }
  # print(sel)
  i <- min(sel)
  j <- max(sel)
  # print(paste("sel=",c(i,j)))
  sel <- c(i, j)
  N <- O[, c(i, j)]
  # I=diag(replicate(dim(P_X)[1],1))



  ##### case 1: i in [1,u] and j in [u+1,r] #####
  if ((i <= u) * (j > u)) {
    A <- .5 * t(N) %*% (G_tilde / omega[i] + G_prior / D_prior[i]) %*% N # D_prior[i] new prior
    B <- .5 * t(N) %*% (t(U) %*% U / omega0[j - u] + G_prior / D_prior[j]) %*% N
    # print(" ###### case 1: i in [1,u] #### and ##### j in [u+1,r] ######")
  }


  #### case 2: i,j in [1,u]  #######
  if (j <= u) {
    A <- .5 * t(N) %*% (G_tilde / omega[i] + G_prior / D_prior[i]) %*% N
    B <- .5 * t(N) %*% (G_tilde / omega[j] + G_prior / D_prior[j]) %*% N
    # print(" ##########case 2: i,j in [1,u] ###")
  }


  #### case 3: i,j in [u+1,r]  ######
  if (i > u) {
    A <- .5 * t(N) %*% (t(U) %*% U / omega0[i - u] + G_prior / D_prior[i]) %*% N
    B <- .5 * t(N) %*% (t(U) %*% U / omega0[j - u] + G_prior / D_prior[j]) %*% N
    # print(" ###case 3: i,j in [u+1,r] ############")
  }


  # print(A)
  # print(B)
  Z <- rGbing.O2(A, B, discrete, MAX_rejection, partition) #  provides 2 cross 2 orthogonal matrix

  OijNew <- signUpdate(N %*% Z)
  # Next step is to update the ith col and ,j th col
  O[, sel] <- OijNew
  return(O)
}

#####################
