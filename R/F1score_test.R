
## Takahashi et al. ##
## Computation of test statistics of binary and multi-class F1 scores ##

#install.packages("nleqslv") # install package if required
library(nleqslv)

f1scores_test <- function(test1,test2,actual, binary1=c("A")){
  ## This function computes test statistics
  ## for binaryF1, microF1, macroF1, and macroF1* scores.
  
  ## test1 is a numeric vector of data for predicted condition of test 1.
  ## test2 is a numeric vector of data for predicted condition of test 2.
  ## actual is a numeric vector of data for actual condition.
  ## binary1 is conditions included in class 1 when the data is dichotomized to calculate binaryF1.
  
  ## biF1 is binary F1
  ## miF1 is micro F1
  ## maF1 is macro F1
  ## maF2 is macro F1* (Sokolova and Lapalme)
  
  
  ## ###### ##
  ## Set up ##
  ## ###### ##
  
  tab=ftable(test1,test2,actual)
  tab=as.data.frame.matrix(tab)
  mat=list()
  for(i in 1:ncol(tab)){
    mat[[i]]=t(matrix(tab[,i],nrow=ncol(tab)))
  }
  
  r <- length(mat) # Number of classes
  N <- sum(sapply(mat,sum)) ## Total sample size
  p <- lapply(mat,function(x){x/N}) ## Probabilities
  
  
  matv = NA
  for(i in 1:r){
    matvb = as.vector(t(mat[[i]]))
    matv = c(matv,matvb)
  }
  matv = matv[-1]
  matvp = matv/N
  
  
  pi.j <- p.ij <- matrix(numeric(r*r),nrow=r)
  pi.. <- p.i. <- piii <- 0
  for(i in 1:r){
    for(j in 1:r){
      pi.j[i,j]=sum(p[[j]][i,])
      p.ij[i,j]=sum(p[[j]][,i])
    }
    pi..[i]=sum(pi.j[i,])
    p.i.[i]=sum(p.ij[i,])
    piii[i]=p[[i]][i,i]
  }
  p..i <- apply(pi.j,2,sum)
  pi.i <- diag(pi.j)
  p.ii <- diag(p.ij)
  
  
  test1_bi = ifelse(test1 %in% binary1,1,2)
  test2_bi = ifelse(test2 %in% binary1,1,2)
  actual_bi = ifelse(actual %in% binary1,1,2)
  tab_bi=ftable(test1_bi,test2_bi,actual_bi)
  tab_bi=as.data.frame.matrix(tab_bi)
  mat_bi=list()
  for(i in 1:ncol(tab_bi)){
    mat_bi[[i]]=t(matrix(tab_bi[,i],nrow=ncol(tab_bi)))
  }
  p_bi <- lapply(mat_bi,function(x){x/N}) ## Probabilities
  
  
  p111 <- p_bi[[1]][1,1]
  p121 <- p_bi[[1]][1,2]
  p211 <- p_bi[[1]][2,1]
  p221 <- p_bi[[1]][2,2]
  p112 <- p_bi[[2]][1,1]
  p122 <- p_bi[[2]][1,2]
  p212 <- p_bi[[2]][2,1]
  p222 <- p_bi[[2]][2,2]
  
  p1.1 <- p111+p121
  p1.2 <- p112+p122
  p2.1 <- p211+p221
  p.11 <- p111+p211
  p.12 <- p112+p212
  p.21 <- p121+p221
  p1.. <- p1.1+p1.2
  p.1. <- p.11+p.12
  p..1 <- p1.1+p2.1
  
  
  
  ## Binary F1 Score ---------------
  
  
  ## --------------- ##
  ## Point estimates ##
  ## --------------- ##
  
  biP1 <- p1.1/p1..
  biR1 <- p1.1/p..1
  biF1 <- 2*biP1*biR1/(biP1+biR1) ## binary F1 for test 1
  
  biP2 <- p.11/p.1.
  biR2 <- p.11/p..1
  biF2 <- 2*biP2*biR2/(biP2+biR2) ## binary F1 for test 2
  biF.d <- biF1 - biF2  ## Difference of binary F1 between test 1 and test 2
  
  ## --------------------------------- ##
  ## Test statistics & p values (Wald) ##
  ## --------------------------------- ##
  
  biF1.v <- 1/(p1..+p..1)^2*(p1.1*(2*(1-biF1))^2 + (p1.2+p2.1)*biF1^2)
  biF2.v <- 1/(p.1.+p..1)^2*(p.11*(2*(1-biF2))^2 + (p.12+p.21)*biF2^2)
  biF.cov <- 1/(p.1.+p..1)/(p.1.+p..1)*(2^2*p111*(1-biF1)*(1-biF2) - 2*p121*(1-biF1)*biF2
                                        -2*p211*biF1*(1-biF2) + (p221+p112)*biF1*biF2)
  biF.v <- (biF1.v+biF2.v-2*biF.cov)/N ## Variance of binary F1 (Wald)
  
  biF.T <- (biF.d)^2 / biF.v
  biF.p <- pchisq(biF.T, 1, lower.tail=F) ## p-value of binary F1 (Wald)
  
  ## ---------------------------------- ##
  ## Test statistics & p values (Score) ##
  ## ---------------------------------- ##
  
  f1binary_s <- function(x){
    y <- 0
    
    x1.1=x[1]+x[2]
    x1.2=x[5]+x[6]
    x2.1=x[3]+x[4]
    x.11=x[1]+x[3]
    x.12=x[5]+x[7]
    x.21=x[2]+x[4]
    x.22=x[6]+x[8]
    x1..=x1.1+x1.2
    x.1.=x.11+x.12
    x..1=x1.1+x2.1
    biF1s=2*x1.1/(x1..+x..1)
    biF2s=2*x.11/(x.1.+x..1)
    
    y[1] <- mat_bi[[1]][1,1]-N*x[1]-2*x[9]*x[1]*((1-biF1s)/(x1..+x..1)-(1-biF2s)/(x.1.+x..1))
    y[2] <- mat_bi[[1]][1,2]-N*x[2]-2*x[9]*x[2]*((1-biF1s)/(x1..+x..1)+biF2s/(2*(x.1.+x..1)))
    y[3] <- mat_bi[[1]][2,1]-N*x[3]+2*x[9]*x[3]*(biF1s/(2*(x1..+x..1))+(1-biF2s)/(x.1.+x..1))
    y[4] <- mat_bi[[1]][2,2]-N*x[4]+2*x[9]*x[4]*(biF1s/(2*(x1..+x..1))-biF2s/(2*(x.1.+x..1)))
    y[5] <- mat_bi[[2]][1,1]-N*x[5]+2*x[9]*x[5]*(biF1s/(2*(x1..+x..1))-biF2s/(2*(x.1.+x..1)))
    y[6] <- mat_bi[[2]][1,2]-N*x[6]+2*x[9]*x[6]*(biF1s/(2*(x1..+x..1)))
    y[7] <- mat_bi[[2]][2,1]-N*x[7]-2*x[9]*x[7]*(biF2s/(2*(x.1.+x..1)))
    y[8] <- mat_bi[[2]][2,2]-N*x[8]
    y[9] <- x1.1/(x1..+x..1) - x.11/(x.1.+x..1)
    
    return(y)
  }
  
  pv = c(p111,p121,p211,p221,p112,p122,p212,p222)
  biF.nr=nleqslv(c(pv,1),
                 f1binary_s,jacobian=TRUE,method="Newton",control=list(maxit=1000))
  phat=biF.nr$x
  
  ps111 <- phat[1]
  ps121 <- phat[2]
  ps211 <- phat[3]
  ps221 <- phat[4]
  ps112 <- phat[5]
  ps122 <- phat[6]
  ps212 <- phat[7]
  
  ps1.1 <- ps111+ps121
  ps1.2 <- ps112+ps122
  ps2.1 <- ps211+ps221
  ps.11 <- ps111+ps211
  ps.12 <- ps112+ps212
  ps.21 <- ps121+ps221
  ps1.. <- ps1.1+ps1.2
  ps.1. <- ps.11+ps.12
  ps..1 <- ps1.1+ps2.1
  
  biP1s <- ps1.1/ps1..
  biR1s <- ps1.1/ps..1
  biF1s <- 2*biP1s*biR1s/(biP1s+biR1s)
  
  biP2s <- ps.11/ps.1.
  biR2s <- ps.11/ps..1
  biF2s <- 2*biP2s*biR2s/(biP2s+biR2s)
  
  ## --------------------------------- ##
  ## Test statistics & p values (Wald) ##
  ## --------------------------------- ##
  
  biF1s.v <- 1/(ps1..+ps..1)^2*(ps1.1*(2*(1-biF1s))^2 + (ps1.2+ps2.1)*biF1s^2)
  biF2s.v <- 1/(ps.1.+ps..1)^2*(ps.11*(2*(1-biF2s))^2 + (ps.12+ps.21)*biF2s^2)
  biFs.cov <- 1/(ps.1.+ps..1)/(ps.1.+ps..1)*(2^2*ps111*(1-biF1s)*(1-biF2s) - 2*ps121*(1-biF1s)*biF2s
                                             -2*ps211*biF1s*(1-biF2s) + (ps221+ps112)*biF1s*biF2s)
  biFs.v <- (biF1s.v+biF2s.v-2*biFs.cov)/N ## Variance of binary F1 (Score)
  
  biFs.T <- (biF.d)^2 / biFs.v
  biFs.p <- pchisq(biFs.T, 1, lower.tail=F) ## p-value of binary F1 (Score)
  
  
  
  ## Micro F1 Score ---------------
  
  
  ## --------------- ##
  ## Point estimates ##
  ## --------------- ##
  
  miF1 <- sum(pi.i) ## Micro F1 for test 1
  miF2 <- sum(p.ii) ## Micro F1 for test 2
  miF.d <- miF1 - miF2  ## Difference of Micro F1 between test 1 and test 2
  
  ## --------------------------------- ##
  ## Test statistics & p values (Wald) ##
  ## --------------------------------- ##
  
  miF1.v <- miF1*(1-miF1)
  miF2.v <- miF2*(1-miF2)
  miF.cov <- sum(piii)-miF1*miF2
  miF.v <- (miF1.v+miF2.v-2*miF.cov)/N  ## Variance of Micro F1 (Wald)
  
  miF.T <- (miF.d)^2 / miF.v
  miF.p <- pchisq(miF.T, 1, lower.tail=F)  ## p-value of Micro F1 (Wald)
  
  ## ---------------------------------- ##
  ## Test statistics & p values (Score) ##
  ## ---------------------------------- ##
  
  f1micro_s <- function(x){
    y <- 0
    xi.i <- x.ii <- 0
    
    for(i in 1:r){
      y[r^2*(i-1)+r*(i-1)+i]=mat[[i]][i,i]-N*x[r^2*(i-1)+r*(i-1)+i] #piii
      for(j in c(1:r)[-i]){
        y[r^2*(i-1)+r*(i-1)+j]=mat[[i]][i,j]-N*x[r^2*(i-1)+r*(i-1)+j]-
          x[r^3+1]*x[r^2*(i-1)+r*(i-1)+j]  #piji
        y[r^2*(i-1)+r*(j-1)+i]=mat[[i]][j,i]-N*x[r^2*(i-1)+r*(j-1)+i]+
          x[r^3+1]*x[r^2*(i-1)+r*(j-1)+i] #pjii
        for(k in c(1:r)[-i]){
          y[r^2*(i-1)+r*(j-1)+k]=mat[[i]][j,k]-N*x[r^2*(i-1)+r*(j-1)+k] #pjki
        }}
      xi.i[i]=sum(x[(r^2*(i-1)+r*(i-1)+1):(r^2*(i-1)+r*(i-1)+r)])
      x.ii[i]=sum(x[seq(r^2*(i-1)+i,r^2*(i-1)+r*(r-1)+i,by=r)])
    }
    y[r^3+1] <- sum(xi.i)-sum(x.ii)
    
    return(y)
  }
  
  miF.nr=nleqslv(c(matvp,1),
                 f1micro_s,jacobian=TRUE,method="Newton",control=list(maxit=1000))
  phat=miF.nr$x
  
  ps <- list(0)
  for(i in 1:r){
    ps[[i]] <- t(matrix(phat[(r^2*(i-1)+1):(r^2*i)], ncol=r, nrow=r))
  }
  
  psi.j <- ps.ij <- matrix(numeric(r*r),nrow=r)
  psi.. <- ps.i. <- psiii <- 0
  for(i in 1:r){
    for(j in 1:r){
      psi.j[i,j]=sum(ps[[j]][i,])
      ps.ij[i,j]=sum(ps[[j]][,i])
    }
    psi..[i]=sum(psi.j[i,])
    ps.i.[i]=sum(ps.ij[i,])
    psiii[i]=sum(ps[[i]][i,i])
  }
  ps..i <- apply(psi.j,2,sum)
  psi.i <- diag(psi.j)
  ps.ii <- diag(ps.ij)
  
  miF1s <- sum(psi.i)
  miF2s <- sum(ps.ii)
  
  miF1s.v <- miF1s*(1-miF1s)
  miF2s.v <- miF2s*(1-miF2s)
  miFs.cov <- sum(psiii)-miF1s*miF2s
  miFs.v <- (miF1s.v+miF2s.v-2*miFs.cov)/N ## Variance of Micro F1 (Score)
  
  miFs.T <- (miF.d)^2 / miFs.v
  miFs.p <- pchisq(miFs.T, 1, lower.tail=F) ## p-value of Micro F1 (Score)
  
  
  ## Macro F1 Score ---------------
  
  
  ## --------------- ##
  ## Point estimates ##
  ## --------------- ##
  
  P1 <- pi.i/pi..
  R1 <- pi.i/p..i
  F1 <- 2*P1*R1/(P1+R1)
  maF1 <- sum(F1)/r ## Macro F1 for test 1
  
  P2 <- p.ii/p.i.
  R2 <- p.ii/p..i
  F2 <- 2*P2*R2/(P2+R2)
  maF2 <- sum(F2)/r ## Macro F1 for test 2
  maF.d <- maF1 - maF2  ## Difference of Macro F1 between test 1 and test 2
  
  ## --------------------------------- ##
  ## Test statistics & p values (Wald) ##
  ## --------------------------------- ##
  
  v1 <- v2 <- matrix(numeric(r*r),nrow=r)
  v3 <- list(0)
  for(i in 1:r){
    v1[i,i] <- pi.i[i]*(2*(1-F1[i])/(pi..[i]+p..i[i]))^2
    v2[i,i] <- p.ii[i]*(2*(1-F2[i])/(p.i.[i]+p..i[i]))^2
    v3[[i]] <- matrix(numeric(r*r),nrow=r)
    v3[[i]][i,i] <- p[[i]][i,i]*2^2*(1-F1[i])*(1-F2[i])/(pi..[i]+p..i[i])/(p.i.[i]+p..i[i])
    for(j in c(1:r)[-i]){
      v1[i,j] <- pi.j[i,j]*(F1[i]/(pi..[i]+p..i[i]) + F1[j]/(pi..[j]+p..i[j]))^2
      v2[i,j] <- p.ij[i,j]*(F2[i]/(p.i.[i]+p..i[i]) + F2[j]/(p.i.[j]+p..i[j]))^2
      v3[[i]][i,j] <- -p[[i]][i,j]*2*(1-F1[i])/(pi..[i]+p..i[i])*
        (F2[i]/(p.i.[i]+p..i[i]) + F2[j]/(p.i.[j]+p..i[j]))
      v3[[i]][j,i] <- -p[[i]][j,i]*2*(1-F2[i])/(p.i.[i]+p..i[i])*
        (F1[i]/(pi..[i]+p..i[i]) + F1[j]/(pi..[j]+p..i[j]))
      for(k in c(1:r)[-i]){
        v3[[i]][j,k] <- p[[i]][j,k]*(F1[i]/(pi..[i]+p..i[i]) + F1[j]/(pi..[j]+p..i[j]))*
          (F2[i]/(p.i.[i]+p..i[i]) + F2[k]/(p.i.[k]+p..i[k]))
      }
    }
  }
  
  maF1.v <- 1/r^2*sum(v1)
  maF2.v <- 1/r^2*sum(v2)
  maF.cov <- 1/r^2*sum(sapply(v3,sum))
  maF.v <- (maF1.v+maF2.v-2*maF.cov)/N  ## Variance of Macro F1 (Wald)
  
  maF.T <- (maF.d)^2 / maF.v
  maF.p <- pchisq(maF.T, 1, lower.tail=F)  ## p-value of Macro F1 (Wald)
  
  ## ---------------------------------- ##
  ## Test statistics & p values (Score) ##
  ## ---------------------------------- ##
  
  f1macro_s <- function(x){
    y <- 0
    xj.i <- x.ji <- xi.i <- x.ii <- xi.. <- x.i. <- x..i <- 0
    P1s <- R1s <- F1s <- P2s <- R2s <- F2s <- 0
    
    for(i in 1:r){
      for(j in 1:r){
        xj.i[r*(i-1)+j]=sum(x[(r^2*(i-1)+r*(j-1)+1):(r^2*(i-1)+r*j)])
        x.ji[r*(i-1)+j]=sum(x[seq(r^2*(i-1)+j,r^2*(i-1)+r*(r-1)+j,by=r)])
      }
      xi.i[i]=xj.i[r*(i-1)+i]
      x.ii[i]=x.ji[r*(i-1)+i]
    }
    for(i in 1:r){
      xi..[i]=sum(xj.i[seq(i,r*(r-1)+i,by=r)])
      x.i.[i]=sum(x.ji[seq(i,r*(r-1)+i,by=r)])
      x..i[i]=sum(xj.i[(r*(i-1)+1):(r*i)])
      
      F1s[i] <- 2*xi.i[i]/(xi..[i]+x..i[i])
      F2s[i] <- 2*x.ii[i]/(x.i.[i]+x..i[i])
    }
    
    for(i in 1:r){
      y[r^2*(i-1)+r*(i-1)+i]=mat[[i]][i,i]-N*x[r^2*(i-1)+r*(i-1)+i]-
        2/r*x[r^3+1]*x[r^2*(i-1)+r*(i-1)+i]*
        ((1-F1s[i])/(xi..[i]+x..i[i])-(1-F2s[i])/(x.i.[i]+x..i[i])) #piii
      for(j in c(1:r)[-i]){
        y[r^2*(i-1)+r*(i-1)+j]=mat[[i]][i,j]-N*x[r^2*(i-1)+r*(i-1)+j]-
          1/r*x[r^3+1]*x[r^2*(i-1)+r*(i-1)+j]*
          (2*(1-F1s[i])/(xi..[i]+x..i[i])+
             F2s[i]/(x.i.[i]+x..i[i])+F2s[j]/(x.i.[j]+x..i[j]))  #piji
        y[r^2*(i-1)+r*(j-1)+i]=mat[[i]][j,i]-N*x[r^2*(i-1)+r*(j-1)+i]+
          1/r*x[r^3+1]*x[r^2*(i-1)+r*(j-1)+i]*
          (F1s[i]/(xi..[i]+x..i[i])+F1s[j]/(xi..[j]+x..i[j])+
             2*(1-F2s[i])/(x.i.[i]+x..i[i]))  #pjii
        for(k in c(1:r)[-i]){
          y[r^2*(i-1)+r*(j-1)+k]=mat[[i]][j,k]-N*x[r^2*(i-1)+r*(j-1)+k]+
            1/r*x[r^3+1]*x[r^2*(i-1)+r*(j-1)+k]*
            (F1s[i]/(xi..[i]+x..i[i])+F1s[j]/(xi..[j]+x..i[j])-
               F2s[i]/(x.i.[i]+x..i[i])-F2s[k]/(x.i.[k]+x..i[k])) #pjki
        }}}
    y[r^3+1] <- sum(F1s)-sum(F2s)
    
    return(y)
  }
  
  
  maF.nr=nleqslv(c(matvp,1),
                 f1macro_s,jacobian=TRUE,method="Newton",control=list(maxit=1000))
  phat=maF.nr$x
  
  ps <- list(0)
  for(i in 1:r){
    ps[[i]] <- t(matrix(phat[(r*r*(i-1)+1):(r*r*i)], ncol=r, nrow=r))
  }
  
  psi.j <- ps.ij <- matrix(numeric(r*r),nrow=r)
  psi.. <- ps.i. <- psiii <- 0
  for(i in 1:r){
    for(j in 1:r){
      psi.j[i,j]=sum(ps[[j]][i,])
      ps.ij[i,j]=sum(ps[[j]][,i])
    }
    psi..[i]=sum(psi.j[i,])
    ps.i.[i]=sum(ps.ij[i,])
    psiii[i]=ps[[i]][i,i]
  }
  ps..i <- apply(psi.j,2,sum)
  psi.i <- diag(psi.j)
  ps.ii <- diag(ps.ij)
  
  P1s <- psi.i/psi..
  R1s <- psi.i/ps..i
  F1s <- 2*P1s*R1s/(P1s+R1s)
  maF1s <- sum(F1s)/r
  
  P2s <- ps.ii/ps.i.
  R2s <- ps.ii/ps..i
  F2s <- 2*P2s*R2s/(P2s+R2s)
  maF2s <- sum(F2s)/r
  
  v1s <- v2s <- matrix(numeric(r*r),nrow=r)
  v3s <- list(0)
  for(i in 1:r){
    v1s[i,i] <- psi.i[i]*(2*(1-F1s[i])/(psi..[i]+ps..i[i]))^2
    v2s[i,i] <- ps.ii[i]*(2*(1-F2s[i])/(ps.i.[i]+ps..i[i]))^2
    v3s[[i]] <- matrix(numeric(r*r),nrow=r)
    v3s[[i]][i,i] <- ps[[i]][i,i]*2^2*(1-F1s[i])*(1-F2s[i])/(psi..[i]+ps..i[i])/(ps.i.[i]+ps..i[i])
    for(j in c(1:r)[-i]){
      v1s[i,j] <- psi.j[i,j]*(F1s[i]/(psi..[i]+ps..i[i]) + F1s[j]/(psi..[j]+ps..i[j]))^2
      v2s[i,j] <- ps.ij[i,j]*(F2s[i]/(ps.i.[i]+ps..i[i]) + F2s[j]/(ps.i.[j]+ps..i[j]))^2
      v3s[[i]][i,j] <- -ps[[i]][i,j]*2*(1-F1s[i])/(psi..[i]+ps..i[i])*
        (F2s[i]/(ps.i.[i]+ps..i[i]) + F2s[j]/(ps.i.[j]+ps..i[j]))
      v3s[[i]][j,i] <- -ps[[i]][j,i]*2*(1-F2s[i])/(ps.i.[i]+ps..i[i])*
        (F1s[i]/(psi..[i]+ps..i[i]) + F1s[j]/(psi..[j]+ps..i[j]))
      for(k in c(1:r)[-i]){
        v3s[[i]][j,k] <- ps[[i]][j,k]*(F1s[i]/(psi..[i]+ps..i[i]) + F1s[j]/(psi..[j]+ps..i[j]))*
          (F2s[i]/(ps.i.[i]+ps..i[i]) + F2s[k]/(ps.i.[k]+ps..i[k]))
      }
    }
  }
  
  maF1s.v <- 1/r^2*sum(v1s)
  maF2s.v <- 1/r^2*sum(v2s)
  maFs.cov <- 1/r^2*sum(sapply(v3s,sum))
  maFs.v <- (maF1s.v+maF2s.v-2*maFs.cov)/N ## Variance of Macro F1 (Score)
  
  maFs.T <- (maF.d)^2 / maFs.v
  maFs.p <- pchisq(maFs.T, 1, lower.tail=F) ## p-value of Macro F1 (Score)
  
  
  ## Macro F1* Score ---------------
  
  
  ## --------------- ##
  ## Point estimates ##
  ## --------------- ##
  
  maP1 <- sum(P1)/r
  maR1 <- sum(R1)/r
  maF1b <- 2*maP1*maR1/(maP1+maR1)  ## Macro F1* for test 1
  
  maP2 <- sum(P2)/r
  maR2 <- sum(R2)/r
  maF2b <- 2*maP2*maR2/(maP2+maR2)  ## Macro F1* for test 2
  maFb.d <- maF1b - maF2b  ## Difference of Macro F1* between test 1 and test 2
  
  ## --------------------------------- ##
  ## Test statistics & p values (Wald) ##
  ## --------------------------------- ##
  
  v1 <- v2 <- matrix(numeric(r*r),nrow=r)
  v3 <- list(0)
  for(i in 1:r){
    v1[i,i] <- pi.i[i]*((pi..[i]-pi.i[i])*maR1^2/pi..[i]^2 + (p..i[i]-pi.i[i])*maP1^2/p..i[i]^2)^2
    v2[i,i] <- p.ii[i]*((p.i.[i]-p.ii[i])*maR2^2/p.i.[i]^2 + (p..i[i]-p.ii[i])*maP2^2/p..i[i]^2)^2
    v3[[i]] <- matrix(numeric(r*r),nrow=r)
    v3[[i]][i,i] <- p[[i]][i,i]*((pi..[i]-pi.i[i])*maR1^2/pi..[i]^2 +
                                   (p..i[i]-pi.i[i])*maP1^2/p..i[i]^2)*
      ((p.i.[i]-p.ii[i])*maR2^2/p.i.[i]^2 +
         (p..i[i]-p.ii[i])*maP2^2/p..i[i]^2)
    for(j in c(1:r)[-i]){
      v1[i,j] <- pi.j[i,j]*(pi.j[i,i]*maR1^2/pi..[i]^2 + pi.j[j,j]*maP1^2/p..i[j]^2)^2
      v2[i,j] <- p.ij[i,j]*(p.ij[i,i]*maR2^2/p.i.[i]^2 + p.ij[j,j]*maP2^2/p..i[j]^2)^2
      v3[[i]][i,j] <- -p[[i]][i,j]*((pi..[i]-pi.i[i])*maR1^2/pi..[i]^2 +
                                      (p..i[i]-pi.i[i])*maP1^2/p..i[i]^2)*
        (p.ij[j,j]*maR2^2/p.i.[j]^2 + p.ij[i,i]*maP2^2/p..i[i]^2)
      v3[[i]][j,i] <- -p[[i]][j,i]*((p.i.[i]-p.ii[i])*maR2^2/p.i.[i]^2 +
                                      (p..i[i]-p.ii[i])*maP2^2/p..i[i]^2)*
        (pi.j[j,j]*maR1^2/pi..[j]^2 + pi.j[i,i]*maP1^2/p..i[i]^2)
      for(k in c(1:r)[-i]){
        v3[[i]][j,k] <- p[[i]][j,k]*(pi.j[j,j]*maR1^2/pi..[j]^2 + pi.j[i,i]*maP1^2/p..i[i]^2)*
          (p.ij[k,k]*maR2^2/p.i.[k]^2 + p.ij[i,i]*maP2^2/p..i[i]^2)
      }
    }
  }
  
  maF1b.v <- 2^2/r^2/(maP1+maR1)^4*sum(v1)
  maF2b.v <- 2^2/r^2/(maP2+maR2)^4*sum(v2)
  maFb.cov <- 2^2/r^2/(maP1+maR1)^2/(maP2+maR2)^2*sum(sapply(v3,sum))
  maFb.v <- (maF1b.v+maF2b.v-2*maFb.cov)/N ## Variance of Macro F1* (Wald)
  
  maFb.T <- (maFb.d)^2 / maFb.v
  maFb.p <- pchisq(maFb.T, 1, lower.tail=F) ## p-value of Macro F1* (Wald)
  
  ## ---------------------------------- ##
  ## Test statistics & p values (Score) ##
  ## ---------------------------------- ##
  
  f1macrob_s <- function(x){
    y <- 0
    xj.i <- x.ji <- xi.i <- x.ii <- xi.. <- x.i. <- x..i <- 0
    P1s <- R1s <- P2s <- R2s <- 0
    
    for(i in 1:r){
      for(j in 1:r){
        xj.i[r*(i-1)+j]=sum(x[(r^2*(i-1)+r*(j-1)+1):(r^2*(i-1)+r*j)])
        x.ji[r*(i-1)+j]=sum(x[seq(r^2*(i-1)+j,r^2*(i-1)+r*(r-1)+j,by=r)])
      }
      xi.i[i]=xj.i[r*(i-1)+i]
      x.ii[i]=x.ji[r*(i-1)+i]
    }
    for(i in 1:r){
      xi..[i]=sum(xj.i[seq(i,r*(r-1)+i,by=r)])
      x.i.[i]=sum(x.ji[seq(i,r*(r-1)+i,by=r)])
      x..i[i]=sum(xj.i[(r*(i-1)+1):(r*i)])
      
      P1s[i] <- xi.i[i]/xi..[i]
      R1s[i] <- xi.i[i]/x..i[i]
      P2s[i] <- x.ii[i]/x.i.[i]
      R2s[i] <- x.ii[i]/x..i[i]
    }
    
    maP1s <- sum(P1s)/r
    maR1s <- sum(R1s)/r
    maP2s <- sum(P2s)/r
    maR2s <- sum(R2s)/r
    
    for(i in 1:r){
      y[r^2*(i-1)+r*(i-1)+i]=mat[[i]][i,i]-N*x[r^2*(i-1)+r*(i-1)+i]-
        2/r*x[r^3+1]*x[r^2*(i-1)+r*(i-1)+i]*
        (((xi..[i]-xi.i[i])/xi..[i]^2*maR1s^2+(x..i[i]-xi.i[i])/x..i[i]^2*maP1s^2)/(maP1s+maR1s)^2
         -((x.i.[i]-x.ii[i])/x.i.[i]^2*maR2s^2-(x..i[i]-x.ii[i])/x..i[i]^2*maP2s^2)/(maP2s+maR2s)^2) #piii
      for(j in c(1:r)[-i]){
        y[r^2*(i-1)+r*(i-1)+j]=mat[[i]][i,j]-N*x[r^2*(i-1)+r*(i-1)+j]-
          2/r*x[r^3+1]*x[r^2*(i-1)+r*(i-1)+j]*
          (((xi..[i]-xi.i[i])/xi..[i]^2*maR1s^2+(x..i[i]-xi.i[i])/x..i[i]^2*maP1s^2)/(maP1s+maR1s)^2
           +(x.ii[j]/x.i.[j]^2*maR2s^2+x.ii[i]/x..i[i]^2*maP2s^2)/(maP2s+maR2s)^2)  #piji
        y[r^2*(i-1)+r*(j-1)+i]=mat[[i]][j,i]-N*x[r^2*(i-1)+r*(j-1)+i]+
          2/r*x[r^3+1]*x[r^2*(i-1)+r*(j-1)+i]*
          ((xi.i[j]/xi..[j]^2*maR1s^2+xi.i[i]/x..i[i]^2*maP1s^2)/(maP1s+maR1s)^2
           +((x.i.[i]-x.ii[i])/x.i.[i]^2*maR2s^2+(x..i[i]-x.ii[i])/x..i[i]^2*maP2s^2)/(maP2s+maR2s)^2)  #pjii
        for(k in c(1:r)[-i]){
          y[r^2*(i-1)+r*(j-1)+k]=mat[[i]][j,k]-N*x[r^2*(i-1)+r*(j-1)+k]+
            2/r*x[r^3+1]*x[r^2*(i-1)+r*(j-1)+k]*
            ((xi.i[j]/xi..[j]^2*maR1s^2+xi.i[i]/x..i[i]^2*maP1s^2)/(maP1s+maR1s)^2
             -(x.ii[k]/x.i.[k]^2*maR2s^2-x.ii[i]/x..i[i]^2*maP2s^2)/(maP2s+maR2s)^2) #pjki
        }}}
    y[r^3+1] <- maP1s*maR1s/(maP1s+maR1s) - maP2s*maR2s/(maP2s+maR2s)
    
    return(y)
  }
  
  
  maFb.nr=nleqslv(c(matvp,1),
                  f1macrob_s,jacobian=TRUE,method="Newton",control=list(maxit=1000))
  phat=maFb.nr$x
  
  ps <- list(0)
  for(i in 1:r){
    ps[[i]] <- t(matrix(phat[(r*r*(i-1)+1):(r*r*i)], ncol=r, nrow=r))
  }
  
  psi.j <- ps.ij <- matrix(numeric(r*r),nrow=r)
  psi.. <- ps.i. <- psiii <- 0
  for(i in 1:r){
    for(j in 1:r){
      psi.j[i,j]=sum(ps[[j]][i,])
      ps.ij[i,j]=sum(ps[[j]][,i])
    }
    psi..[i]=sum(psi.j[i,])
    ps.i.[i]=sum(ps.ij[i,])
    psiii[i]=sum(ps[[i]][i,i])
  }
  ps..i <- apply(psi.j,2,sum)
  psi.i <- diag(psi.j)
  ps.ii <- diag(ps.ij)
  
  P1s <- psi.i/psi..
  R1s <- psi.i/ps..i
  P2s <- ps.ii/ps.i.
  R2s <- ps.ii/ps..i
  
  maP1s <- sum(P1s)/r
  maR1s <- sum(R1s)/r
  maF1bs <- 2*maP1s*maR1s/(maP1s+maR1s)
  
  maP2s <- sum(P2s)/r
  maR2s <- sum(R2s)/r
  maF2bs <- 2*maP2s*maR2s/(maP2s+maR2s)
  
  
  v1s <- v2s <- matrix(numeric(r*r),nrow=r)
  v3s <- list(0)
  for(i in 1:r){
    v1s[i,i] <- psi.i[i]*((psi..[i]-psi.i[i])*maR1s^2/psi..[i]^2 + (ps..i[i]-psi.i[i])*
                            maP1s^2/ps..i[i]^2)^2
    v2s[i,i] <- ps.ii[i]*((ps.i.[i]-ps.ii[i])*maR2s^2/ps.i.[i]^2 + (ps..i[i]-ps.ii[i])*
                            maP2s^2/ps..i[i]^2)^2
    v3s[[i]] <- matrix(numeric(r*r),nrow=r)
    v3s[[i]][i,i] <- ps[[i]][i,i]*((psi..[i]-psi.i[i])*maR1s^2/psi..[i]^2 +
                                     (ps..i[i]-psi.i[i])*maP1s^2/ps..i[i]^2)*
      ((ps.i.[i]-ps.ii[i])*maR2s^2/ps.i.[i]^2 +
         (ps..i[i]-ps.ii[i])*maP2s^2/ps..i[i]^2)
    for(j in c(1:r)[-i]){
      v1s[i,j] <- psi.j[i,j]*(psi.j[i,i]*maR1s^2/psi..[i]^2 + psi.j[j,j]*maP1s^2/ps..i[j]^2)^2
      v2s[i,j] <- ps.ij[i,j]*(ps.ij[i,i]*maR2s^2/ps.i.[i]^2 + ps.ij[j,j]*maP2s^2/ps..i[j]^2)^2
      v3s[[i]][i,j] <- -ps[[i]][i,j]*((psi..[i]-psi.i[i])*maR1s^2/psi..[i]^2 +
                                        (ps..i[i]-psi.i[i])*maP1s^2/ps..i[i]^2)*
        (ps.ij[j,j]*maR2s^2/ps.i.[j]^2 + ps.ij[i,i]*maP2s^2/ps..i[i]^2)
      v3s[[i]][j,i] <- -ps[[i]][j,i]*((ps.i.[i]-ps.ii[i])*maR2s^2/ps.i.[i]^2 +
                                        (ps..i[i]-ps.ii[i])*maP2s^2/ps..i[i]^2)*
        (psi.j[j,j]*maR1s^2/psi..[j]^2 + psi.j[i,i]*maP1s^2/ps..i[i]^2)
      for(k in c(1:r)[-i]){
        v3s[[i]][j,k] <- ps[[i]][j,k]*(psi.j[j,j]*maR1s^2/psi..[j]^2 + psi.j[i,i]*maP1s^2/ps..i[i]^2)*
          (ps.ij[k,k]*maR2s^2/ps.i.[k]^2 + ps.ij[i,i]*maP2s^2/ps..i[i]^2)
      }
    }
  }
  
  maF1bs.v <- 2^2/r^2/(maP1s+maR1s)^4*sum(v1s)
  maF2bs.v <- 2^2/r^2/(maP2s+maR2s)^4*sum(v2s)
  maFbs.cov <- 2^2/r^2/(maP1s+maR1s)^2/(maP2s+maR2s)^2*sum(sapply(v3s,sum))
  maFbs.v <- (maF1bs.v+maF2bs.v-2*maFbs.cov)/N ## Variance of Macro F1* (Score)
  
  maFbs.T <- (maFb.d)^2 / maFbs.v
  maFbs.p <- pchisq(maFbs.T, 1, lower.tail=F) ## p-value of Macro F1* (Score)
  
  
  
  ## ################# ##
  ## Formatting output ##
  ## ################# ##
  
  ftest <- data.frame(rbind(biF=c(biF1, biF2, biF.d, biF.v, biF.T, biF.p),
                            biFs=c(biF1, biF2, biF.d, biFs.v, biFs.T, biFs.p),
                            miF=c(miF1, miF2, miF.d, miF.v, miF.T, miF.p),
                            miFs=c(miF1, miF2, miF.d, miFs.v, miFs.T, miFs.p),
                            maF=c(maF1, maF2, maF.d, maF.v, maF.T, maF.p),
                            maFs=c(maF1, maF2, maF.d, maFs.v, maFs.T, maFs.p),
                            maF.star=c(maF1b, maF2b, maFb.d, maFb.v, maFb.T, maFb.p),
                            maFs.star=c(maF1b, maF2b, maFb.d, maFbs.v, maFbs.T, maFbs.p)))
  names(ftest) <- c("Test1_Est","Test2_Est","Diff_Est","Variance","TestStatistics","p-value")
  return(ftest)
}

## Example -------------

library(MASS)

set.seed(1234567)
mu<-c(0, 0, 0)
sigma<-rbind(
  c(1,0.6,0.6),
  c(0.6,1,0.3),
  c(0.6,0.3,1)
)

mvr <- mvrnorm(5000, mu, sigma)
act = ifelse(mvr[,1]<=-0.5,"A",
             ifelse(mvr[,1]<=-0.2,"B",
                    ifelse(mvr[,1]<=0,"C",
                           ifelse(mvr[,1]<=0.2,"D",
                                  ifelse(mvr[,1]<=0.5,"E","F")))))
t1 = ifelse(mvr[,2]<=-0.5,"A",
            ifelse(mvr[,2]<=-0.2,"B",
                   ifelse(mvr[,2]<=0,"C",
                          ifelse(mvr[,2]<=0.2,"D",
                                 ifelse(mvr[,2]<=0.5,"E","F")))))
t2 = ifelse(mvr[,3]<=-0.4,"A",
            ifelse(mvr[,3]<=-0.1,"B",
                   ifelse(mvr[,3]<=0.1,"C",
                          ifelse(mvr[,3]<=0.3,"D",
                                 ifelse(mvr[,3]<=0.6,"E","F")))))

exds = data.frame(act,t1,t2)

f1scores_test(test1=exds$t1,test2=exds$t2,actual=exds$act,binary1=c("A","B"))

## End ##
