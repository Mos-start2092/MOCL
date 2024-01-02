#  The "scoop" package from "http://www.math-evry.cnrs.fr/logiciels/scoop" is necessary for MSLasso.
library(scoop)

# Monotone spline basis generator
## functions for computation of I-Splines of degree 2,cf Ramsay 88
## It's from "mslasso2.r" package

isplineDesign <- function(x,knot.vec){
  I <- numeric()
  deg <- 2
  for(i in 1:(length(knot.vec)-deg)){
    if(x<knot.vec[i]){
      I[i] <- 0
    }
    if (x>knot.vec[i+2]) {
      I[i] <- 1
    }else{
      if(x>=knot.vec[i] && x<=knot.vec[i+1]){
        I[i] <- (x-knot.vec[i])^2 / ((knot.vec[i+1]-knot.vec[i])*(knot.vec[i+2]-knot.vec[i]))
      }
      if(x>=knot.vec[i+1] && x<=knot.vec[i+2]){
        I[i] <- 1-(knot.vec[i+2]-x)^2 / ((knot.vec[i+2]-knot.vec[i])*(knot.vec[i+2]-knot.vec[i+1]))
      }
    }
  }
  return(I)
}

monotone.basis <- function(X,k,spline,xl,xr,intercept){
  n <- length(X)
  dx <- (xr-xl)/(k-1)
  num.knots <- k-2
  
  if(spline=="ispline"){
    deg <- 2
    knot.vec <- seq(xl-(deg-1)*dx,xr+(deg-1)*dx,by=dx) 
    knot.vec[(deg+1):(num.knots+2)] <- quantile(unique(X), seq(0,1,length = (num.knots+2)))[-c(1,(num.knots+2))]
    B <- t(apply(as.matrix(X),1,"isplineDesign",knot.vec=knot.vec))
  } 
  res <- list("B"=B, "knot.vec"=knot.vec)
  return(res)
}

monotone.splines <- function(Xf,num.knots){
  pf <- ncol(Xf)
  Z <- NULL
  for(j in 1:pf){
    x <- Xf[,j]
    spline <- monotone.basis(x, (num.knots+2), "ispline", xl = min(x), xr = max(x), FALSE)
    z1 <- spline$B
    Z=cbind(Z,z1)
  }
  Z <- scale(Z, center = TRUE, scale = FALSE)
  return(Z)
}


########################################################################################
## The ridge regression coefficient vector is needed to calculate degrees of freedom of Ms-lasso.
## Since we choose penalty parameter lambda based on BIC, d.f. must be calculated.
ridge_spline<-function(Zf,Yf,ridge_a){
  r_a=ridge_a
  Yc=Yf-mean(Yf)
  ridge_beta=solve(t(Zf)%*%Zf+r_a*diag(ncol(Zf)))%*%t(Zf)%*%Yc
  
  ridge_L=Zf%*%solve(t(Zf)%*%Zf+r_a*diag(ncol(Zf)))%*%t(Zf) # to calculate df
  ridge_df=2*sum(diag(ridge_L))-sum(diag(t(ridge_L)%*%ridge_L))
  
  ridge_yhat=ridge_L%*%Yc # centered yhat
  SSE_ridge=t(Yc-ridge_yhat)%*%(Yc-ridge_yhat)
  var_ridge=SSE_ridge/(nrow(Zf)-ridge_df)
  
  return(list(var=c(var_ridge),coef=ridge_beta,r_a=r_a))
}

########################################################################################
# MS-lasso function using "scoop" package from "http://www.math-evry.cnrs.fr/logiciels/scoop".
#library(scoop)

## MS-lasso
ms_lasso_coeff<-function(Xf,Zf,Yf,num.knotsf,max.lambda,r_a=0.01){
  
  groups <- as.vector(t(matrix(rep((1):(ncol(Xf)),(num.knotsf+2)),ncol(Xf),(num.knotsf+2)))) # group list for splines
  Yc<-Yf-mean(Yf) # centered cluster score
  b0=mean(Yf)
  lambda_length=100 # number of penalty parameter list
  m_lambda=log(max.lambda) # max lambda for penalty parameter list
  #########################
  # fitting until having large lambda enough to make sign-coherent coefficients 
  coef_0=0
  while(coef_0==0){
    lambda_seq=seq.default(from=m_lambda,to=log(1),length=lambda_length)
    coopfit_model <- coop.lasso(Zf,Yc, groups, family = "gaussian", intercept = FALSE, lambda = lambda_seq)
    msfit <- coopfit_model
    coef_coop<-slot(msfit,"coefficients")
    coef_0=sum(apply(abs(coef_coop),1,sum)==0) # If there exist lambda guarantee monotonocity, coef_0 > 0
    m_lambda=m_lambda*2 # If coef_0 = 0, we refit with lambda list generated with larger max labda value.
  }
  #########################
  ### To get lambda having at least 3 row with no only 0 coefficients 
  ### If max lambda is too large, the interval between lambdas in lambda seq is too large.
  ### It can make all coefficient vectors corresponding to lambda only 0 or sign-incoherent.
  coef_no_0=sum(apply(abs(coef_coop),1,sum)>0)
  while(coef_no_0<3){
    m_lambda=lambda_seq[length(lambda_seq)-2] 
    lambda_seq=seq.default(from=m_lambda,to=log(1),length=lambda_length)
    coopfit_model <- coop.lasso(Zf,Yc, groups, family = "gaussian", intercept = FALSE, lambda = lambda_seq)
    msfit <- coopfit_model
    coef_coop<-slot(msfit,"coefficients")
    coef_no_0=sum(apply(abs(coef_coop),1,sum)>0)
  }
  
  m_lambda_new=lambda_seq[(length(lambda_seq))-coef_no_0]
  lambda_seq_new=lambda_seq
  
  ### The lambda values corresponding to sign coherent coefficients 
  coef_coop_no0=coef_coop[c((nrow(coef_coop)-coef_no_0):nrow(coef_coop)),]
  coef_sign_mat=matrix(0,nrow(coef_coop_no0),ncol(Xf))
  for(i in 1:ncol(Xf)){
    coef_g=coef_coop_no0[,groups==i]
    a=apply(abs(coef_g)+coef_g,1,sum)*apply(abs(-coef_g)-coef_g,1,sum)
    coef_sign_mat[,i]=a # incoherent -> +, sign-coherent ->0
  }
  sign_co_lst= (apply(coef_sign_mat,1,sum)!=0)+(apply(abs(coef_coop_no0),1,sum)==0) 
  sign_co=sum(sign_co_lst==0) # number of lambda with sign-coherent coeff 
  length_new=100
  
  ## If there is no more than 2 lambda values corresponding to ign-coherent coeff, make the length of lambda seq large.
  while(sign_co<2){
    
    lambda_seq_new=seq.default(from=m_lambda_new,to=log(1),length=length_new)
    coopfit_model <- coop.lasso(Zf,Yc, groups, family = "gaussian", intercept = FALSE,
                                lambda = lambda_seq_new)
    msfit <- coopfit_model
    coef_coop<-slot(msfit,"coefficients")
    coef_no_0=sum(apply(abs(coef_coop),1,sum)>0)
    coef_coop_no0=coef_coop[c((nrow(coef_coop)-coef_no_0):nrow(coef_coop)),]
    coef_sign_mat=matrix(0,nrow(coef_coop_no0),ncol(Xf))
    for(i in 1:ncol(Xf)){
      coef_g=coef_coop_no0[,groups==i]
      a=apply(abs(coef_g)+coef_g,1,sum)*apply(abs(-coef_g)-coef_g,1,sum)
      coef_sign_mat[,i]=a # incoherent -> +, sign-coherent ->0
    }
    ##################### not sign incoherent  +  not 0 => if 0, sign coherent      
    sign_co_lst= (apply(coef_sign_mat,1,sum)!=0)+(apply(abs(coef_coop_no0),1,sum)==0) 
    
    sign_co=sum(sign_co_lst==0) # number of lambda with sign-coherent coeff 
    length_new=length_new*2
  }
  
  #########################
  #### lambda vector and coef mat giving sign coherent coeff
  coef_coop_m=coef_coop_no0[(sign_co_lst==0),]  
  lambda_seq_new_m=lambda_seq_new[(sign_co_lst==0)]
  
  #########################
  # yhat(estimated cluster score vector) of each lambda & SSE
  yhat_mat=Zf%*%t(as.matrix(coef_coop_m)) # nrow(zf) x w, w : length of lambda_seq_new_m
  diff_mat=(yhat_mat-t(matrix(rep(Yc,ncol(yhat_mat)),ncol(yhat_mat),nrow(yhat_mat),by=T) ) )^2
  SSE_list=apply(diff_mat,2,sum) # length=w  list
  
  #########################
  # ridge for calculation of d.f. (ridge parameter # r_a=0.01 is default)
  ridge_m=ridge_spline(Zf=Zf,Yf=Yc,ridge_a=r_a)
  ridge_beta=ridge_m$coef # ridge coef
  r_a=ridge_m$r_a 
  var_ridge=ridge_m$var # estimated variance by ridge regression
  
  #########################
  # df up to lambda
  w=nrow(coef_coop_m)
  df_list=c(rep(0,w))
  ridge_mat=t(matrix(rep(ridge_beta,w),ncol(Zf)))
  # sum df up to group by group 
  for(j in 1:ncol(Xf)){
    coef_g=coef_coop_m[,groups==j] # group coef matrix
    ridge_g=ridge_mat[c(1:w),groups==j] # group coef in ridge
    # identify matrix
    coef_g_p_i=(coef_g>0)*1
    coef_g_m_i=(coef_g<0)*1
    # number of sign group
    p_g_p=apply(coef_g_p_i,1,sum)
    p_g_m=apply(coef_g_m_i,1,sum)
    # identity function of sign group 
    i_g_p=(p_g_p>0)*1
    i_g_m=(p_g_m>0)*1
    # ridge sign coef mat
    ridge_g_p=coef_g_p_i*ridge_g
    ridge_g_m=coef_g_m_i*ridge_g
    # L2 norm of each sign coefficients
    ridge_p_norm=sqrt(apply(ridge_g_p^2,1,sum))
    ridge_m_norm=sqrt(apply(ridge_g_m^2,1,sum))
    coop_p_norm=sqrt(apply((coef_g_p_i*coef_g)^2,1,sum))
    coop_m_norm=sqrt(apply((coef_g_m_i*coef_g)^2,1,sum))
    
    # division of norms
    norm_p_d=coop_p_norm/ridge_p_norm; norm_p_d[is.na(norm_p_d)] <- 0 # na -> 0 to calculate
    norm_m_d=coop_m_norm/ridge_m_norm; norm_m_d[is.na(norm_m_d)] <- 0
    
    # df
    df_p=i_g_p*(norm_p_d*(p_g_p-1)/(1+r_a)+1) # r_a : ridge parameter
    df_m=i_g_m*(norm_m_d*(p_g_m-1)/(1+r_a)+1)
    
    df_list=df_list+(df_p+df_m)
  }
  
  ##########################
  ## Select lambda bound to  satisfy the conditions (monotonicity, df > 0) 
  # Monotone
  BIC_list=SSE_list/var_ridge+log(nrow(Zf))*df_list
  BIC_list_m=BIC_list
  SSE_list_m=SSE_list
  
  # Select  lambda only satisfying condition(df>0) 
  condition_index=(apply(abs(coef_coop_m),1,sum)>0)
  coef_coop_condition=coef_coop_m[condition_index,]
  BIC_list_condition=BIC_list_m[condition_index]
  SSE_list_condition=SSE_list_m[condition_index]
  
  #########################
  # choose lambda based on BIC and MSE # default of our method is BIC
  mw_bic=as.integer(which.min(BIC_list_condition))
  mw_MSE=as.integer(which.min(SSE_list_condition))
  if(sum(condition_index)==1){
    coef_m_bic=as.matrix(coef_coop_condition,1)
    coef_m_mse=as.matrix(coef_coop_condition,1)
  }else{
    coef_m_bic=as.matrix(coef_coop_condition[mw_bic,],1)
    coef_m_mse=as.matrix(coef_coop_condition[mw_MSE,],1)
  }
  
  return(list(coef=coef_m_bic,coef_mse=coef_m_mse,coef_mat=coef_coop,bic=BIC_list,df=df_list,b0=b0))
  # coef : coefficients corresponding lambda selected based on BIC
  # coef_mse : coefficients corresponding lambda selected based on mse
  # coef_mat : matrix of all coefficients of all lambda
  # bic : bic list / df : degrees of freedom list / b0 = mean of y (cluster score vector)
}


