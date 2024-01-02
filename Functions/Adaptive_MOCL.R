library(ClusterR) # For K-means++

###############################################################################
# reassign cluster function -> to reorder clusters
## we give clusters the order corresponding to y values.
cluster_reassign<-function(y,new_cluster,k){
  original_cluster=y
  n=length(original_cluster)
  c_mat=cbind(original_cluster,new_cluster)
  
  # reassign by mean of original_cluster up to new_cluster
  cc_mat=matrix(0,k,2)
  cc_mat[,1]=c(1:k)
  
  for(i in 1:k){
    cc_mat[i,2]=mean(c_mat[(c_mat[,2]==i),1])
  }
  
  colnames(cc_mat)=c("new_cluster","mean_original_data")
  cc_mat=data.frame(cc_mat)
  cc_mat=cc_mat[order(cc_mat$mean_original_data),]
  
  new_cc=c(1:k)
  cc_mat=cbind(cc_mat,new_cc)
  cc_mat=cc_mat[,c(1,3)]     # remove mean
  
  # merge with new_cluster(non_ordinal_cluster)
  index=c(1:n)
  ordinal_clu_mat=cbind(new_cluster,index)
  ordinal_clu_mat=merge(ordinal_clu_mat,cc_mat,by="new_cluster")
  ordinal_clu_mat=ordinal_clu_mat[order(ordinal_clu_mat$index),]
  
  # new_ordinal_cluster
  new_ordinal_cluster=ordinal_clu_mat$new_cc
  return(new_ordinal_cluster)
}

###############################################################
# minmax scaler for matrix
nor_minmax = function(x){
  xx=matrix(0,nrow(x),ncol(x))
  colnames(xx)=colnames(x)
  for(i in 1:ncol(x)){
    xx[,i]=(x[,i] - min(x[,i])) / (max(x[,i]) - min(x[,i]))
  }
  return(xx)
}


##################################################################################################################
# Monotone Clustering function
## MSLasso_bic.R function is necessary.
Ad_MOCL<-function(Xf,k,max.iter=500,tr_ratio=0.8,num.knotsf=6,max.lambda=2,r_a=0.01,seed=100,delta=0.01){
  
  # Xf : data
  # k : Number of clusters
  # iter : max # of iterations
  # tr_ratio : Proportion of subsample for MSLasso
  # num.knotsf : Number of knots
  # max.lambda : Max lambda value for lambda sequence used for MSLasso. log(max.lambda) is max value.
  # r_a : ridge penalty parameter used to calculate degrees of freedom of MSLasso
  # delta : for stopping rule.
  
  Xf=nor_minmax(Xf);n=nrow(Xf);p=ncol(Xf);data_x=as.matrix(Xf)
  iter=max.iter
  kk=k
  m_lambda=max.lambda
  
  #########################################################################
  #########################################################################
  # PC1 -> K-means++ -> initial cluster score
  set.seed(seed*13)
  pca_dt <- prcomp(data_x,center = T,scale. = T) # scaling -> PCA
  PC1=pca_dt$x[,1] # PC1
  km_PC1=KMeans_rcpp(as.matrix(PC1), clusters = kk, num_init = 5, max_iters = 200, 
                     initializer = 'kmeans++')### apply K-means++ to PC1
  
  yhat_init=PC1 # initial estimated cluster score (continuous)
  init_clu=km_PC1$cluster # initial cluster assignments
  y=yhat_init*0  # for discrete initial cluster score 
  for(ks in 1:kk){
    y[init_clu==ks]=mean(yhat_init[init_clu==ks])
  }
  init_clu_rank=cluster_reassign(y=y,new_cluster=init_clu,k=kk) # rank of initial ordinal clusters
  
  
  #########################################################################
  #########################################################################
  # iteration to find monotone clusters
  set.seed(seed*20)
  ss_uni=round(runif(iter,1,iter*10)) # for randomization
  
  cluster_matrix=matrix(0,nrow(data_x),iter) # matrix to save cluster rank
  cluster_matrix_order=matrix(0,nrow(data_x),iter) # matrix to save cluster score
  stop=0 # to check the number of iteration
  tr_mat=matrix(0,n*tr_ratio,iter) # matrix to save subsample for MSLasso training
  # Monotone spline
  Z<-monotone.splines(data_x, num.knotsf)  
  groups <- as.vector(t(matrix(rep((1):(ncol(data_x)),(num.knotsf+2)),ncol(data_x),(num.knotsf+2))))
  objective=matrix(0,iter,1)
  
  # Initial information
  clu_rank=init_clu_rank
  y=yhat_init 
  
  for(i in 1:iter){
    
    set.seed(ss_uni[i])
    # subsampling
    sub_sam=sample(1:n,n*tr_ratio)
    tr_mat[,i]=sub_sam
    
    Z_tr=Z[sub_sam,]
    y_tr=y[sub_sam]
    
    #################################################################
    #################################################################
    # Fitting MS-lasso (MAM step)
    ms_lasso=ms_lasso_coeff(Xf=data_x,Zf=Z_tr,Yf=y_tr,num.knotsf=num.knotsf,
                            max.lambda=max.lambda,r_a=r_a)
    
    ms_lasso_coef=ms_lasso$coef
    y_mean=ms_lasso$b0
    
    # prediction
    y_hat=Z%*%ms_lasso_coef+y_mean # estimated clusters score for all observation
    
    # Matrix of estimated monotone function of each variable
    function_mat=matrix(0,nrow(data_x),ncol(data_x))
    for(j in 1:ncol(data_x)){
      coef_j=ms_lasso_coef[groups==j,]
      x_j=Z[,groups==j]
      function_mat[,j]=x_j%*%coef_j
    }
    
    #################################################################
    #################################################################
    # for initial centers of UCF step
    km_mean=rep(0,kk)
    for(r in 1:kk){
      y_hat_r=y_hat[clu_rank==r,]
      km_mean[r]=mean(y_hat_r)
    }
    # K-means
    km_new=kmeans(y_hat,centers=c(km_mean))
    o_c=km_new$cluster # cluster assignments
    clu_rank=cluster_reassign(y=y,new_cluster=o_c,k=kk) # cluster rank
    # cluster score vector
    oc_order=o_c*0 
    for(ks in 1:kk){
      oc_order[clu_rank==ks]=mean(y_hat[clu_rank==ks])
    }
    
    cluster_matrix[,i]=clu_rank # save cluster
    cluster_matrix_order[,i]=oc_order # save cluster score vector
    
    #################################################################
    #################################################################
    # Stopping rule
    if(i>1){
      if(sum(cluster_matrix[,i-1]!=cluster_matrix[,i])<(max(1,floor(delta*n))+1)){
      break
      }
    }
    #update cluster score
    y=oc_order
    # save number of iteration
    stop=stop+1
  }
  
  ########################### # retrun
  cluster_matrix=cluster_matrix[,c(1:(stop+1))]
  cluster_matrix_score=cluster_matrix_order[,c(1:(stop+1))]
  cluster=clu_rank
  
  if(stop>=max.iter){
    return("no convergence")
  }else{
    return(list(cluster=cluster, cluster_matrix=cluster_matrix,cluster_matrix_score=cluster_matrix_score,m_lambda=m_lambda,
                iteration=(stop+1),fx=function_mat,coef=ms_lasso_coef,y_tr_mean=y_mean))
  }
}
