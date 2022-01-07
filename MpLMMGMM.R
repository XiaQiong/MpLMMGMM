
Omics_GMMLasso<-function(y,Dem=NULL,X,Gen,index,K=NULL,returnK=F, maxiter=100, tol=1e-4, standardize=FALSE){
  
  test_index <- index  ##test index
  train_index <- c(1:length(y))[-index]  ##train index
  
  if(standardize){
    Pheno <- y[train_index]  ##training Phenotypes
    mean_adjust=mean(Pheno)
    sd_adjust=sd(Pheno)
    Pheno=scale(Pheno)
  } else {
    Pheno <- y[train_index]  ##training Phenotypes
  }
  
  
  
  Xtrain <- X[train_index,]
  
  if(is.null(K)){
    K=lapply(Gen,Fun_Klinear) ##get kernel function
  }
  list_K=lapply(K,function(kk){
    kk[train_index,train_index]
  })
  
  list_K=lapply(list_K, scaleK)
  
  
  numK=length(list_K)
  list_K[[numK+1]]=diag(length(train_index))
  
  if(is.null(Dem)){
    n=length(Pheno);   
    V=diag(n)
    E=eigen(V)$vectors
    A=E                       ##get A matrix
  }else{
    n=length(Pheno);   
    V=diag(n)-Dem %*% MASS::ginv(t(Dem) %*% Dem) %*% t(Dem)
    E=eigen(V)$vectors
    nc=ncol(Dem)
    A=E[,1:(n-nc)]                       ##get A matrix
  }
  
  Tlist=lapply(list_K,function(kk) c(t(A)%*%kk%*%A))
  T=matrix(unlist(Tlist),nrow=ncol(A)*ncol(A),ncol=length(list_K))
  
  run=TRUE; iter=0;    
  betaold=0
  sigmaold=0
  
  fitX_cv=cv.glmnet(Xtrain, Pheno,nfold = 10,alpha = 1,intercept=FALSE)
  lambda0=chooselambda(fitX_cv,1)
  beta <- as.numeric(coef(fitX_cv, s = lambda0))[-1]
  
  while(run & iter<= maxiter){
    iter=iter+1;
    #print(iter)
    
    gPheno=Pheno-Xtrain %*% beta
    
    lasso_Pheno = c(t(A)%*%gPheno%*%t(gPheno)%*%A)   ##training Phenotypes used in GMMLasso
    
    p.fac = c(rep(1,numK),0)
    fit_cv <- cv.glmnet(T, lasso_Pheno, alpha=1,lower.limits =0,penalty.factor = p.fac,intercept=FALSE)     ## fit the model
    lambda1=chooselambda(fit_cv,1.2)     ##search the optimal lambda
    fit_best <- glmnet(T, lasso_Pheno, alpha = 1, lambda = lambda1, lower.limits =0, penalty.factor = p.fac,intercept=FALSE)       ## fit the model
    sigma= as.matrix(fit_best$beta)  #get the parameter estimations 
    
    TSigma=fit_best$beta[length(list_K)]*diag(length(Pheno))               
    
    for (j in 1:numK) {
      TSigma=TSigma + fit_best$beta[j] * list_K[[j]]
    }
    Vinv=MASS::ginv(TSigma);
    l=round(KFAS::ldl(Vinv, tol = tol),-log10(tol))
    d <- diag(sqrt(diag(l)))
    diag(l) <- 1
    tmp=t(l %*% d);
    Vinvsqrt=round(tmp,-log10(tol))
    Ynew=Vinvsqrt %*% Pheno;
    Xnew=Vinvsqrt %*% Xtrain;
    lasso_cv <- glmnet::cv.glmnet(x = Xnew, y = Ynew,nfold = 10,alpha = 1,intercept=FALSE);
    lambda2=chooselambda(lasso_cv,1)     ##search the optimal lambda
    beta <- as.numeric(coef(lasso_cv, s = lambda2))[-1];
    if(mean((sigma-sigmaold)^2)<tol & mean((beta-betaold)^2)<tol) run = FALSE

    sigmaold=sigma
    betaold=beta

  }
  
  #####################################################################################################################
  #prediction
  Sigma=fit_best$beta[length(list_K)]*diag(length(y))               
  
  for (j in 1:numK) {
    Sigma=Sigma + fit_best$beta[j] * K[[j]]
  }
  
  
  Sigma11=Sigma[test_index,test_index]
  Sigma12=Sigma[test_index,train_index]
  Sigma21=Sigma[train_index,test_index]
  Sigma22=Sigma[train_index,train_index]
  
  
  
  Pheno_pred=X[test_index,] %*% beta + Sigma12 %*% MASS::ginv(Sigma22) %*%(Pheno-Xtrain %*% beta )
  
  if(standardize){
    Pheno_pred=Pheno_pred*sd_adjust+mean_adjust
  } else {
    Pheno_pred=Pheno_pred
  }
  
  
  
  
  combineout=cbind(pred=Pheno_pred,true=y[test_index])
  
  if(returnK) {
    out=list(out=combineout,beta=beta, sigma=sigma, K=K)
  }else {
    out=list(out=combineout,beta=beta,sigma=sigma)
  } 
  
  out
}









