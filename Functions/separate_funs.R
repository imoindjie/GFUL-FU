## Functions for the estimation of FU GFUL and GL


Kern<-function(Alpha_1, Y_s, group, sub_group, scal_lambda, k_folds, K, is_binary, cross_val, n, reg_type, verbose, m_iter, folds){ 
  Xi=NULL
  U=NULL
  for(i in 1:length(group)){ 
    alpha=Alpha_1[, group[[i]]]
    eg<-prcomp(alpha, scale. = F, center = F)
    Xi=append(Xi, list(eg$x))# La derniere dimensio est souvent \simeq 0: Mais ça pose pas de problem 
    U=append(U, list(eg$rotation) )
  }
  ## Get it in the local env
  big_Xi<-rlist::list.cbind(Xi)
  big_U<-bdiag(U)
  
  message('Estimation du modèle ...')
  
  
  max_lambda=lambdamax(big_Xi, y=Y_s,
                       index=sub_group, penscale = sqrt, 
                       model=reg_type, center=F, standardize = F) 
  
  lambdas=c(max_lambda*scal_lambda)
  if(cross_val){ 
    message('Validation croisée ...')
    if(is.null(folds))folds<-createFolds(1:n, k_folds)
    
    acc_folds<-NULL
    for(k in 1:k_folds){ 
      
      if(verbose){
        k_lasfit <- grplasso(big_Xi[-folds[[k]], ], y=Y_s[-folds[[k]] ], index=sub_group, 
                             lambda = lambdas, model= reg_type, penscale=sqrt, control = grpl.control(update.hess = "lambda", trace=0, max.iter = m_iter), center=F, standardize = F)
        
      }else{ 
        k_lasfit <- grplasso(big_Xi[-folds[[k]], ], y=Y_s[-folds[[k]] ], index=sub_group, 
                             lambda = lambdas, model= reg_type, penscale=sqrt, control = grpl.control(update.hess = "lambda", trace=0, max.iter =  m_iter), center=F, standardize = F)
        
      }
      if(is_binary){ 
        if(length(lambdas)==1){ 
          acc_folds<-rbind(acc_folds, acc(predict(k_lasfit, big_Xi[folds[[k]], ], type='response')>.5, Y_s[folds[[k]]]) )
        }else{ 
          acc_folds<-rbind(acc_folds, apply(predict(k_lasfit, big_Xi[folds[[k]], ], type='response')>.5, function(x)acc(x, Y_s[folds[[k]]]), MARGIN = 2) )
        }
      }else{
        if(length(lambdas)==1){ 
          acc_folds<-rbind(acc_folds, RMSE(predict(k_lasfit, big_Xi[folds[[k]], ], type='response')>.5, Y_s[folds[[k]]]) )
          
        }else{
          acc_folds<-rbind(acc_folds, apply(predict(k_lasfit, big_Xi[folds[[k]], ]), function(x)RMSE(x, Y_s[folds[[k]]]), MARGIN = 2) )
          }
      
      }
      
      
    }
    if(is_binary){ 
      opt_lam<-lambdas[which.max(round(apply(acc_folds, mean, MARGIN = 2), 2) )] ## Difference de .8
    }else{ 
      opt_lam<-lambdas[which.min(round(apply(acc_folds, mean, MARGIN = 2)))]
    }
    message('Estimation du modèle final ...')
    
    ## Last training ! 
    
    if(verbose){
      lasfit <- grplasso(big_Xi, y=Y_s, index=sub_group, 
                          lambda = opt_lam, model= reg_type, penscale=sqrt, control = grpl.control(update.hess = "lambda", trace=0, max.iter =  m_iter), center=F, standardize = F)
      
    }else{ 
      lasfit <- suppressWarnings(grplasso(big_Xi, y=Y_s, index=sub_group, 
                                           lambda = opt_lam, model= reg_type, penscale=sqrt, 
                                           control = grpl.control(update.hess = "lambda", trace=0, max.iter =  m_iter), center=F, standardize = F)
      )
    }
    
  }else{ 
    
    
    if(verbose){
      lasfit <- grplasso(big_Xi, y=Y_s, index=sub_group, 
                         lambda = lambdas, model= reg_type, penscale=sqrt, control = grpl.control(update.hess = "lambda", trace=0), center=F, standardize = F)
      
    }else{ 
      lasfit <- suppressWarnings(grplasso(big_Xi, y=Y_s, index=sub_group, 
                                          lambda = lambdas, model= reg_type, penscale=sqrt, 
                                          control = grpl.control(update.hess = "lambda", trace=0), center=F, standardize = F)
      )
    }
    opt_lam<-NULL
    acc_folds<-NULL
  }
  
  beta_0<-big_U%*%lasfit$coefficients
  
  list(beta_0,acc_folds, opt_lam, lasfit, big_Xi, big_U )
  
}

# Group fused lasso for non-fully graph 

GFL<-function(X_t, Y, L,pen_T=T, K=10, reg_type=LinReg(),scal_lambda=(9/10)^(seq(0, 100, length.out=50 ) ), cross_val=F, k_folds=10, m_iter=1e4, verbose=T, folds=NULL){ 
  mu<-meanFunction(X_t)
  n<-nObs(X_t)
  u_B=get_u_B(X_t-mu, K)
  Y_s=Y
  is_binary=length(unique(Y)) == 2
  m_Y<-0
  if( !is_binary ){ 
    m_Y<-mean(Y)
    Y_s=scale(Y, center = T, scale = F)
  }
  
  p=length(X_t)
  Alpha=rlist::list.cbind(lapply(u_B, function(x)x$scores) )
  L_1<-get_L0(L)
  r_1=nrow(L_1)
  
  D=rbind(L_1, t(MASS::Null(t(L_1) ) ) )
  
  ## D est d'intêret seulement la partie jusqu'à r  
  
  group=NULL
  for(i in seq(1, ncol(Alpha), by=K)){ 
    group=append(group, list(seq(i, i+(K-1) ) ) )
  }
  
  sub_group=NULL
  for(i in 1:(r_1)){ 
    sub_group=c(sub_group,rep(i, min(K, n) ))
  }
  
  ## si r<< p, ça peut donner des valeurs abérantes si pen_T n'est pas true 
  if(pen_T){ 
    sub_group=c(sub_group, rep(r_1+1, (p-r_1)*min(c(K,n) ) )) ## On penalise le residu ||(T\beta)||_2 \leq s' 
  }else{ 
    sub_group=c(sub_group, rep(NA, (p-r_1)*min(c(K,n) ) )) ## ==> Les valeurs non pénalisées impliquent des gros valeurs. 
    
  }
  
  R_M<-solve(D) ## Retrouver le beta 
  F_sqrt<-bdiag(lapply(u_B, function(x)expm::sqrtm(x$B) ) )
  F_sqrt_inv<-solve(F_sqrt)
  
  Alpha_1=Alpha%*%((R_M%x%diag(1, nrow=K)))%*%F_sqrt
  
  ker_out<-Kern(Alpha_1,Y_s,  group, sub_group, scal_lambda, k_folds, K, is_binary, cross_val, n, reg_type, verbose, m_iter, folds)
  
  beta_0<-ker_out[[1]]; acc_folds<-ker_out[[2]]
  opt_lam<-ker_out[[3]]; lasfit<-ker_out[[4]]
  big_Xi<-ker_out[[5]]; big_U<-ker_out[[6]]
  
  
  beta_coef<-(R_M%x%diag(1,nrow = K))%*%F_sqrt_inv%*%beta_0
  g_index=NULL
  for(i in 1:p) g_index=c(g_index, rep(i, K))
  
  g_index=NULL
  for(i in 1:p){ 
    g_index=c(g_index, rep(i, K))
  }
  
  b_lasso=NULL
  
  for(i in 1:p){ 
    b_lasso=append(b_lasso, list( univExpansion(
      scores = t(as.matrix(beta_coef[which(g_index==i)[1:K],  ])),type='splines1D', 
      argvals = argvals(X_t)[[1]],functions=NULL, params = list(bs = "ps", m = 2, k = K)
    )
    ))
  }
  b_lasso=multiFunData(b_lasso)
  
  res=list(model=lasfit, pred_obj=list(mean=mu, Inv_M=R_M, U=big_U, K=K, mu_Y=m_Y), 
           coef=beta_coef, b_estim=b_lasso, opt_lambda=opt_lam, v_prec=acc_folds)
  
  class(res)<-'GFL'
  res
}

## Grouped fused lasso 

GdFL<-function(X_t, Y,Group,K=10, reg_type=LinReg(),al=0, scal_lambda=(9/10)^(seq(0, 100, length.out=50 ) ), cross_val=F, k_folds=10, m_iter=1e4, verbose=F, folds=NULL ){ 
  
  if((!is.null(Group)) & any(table(Group)<2) ) stop('Tous les groupes doivent avoir au moins deux elements')
  message('Projection base ...')
  mu<-meanFunction(X_t)
  n<-nObs(X_t)
  u_B=get_u_B(X_t-mu, K)
  Y_s=Y
  is_binary=length(unique(Y)) == 2
  m_Y<-0
  if( !is_binary ){ 
    m_Y<-mean(Y)
    Y_s=scale(Y, center = T, scale = F)
  }
  
  p=length(X_t)
  Alpha=rlist::list.cbind(lapply(u_B, function(x)x$scores) )
  S=diag(1, length(Group))[order(Group), ] ## Matrice pour ordonner XS=(X[1], X[2], ...., X[g])
  
  p_i=NULL
  for( i in unique(Group[order(Group)]) ){ 
    p_i=c(p_i, sum(i==Group) )
  }
  
  
  ## nombre de groupes 
  n_g=length(p_i)
    ## M_G and R
  M_G=NULL
  R=NULL
  for(i in 1:n_g){ 
    R_i=matrix(
      qr.R(
        qr(diag(1, nrow=p_i[i])-1/p_i[i]*matrix(1, nrow=p_i[i], ncol=p_i[i]))
      )[1:(p_i[i]-1), ]
      , nrow=(p_i[i]-1))
    
    R=append(R, list(R_i))
    M_G=append(M_G, list(matrix(1/p_i[i], ncol=p_i[i]) )   )  
  }
  
  R=bdiag(R)
  M_G=bdiag(M_G)
  ## introduire alpha 
  if(!al%in%c(0, 1) ){ 
    G=rbind((1-al)*R, al*M_G)
  }else{ 
    
    G=rbind(R, M_G)
    
    }
  Inv_G=solve(G)
  ## 
  R_M=(t(S) %*%Inv_G) ## --> Avoir les beta ## 
  
  #Alpha_1=Alpha_s%*%((S%*%Inv_G)%x%diag(1, K))
  
  F_sqrt<-bdiag(lapply(u_B, function(x)expm::sqrtm(x$B ) ) )
  F_sqrt_inv<-solve(F_sqrt)
  
  Alpha_1=Alpha%*%((R_M%x%diag(1, nrow=K)))%*%F_sqrt
  
  p_i_2=c(p_i-1, rep(1, n_g))
  
  p_K=cumsum(p_i_2*K)
  group=NULL
  for(i in 1:(2*n_g) ){ 
    if(i ==1){ 
      group=append(group, 
                   list(c(1:p_K[i]) ))
    }else{ 
      group=append(group, 
                   list((p_K[i-1]+1):p_K[i] )
      )  
    }
  }
  ## Expetions 
  t_1<-t_0<-1
  if(al==1)t_1<-NA
  if(al==0)t_0<-NA
  
  sub_group=NULL
  
  for(i in 1:n_g ){ 
      sub_group=c(sub_group, rep(i*t_1 , min(p_i_2[i]*K, n ) ))
    }

  for(i in (n_g+1):(2*n_g) ){ 
    sub_group=c(sub_group, rep(i*t_0 , min(p_i_2[i]*K, n ) ))
  }
  
  ker_out<-Kern(Alpha_1,Y_s,  group, sub_group, scal_lambda, k_folds, K, is_binary, cross_val, n, reg_type, verbose, m_iter, folds)
  
  beta_0<-ker_out[[1]]; acc_folds<-ker_out[[2]]
  opt_lam<-ker_out[[3]]; lasfit<-ker_out[[4]]
  big_Xi<-ker_out[[5]]; big_U<-ker_out[[6]]
  
  
  #beta_coef<-(R_M%x%diag(1,nrow = K))%*%beta_0
  beta_coef<-(R_M%x%diag(1,nrow = K))%*%F_sqrt_inv%*%beta_0
  
  g_index=NULL
  for(i in 1:p) g_index=c(g_index, rep(i, K))
  
  g_index=NULL
  for(i in 1:p){ 
    g_index=c(g_index, rep(i, K))
  }
  
  b_lasso=NULL
  
  for(i in 1:p){ 
    b_lasso=append(b_lasso, list( univExpansion(
      scores = t(as.matrix(beta_coef[which(g_index==i)[1:K],  ])),type='splines1D', 
      argvals = argvals(X_t)[[1]],functions=NULL, params = list(bs = "ps", m = 2, k = K)
    )
    ))
  }
  b_lasso=multiFunData(b_lasso)
  
  res=list(model=lasfit, pred_obj=list(mean=mu, Inv_M=R_M, U=big_U, K=K, mu_Y=m_Y), 
           coef=beta_coef, b_estim=b_lasso, is_group=!is.null(Group), opt_lambda=opt_lam, v_prec=acc_folds, al=al)
  
  class(res)<-'GFL'
  res
  
}







## Classical Group lasso 

GLasso<-function(X_t, Y,Group,K=10, reg_type=LinReg(),scal_lambda=(9/10)^(seq(0, 100, length.out=50 ) ), cross_val=F, k_folds=10, m_iter=100, verbose=F, folds=NULL ){ 
  
  message('Projection base ...')
  mu<-meanFunction(X_t)
  n<-nObs(X_t)
  u_B=get_u_B(X_t-mu, K)
  Y_s=Y
  is_binary=length(unique(Y)) == 2
  m_Y<-0
  if( !is_binary ){ 
    m_Y<-mean(Y)
    Y_s=scale(Y, center = T, scale = F)
  }
  p=length(X_t)
  Alpha=rlist::list.cbind(lapply(u_B, function(x)x$scores) )
  
  S=diag(1, length(Group))[order(Group), ] ## Matrice pour ordonner XS=(X[1], X[2], ...., X[g])
  
  p_i=NULL
  for( i in unique(Group[order(Group)]) ){ 
    p_i=c(p_i, sum(i==Group) )
  }
  
  ## nombre de groupes 
  n_g=length(p_i)
  
  R_M=t(S)## --> Avoir les beta ## 
  
  #Alpha_1=Alpha_s%*%((S%*%Inv_G)%x%diag(1, K))
  
  F_sqrt<-bdiag(lapply(u_B, function(x)expm::sqrtm(x$B ) ) )
  F_sqrt_inv<-solve(F_sqrt)
  
  Alpha_1=Alpha%*%((R_M%x%diag(1, nrow=K)))%*%F_sqrt
  
  #Alpha_1=Alpha%*%(R_M%x%diag(1, K))%*%bdiag(lapply(u_B, function(x)x$B ))
  p_i_2=p_i
  
  p_K=cumsum(p_i_2*K)
  
  group=NULL
  for(i in 1:(n_g) ){ 
    if(i ==1){ 
      group=append(group, 
                   list(c(1:p_K[i]) ))
    }else{ 
      group=append(group, 
                   list((p_K[i-1]+1):p_K[i] )
      )  
    }
  }
  sub_group=NULL
  for(i in 1:(n_g) ){ 
    sub_group=c(sub_group, rep(i , min(p_i_2[i]*K, n ) ))
  }
 
  
  ker_out<-Kern(Alpha_1,Y_s,  group, sub_group, scal_lambda, k_folds, K, is_binary, cross_val, n, reg_type, verbose, m_iter, folds)
  
  beta_0<-ker_out[[1]]; acc_folds<-ker_out[[2]]
  opt_lam<-ker_out[[3]]; lasfit<-ker_out[[4]]
  big_Xi<-ker_out[[5]]; big_U<-ker_out[[6]]
  
  
  beta_coef<-(R_M%x%diag(1,nrow = K))%*%F_sqrt_inv%*%beta_0
  g_index=NULL
  for(i in 1:p) g_index=c(g_index, rep(i, K))
  
  g_index=NULL
  for(i in 1:p){ 
    g_index=c(g_index, rep(i, K))
  }
  
  b_lasso=NULL
  
  for(i in 1:p){ 
    b_lasso=append(b_lasso, list( univExpansion(
      scores = t(as.matrix(beta_coef[which(g_index==i)[1:K],  ])),type='splines1D', 
      argvals = argvals(X_t)[[1]],functions=NULL, params = list(bs = "ps", m = 2, k = K)
    )
    ))
  }
  b_lasso=multiFunData(b_lasso)
  
  res=list(model=lasfit, pred_obj=list(mean=mu, Inv_M=R_M, U=big_U, K=K, mu_Y= m_Y), 
           coef=beta_coef, b_estim=b_lasso, is_group=!is.null(Group), opt_lambda=opt_lam, v_prec=acc_folds)
  
  class(res)<-'GFL'
  res
  
  
}




