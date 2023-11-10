## Intermediate functions

get_u_B=function(mFData, K){ 
  
  p=length(mFData)
  
  uniExpansions=NULL
  for( j in 1:p){ 
    uniExpansions=c(uniExpansions, list(list(type='splines1D', k = K)) )  ## splines cubic 
  }
  
  uniBasis <- mapply(function(expansion, data) {
    do.call(univDecomp, c(list(funDataObject = data), expansion))
  }, expansion = uniExpansions, data = mFData, SIMPLIFY = FALSE)
  
  uniBasis
  
}

get_dims=function(f_data, ind_i){ 
  if(length(ind_i)>1){
    M_res=list()
    for( i in ind_i){ 
      M_res=append(M_res, f_data[[i]])
    }
    
    multiFunData(M_res)
  }else{ 
    
    #M_res=multiFunData(f_data[[ind_i]])
    M_res=f_data[[ind_i]]
    M_res
  }
}


to_row<-function(fd){ 
  p_dims<-length(fd)
  if(p_dims==1){ 
    fd
  }else{ 
    Xs<-rlist::list.rbind(X(fd))
    funData(argvals(fd)[[1]], Xs)
    }

}

get_L0=function(L, pen=T){ 
  p=ncol(L)
  I=diag(1, nrow=p)
  
  if(pen){ 
    
    C_0=NULL
    for(j in 1:p){ 
      v_i=which(L[j, ]==1)
      C_0=c(C_0, which(L[v_i,  ]==1 )==j) 
    }
    I[C_0, ]=2*I[C_0, ]
    
    
  }
  C_1=NULL
  
  for(i in 1:p){
    v_i=which(L[i, ]==1)
    C_1=c(C_1, i>v_i & which(L[v_i,  ]==1 )==i) 
  }
  
  I[!C_1, ]%*%L
  
}



get_equal<-function(L_0, norm_i){ 
  r=nrow(L_0)
  eq_df=NULL
  for(i in (1:r)){ 
    if( norm_i[i]==0){ 
      eq_df=rbind(eq_df, c(which(L_0[i, ]>0),  which(L_0[i, ]<0)))
      
    }
    
  }
  
  eq_df
  
}

get_eq_b<-function(xt){ 
  if(is.null(xt) ){ 
    NULL
    
  }else{ 
    
    ind_all=NULL
    for( i in 1:nrow(xt)){ 
      
      ij<-xt[i, ]
      
      ind_ij<-apply(xt[-i, ], function(x)any(x%in%ij), MARGIN=1) 
      
      if(any(ind_ij) ){ 
        ij<-c(ij,xt[-i,  ][which(ind_ij), ]  )
        
      }
      
      ind_all=append(ind_all, list(unique(ij)) )    
    }
    unique(sapply(ind_all, function(x) 
      unique(unlist(ind_all[sapply(ind_all, function(y) 
        any(x %in% y))]))))
  }
  
}




predict.GFL<-function(obj_GL, newx){ 
  
  u_new=get_u_B(newx-obj_GL$pred_obj$mean, obj_GL$pred_obj$K)
  
  ## Toutes les transformations d'un coup ! 
  Alpha_new=rlist::list.cbind(
    lapply(u_new, function(x)x$scores)
  )%*%(
    obj_GL$pred_obj$Inv_M%x%diag(1, obj_GL$pred_obj$K)
  )%*%bdiag(lapply(u_new, function(x)expm::sqrtm(x$B) ))
  
  y_hat<-predict(obj_GL$model, Alpha_new%*%obj_GL$pred_obj$U, type='response')
  if(is.null(obj_GL$pred_obj$mu_Y)){ 
    y_hat
  }else{ 
    y_hat+obj_GL$pred_obj$mu_Y
    
    }
}





get_Z<-function(oGp ){ 
  cols<-oGp[[1]]; 
  Gps<-oGp[[2]]
  Z<-NULL
  cs<-NULL
  p<-length(cols)
  for(j in 1:p ){ 
    u_j<-rep(0,p )
    u_j[j]<-1
    c_j<-char2vec(cols[j])
    cs<-c(cs, c_j[c_j %in% Gps  ] )
    n_rep<-sum(c_j %in% Gps)
    Z<-cbind(Z, matrix(rep(u_j, n_rep), ncol=n_rep)   )
    
  }
  
  res<-list(z=Z,ngroups= cs)
  res
}


Xtilde<-function(X_f, Z){ 
  new_cols<-NULL
  p<-length(X_f)
  r_Z<-rowSums(Z)
  for(j in 1:p) new_cols<-c(new_cols, rep(j, r_Z[j]) )
  
  get_dims(X_f, new_cols)
  
}


predict.GoFL<-function(obj_GL, newx){ 
  
  
  u_new=get_u_B(Xtilde(newx, obj_GL$pred_obj$Z)-obj_GL$pred_obj$mean, obj_GL$pred_obj$K)
  
  ## Toutes les transformations d'un coup ! 
  Alpha_new=rlist::list.cbind(
    lapply(u_new, function(x)x$scores)
  )%*%(
    obj_GL$pred_obj$Inv_M%x%diag(1, obj_GL$pred_obj$K)
  )%*%bdiag(lapply(u_new, function(x)expm::sqrtm(x$B) ))
  
  y_hat<-predict(obj_GL$model, Alpha_new%*%obj_GL$pred_obj$U, type='response')
  
  y_hat
  
  
  
  }
char2vec<-function(xx){
  yy<-NULL
  for(i in 1:nchar(xx))yy<-c(yy, substr(xx, i, i))
  yy
}








