## This script gives the implementation of the simulation framework. 
## All methods estimations are done using the scripts in 'Functions' folders. 
## For any information or issues while running this code, feel free to contact me at: issam-ali.moindjie@inria.fr (or issam.moindjie@gmail.com)

library(spdep)
library(MFPCA)
library(caret)
library(grplasso)
library(Matrix)
library(MLmetrics)
library(utils)
library(MASS)


## 
wd=rstudioapi::getActiveDocumentContext()$path
setwd(paste0(head(strsplit(wd, '/')[[1]], -1), collapse = '/'))

dir_scripts='../Functions'
for( i in grep(list.files(dir_scripts), pattern='.R', value = T))source(paste0(dir_scripts, '/', i) ) ## Chargement des fonctions 

## Some useful functions 

tr<-function(temp, orig){ 
  res<-1-1/5*(10*abs(temp-orig))^2 ## The triangle function 
  res[res<0]<-0
  res
}

dist_m<-function(x){ 
  sqrt(rowSums(x^2)) ## Euclidean distance 
  
}

get_dist<-function(f){
  M_dist<-matrix(0, nrow=length(f), ncol=length(f))
  for(i in 1:length(f)){ 
    for(j in 1:length(f)){ 
      M_dist[i, j]<-round(sqrt(scalarProduct(f[[i]]-f[[j]], f[[i]]-f[[j]])), 10)
    }
    
  }
  M_dist ## Distance in L_2
  
}

get_conditions<-function(x, c_1){ 
  matrix(c(cos(x*2*pi), 
           sin(x*2*pi))+c_1, ncol=2, byrow = F) ## Give points around the circle of rayon= 1 and centers =c_1
  
}

### 

N=200 ## Number of observations 
n_g<-4 ## Number of groups 
d_g =3 ## Number of elements in the group (Scenario 1: 3, Scenario 2: 20)
nsims<-100 ## Number of replication (100 in the article)

t_vector<-seq(0, 1, length.out=100) ## the grid observations in [0,1]

F_bases<-matrix(c(tr(t_vector, .1), tr(t_vector, .2),tr(t_vector, .3), 
                  tr(t_vector, .4), tr(t_vector, .5), tr(t_vector, .6), 
                  tr(t_vector, .7), tr(t_vector, .8), tr(t_vector, .9)), byrow = F, ncol=9 ) ## The Delta functions in the text

## For the conditions (locations)

cs<-c(rep(-3, d_g), rep(0, d_g), rep(3, d_g), rep(6, d_g)) ## Circles centers 

## Outputs in specific folder by d_g

dir_f<-paste0('d_g=',d_g)
dir.create(dir_f, showWarnings = F)

p<-n_g*d_g ## Number of dimensions 

## Saved ouptut

Mse<-NULL ## Mean squared error 
Dist<-NULL ## The matrix of differences : Dist=(||\beta^{(j)}-\beta^{(k)}||)_{k,j}
betas<-NULL ## The estimated beta by methods 

for(nn in 1){ 
  print(nn)
  cent_r<-unique(cs)*matrix(1, nrow=n_g, ncol=2)
  
  ## Generation of conditions 
  
  l<-NULL
  for(j in 1:p){ 
    l<-rbind(l, get_conditions((j%%d_g)/d_g,cs[j] ))  
  }
  
  n_obj=knearneigh(l, 1)
  W_t=knn2nb(n_obj)
  
  plot(W_t, l)
  
  z_s<-lapply(1:n_g,function(x)as.numeric(round(dist_m(l-cent_r[rep(x, p), ]))==1) )
  
  Group<-apply(data.frame(lapply(1:n_g,function(x)ifelse(z_s[[x]]==1, x, '')  )), function(x)paste0(x, collapse=''), MARGIN = 1)
  
  h_1<-funData(t_vector, t(0*F_bases[, 1]) )
  h_2<-funData(t_vector, t(sqrt(2)*rowSums(F_bases[, 1:3]) ))
  
  b_coef<-NULL
  
  for(j in 1:p){
    h_3<-funData(t_vector, t(F_bases%*%c(rep(1, 9)) ))
    
    h_4<-funData(t_vector, t(F_bases%*%c(rep(0, 4), rep(1, 5)) ))
    
    b_coef<-c(b_coef,
              z_s[[1]][j]*h_1+z_s[[2]][j]*h_2+(-1)^j*(1+j%%d_g)/d_g *z_s[[3]][j]*h_3-z_s[[4]][j]*h_2)  
    
  }
  
  
  b_coef<-multiFunData(b_coef)
  
  a_j<-matrix(rnorm(9*N*p, 0, 1 ), ncol=p)
  
  X_all<-NULL
  for(j in 1:p){ 
    X_all<-c(X_all, 
             funData(t_vector, t(F_bases%*% matrix(a_j[, j], ncol=N) ) ) 
    )
  }
  
  X_all<-multiFunData(X_all)
  
  Y_0<-(scalarProduct(b_coef, X_all)) ## : \langle X, \beta \rangle_\mathcal{H}
  
  Y_all<-Y_0+rnorm(n=N, sd= ifelse(d_g==20, 3.6, 1.8) ) 
  
  
  ## Split Train and Test (on a 0.8 and 0.2 rate)
  n_train<-sample(1:N, N*.8)
  X_1<-X_all[n_train]; X_v<-X_all[-n_train]
  Y_1<-Y_all[n_train]; Y_v<-Y_all[-n_train]
  folds<-createFolds(1:nObs(X_1), 10)
  
  ## GFUL estimation
  
  # Estimations for alpha in \{0.1, 0.2, \ldots, 1\}
  GdFs<-NULL
  for(al in seq(0.1, 1, by=.1)){ 
    GdFs<-c(GdFs, list(GdFL(X_1, Y_1, K=20, verbose = T, 
                            scal_lambda =c(1, .96^(seq(1, 148, 1) ), 0), cross_val = T,
                            Group = Group, al=al, m_iter=50, folds=folds)) )
    
  }
  
  ## Which alpha gives the best cross-validation performance ? ( minimal prec) 
  prec<-NULL
  for(i in 1:length(GdFs)){ 
    G_i<-GdFs[[i]]
    prec<-c(prec, mean(G_i$v_prec[, match(G_i$opt_lambda, colnames(G_i$v_prec))]) )
    
  }
  
  GF<-GdFs[[which.min(prec)]] ## The "optimal" GFUL in terms of \alpha and \lambda  
  
  ## HG estimation 
  
  G0<-GdFL(X_1, Y_1, K=20, verbose = T, 
           scal_lambda =1 , cross_val = F,
           Group = Group, al=0, m_iter=50)
  
  ## FU estimation 
  
  L=nb2mat(W_t)-diag(1, nrow=nrow(l)) ## Compute the L matrix (W-I)
  
  FL<-GFL(X_1, Y_1,L, pen_T = T,  K=20, 
          scal_lambda =c(1, .96^(seq(1, 148, 1) ), 0), 
          cross_val = T, m_iter=50, folds=folds) ## FU model in the text 
  
  ## GL 1 estimation, each X^{(j)} forms a group 
  LA<-GLasso(X_1, Y_1, Group = 1:length(X_1), 
             scal_lambda =c(1, .96^(seq(1, 148, 1) ), 0), 
             cross_val =T, m_iter=50, folds=folds) ## GL1 in the text
  
  ## GL 2 estimation, each X^{(j)} forms a group 
  GL<-GLasso(X_1, Y_1, Group = Group, 
             scal_lambda =c(1, .96^(seq(1, 148, 1) ), 0), 
             cross_val =T, m_iter=50, folds=folds  ) ## GL2 in the text 
  
  ## Obtained metrics for saving 
  
  mse<-c( MSE(predict(LA, X_v), Y_v), 
          MSE(predict(GL, X_v), Y_v),
          MSE(predict(FL, X_v), Y_v), 
          MSE(predict(G0, X_v), Y_v), 
          MSE(predict(GF, X_v), Y_v)
  )
  
  print(mse)
  Mse<-append(Mse, list(mse))
  
  Dist<-append(Dist,  
               list(list(GL=get_dist(GL$b_estim), 
                         GF=get_dist(GF$b_estim), 
                         FL=get_dist(FL$b_estim), 
                         LA=get_dist(LA$b_estim), 
                         Tr=get_dist(b_coef))) ) # Note that betas$Tr is the true coefficient function 
  
  betas<-append(betas,  
                list(list(GL=GL$b_estim,  
                          GF=GF$b_estim, 
                          FL=FL$b_estim,
                          LA=LA$b_estim, 
                          G0= G0$b_estim,
                          Tr=b_coef)) )
  ## Saving ...
  saveRDS(file=paste0(dir_f, '/mse.RDS'),Mse  )
  saveRDS(file=paste0(dir_f, '/Dist.RDS'), Dist)
  saveRDS(file=paste0(dir_f, '/betas.RDS'), betas)
}
