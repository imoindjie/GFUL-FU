## This script is used to aggregate the results (in folders d_g=...) obtained by 'Framework.R'. 
## It uses also functions defined in the scripts in 'Functions'. 
## For any information or issues while running this code, feel free to contact me at: issam-ali.moindjie@inria.fr (or issam.moindjie@gmail.com)

#Useful Package
library(spdep)
library(MFPCA)
library(caret)
library(grplasso)
library(Matrix)
library(MLmetrics)
library(utils) 
library(MASS)

wd=rstudioapi::getActiveDocumentContext()$path
setwd(paste0(head(strsplit(wd, '/')[[1]], -1), collapse = '/'))


dir_scripts='../Functions'
for( i in grep(list.files(dir_scripts), pattern='.R', value = T))source(paste0(dir_scripts, '/', i) ) ## Chargement des fonctions 

## Some functions : 

Sens_M<-function(M_t, M_p){ 
  mean(M_p[which(M_t==0, arr.ind = T)]==0) ## Get sensitivity from Dist matrices 
}

Spec_M<-function(M_t, M_p){ 
  mean(M_p[which(M_t!=0, arr.ind = T)]!=0) ## Get specificity from Dist matrices
}

upper<-function(M)as.numeric(M[upper.tri(M)]) ## Get the upper matrix of M (without the diagonal)

d_g=3
{ 
  p<-d_g*4
  dir_f<-paste0('d_g=', d_g)
  Dist<-readRDS(paste0(dir_f, '/Dist.RDS'))
  
  Mse<-t(as.data.frame(readRDS(paste0(dir_f, '/mse.RDS'))))
  colnames(Mse)<-c('LA','GL','FL','G0', 'GF')
  
  M_Mse<-sapply(as.data.frame(Mse), mean)
  S_Mse<-sapply(as.data.frame(Mse), sd)
  
  betas<-readRDS(paste0(dir_f, '/betas.RDS'))
  
  Dist_G0<-Dist[[1]]$Tr
  Dist_G0[(2*d_g+1): (3*d_g), (2*d_g+1): (3*d_g)]<-0
  
  Sens<-Spec<-NULL
  for(k in 1:length(Dist)){ 
    st<-c(
    Sens_M(upper(Dist[[k]]$Tr), upper(Dist[[k]]$LA) ), 
    Sens_M(upper(Dist[[k]]$Tr), upper(Dist[[k]]$GL) ), 
    Sens_M(upper(Dist[[k]]$Tr), upper(Dist[[k]]$FL) ), 
    Sens_M(upper(Dist[[k]]$Tr), upper(Dist_G0) ), 
    Sens_M(upper(Dist[[k]]$Tr), upper(Dist[[k]]$GF) )) 
    
    sp<-c(
      Spec_M(upper(Dist[[k]]$Tr), upper(Dist[[k]]$LA) ), 
      Spec_M(upper(Dist[[k]]$Tr), upper(Dist[[k]]$GL) ), 
      Spec_M(upper(Dist[[k]]$Tr), upper(Dist[[k]]$FL) ), 
      Spec_M(upper(Dist[[k]]$Tr), upper(Dist_G0) ), 
      Spec_M(upper(Dist[[k]]$Tr), upper(Dist[[k]]$GF) )
      )
    Sens<-rbind(Sens, st)
    Spec<-rbind(Spec, sp)
    
    }
  Sens<-as.data.frame(Sens)
  Spec<-as.data.frame(Spec)
  Sens[is.na(Sens)]<-0
  Spec[is.na(Spec)]<-0
  
  M_Sens<-sapply(Sens, mean)
  M_Spec<-sapply(Spec, mean)
  S_Sens<-sapply(Sens, sd)
  S_Spec<-sapply(Spec, sd)
  
  
  colnames(Spec)<-colnames(Sens)<-c('GL1','GL2','FU','HG', 'GFUL')
  
  X_table<-as.data.frame(cbind(paste0(round(M_Mse, 2), '(',round(S_Mse,2) ,')'), 
        paste0(round(M_Sens, 2), '(',round(S_Sens,2) ,')'), 
        paste0(round(M_Spec, 2), '(',round(S_Spec,2) ,')')
        ))
  rownames(X_table)<-colnames(Spec)
  colnames(X_table)<-c('MSE', 'Sens', 'Spec' )
  xtable::xtable(X_table[c('GL1', 'GL2', 'FU', 'GFUL', 'HG'), ], rownames=NULL)
 
}  
  
## Estimated functions by groups of conditions

k<-1 ## number of simulations of which we plot 

betas<-readRDS(paste0(dir_f, '/betas.RDS'))

groups=rep(seq(p/d_g), each= d_g)


dir.create('Fig')

for(g in unique(groups)){ 
  png(paste0('Fig/True', g, '.png'), width = 800, height = 600)
  plot(to_row(get_dims(betas[[k]]$Tr, which(groups==g)) ), cex.axis=2, ylab="", xlab='', ylim=c(-4, 4), lwd=2, col='black')
  dev.off()
  
  
  png(paste0('Fig/GFUL', g, '.png'), width = 800, height = 600)
  plot(to_row(get_dims(betas[[k]]$GF, which(groups==g)) ), cex.axis=2, ylab="", xlab='', ylim=c(-4, 4), lwd=2, col='black')
  dev.off()
  
  png(paste0('Fig/FU', g, '.png'), width = 800, height = 600)
  plot(to_row(get_dims(betas[[k]]$FL, which(groups==g)) ), cex.axis=2, ylab="", xlab='', ylim=c(-4, 4), lwd=2, col='black')
  dev.off()
  
  png(paste0('Fig/GL1', g, '.png'), width = 800, height = 600)
  plot(to_row(get_dims(betas[[k]]$LA, which(groups==g)) ), cex.axis=2, ylab="", xlab='', ylim=c(-4, 4), lwd=2, col='black')
  dev.off()
  
  png(paste0('Fig/GL2', g, '.png'), width = 800, height = 600)
  plot(to_row(get_dims(betas[[k]]$GL, which(groups==g)) ), cex.axis=2, ylab="", xlab='', ylim=c(-4, 4), lwd=2, col='black')
  dev.off()
  
  png(paste0('Fig/HG', g, '.png'), width = 800, height = 600)
  plot(to_row(get_dims(betas[[k]]$G0, which(groups==g)) ), cex.axis=2, ylab="", xlab='', ylim=c(-16, 16), lwd=2, col='black')
  dev.off()
  
  
}
