# correlated residual errors, effect sizes are correlated across subgroups

n=100
m = 105
k = 15
p = m*k
q=3

get_cfg<-function(cfg){
   rst = c();
   for (i in 1:q){
       rst = c(rst,cfg %% 2)
       cfg = cfg %/% 2
   }
   return(rst)
}




library(MASS)

#simulation
true.B=matrix(nrow=p, ncol=q,0)
for( i in 1:p){
  bbar = rnorm(1,sd = 1)
  true.B[i,] = rnorm(3,mean=bbar, sd=0.1*abs(bbar))
}


iv = rbinom(m,1,0.03)
indi = matrix( 0, nrow=p, ncol=q)

for (i in which(iv==1)){
  index = (i-1)*k + sample(1:k,1)
  if(runif(1)<0.5){
     indi[index,] = c(1,1,1);
  }else{
     indi[index,] = get_cfg(sample(1:6,1))
  }	      
}
true.B = true.B*indi
true.B = true.B%*%diag(c(1,1.5,2))



X.tr=as.matrix(nrow=n, ncol=p,read.table("sim.geno.region"))
covm = matrix(nrow = q, c(1,0.2,0.8,0.2,1,0.6,0.8,0.6,1))
covm = diag(c(1,1.5,2))%*%covm%*%diag(c(1,1.5,2))

# generate correlated errors
E.tr = mvrnorm(n, rep(0,q), covm)



Y.tr=X.tr%*%true.B+E.tr
mx = apply(X.tr, 2, mean)
my = apply(Y.tr, 2, mean)
X.tr=scale(X.tr, center=mx, scale=FALSE)
Y.tr=scale(Y.tr, center=my, scale=FALSE)

rsv = paste("covariate rs",seq(1,p),sep="")
dX = rbind(rsv,X.tr)
dY = rbind(c("response F","response L","response T"),Y.tr)
write.table(file="sim.dat",t(dY), quote = FALSE, append=FALSE, sep=" ",row.names=FALSE, col.names=FALSE)
write.table(file="sim.dat",t(dX), quote = FALSE, append=TRUE, sep=" ",row.names=FALSE, col.names=FALSE)

write.table(file="sim.r.dat",cbind(Y.tr,X.tr),quote = FALSE, append=FALSE, sep=" ",row.names=FALSE, col.names=FALSE)



V = indi[,1]
for( i in 2:q){
   V = V+indi[,i]*(2**(i-1));
}    

sink(file="sim.truth")
rst = cbind(rsv,V,true.B)
rst[V!=0,]
sink()



