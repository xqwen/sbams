library("lars")


d = read.table("sim.r.dat")
p = dim(d)[2]-3

X = as.matrix(d[,4:(p+3)])
mx = apply(X, 2, mean)
X=scale(X, center=mx, scale=TRUE)
X = diag(c(1,1,1))%x%X



Y = as.matrix(d[,1:3])
my = apply(Y, 2, mean)
Y=scale(Y, center=my, scale=TRUE)
Y = matrix(ncol=1,Y)



rst = lars(X,Y,type="lasso")

v = rst$lambda

sink("lasso.out");

#out = cbind(v%%p, round(2^(v%/%p)))
for (i in 1:length(v)){
   if(length(rst$act[[i]])==1 && rst$act[[i]]>0){
     out = rst$act[[i]]
     id = out%%p
     cfg = out%/%p
     if(id==0){
       id = p
       cfg = cfg-1
     }
     cat(sprintf("  %7.5e   rs%d  %d\n",v[i],id,2^(cfg)))
   }   
 

}


sink()