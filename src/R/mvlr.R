library(MASS)	

compute_log10_ABF<-function(Y, Xg, Xc, Wg, alpha=0, H=0, nu=0, ES=1){
			       
   n = dim(Y)[1]
   s = dim(Y)[2]
   q = dim(Xc)[2]
   p = dim(Xg)[2]   

   X = cbind(Xc,Xg)
   
   Sigma_hat = Sigma_tilde = diag(rep(0,s))

   if(H == 0){ 
      H = diag(rep(1e-8,s))
   }


   if(alpha > 0){
       Sigma_hat = (t(Y)%*%(diag(rep(1,n)) - X%*%ginv(t(X)%*%X)%*%t(X))%*%Y)/n
   }

   if(alpha<1){
       Sigma_tilde = (t(Y)%*%(diag(rep(1,n)) - Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc))%*%Y)/n
   }
 



   Sigma = (nu/(nu+n))*H + (n/(n+nu))*(alpha*Sigma_hat + (1-alpha)*Sigma_tilde)
   Sigma_inv = solve(Sigma)

   
   # if ES == 1, use the scale-invariant prior formulation
   if(ES == 1){
      S = diag(rep(1,p))%x%diag(sqrt(diag(Sigma)))
      Wg = S%*%Wg%*%S
   }

   Vg_inv = (t(Xg)%*%Xg - t(Xg)%*%Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc)%*%Xg)%x%Sigma_inv
   vec = matrix(ncol=1, as.vector(t(Y-Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc)%*%Y)))
   bVi = t(vec)%*%(Xg%x%Sigma_inv)
   ivw = diag(rep(1,p*s))+Vg_inv%*%Wg
   log10_abf = (.5*bVi%*%Wg%*%solve(ivw)%*%t(bVi)-0.5*determinant(ivw)$modulus[[1]])/log(10)

   return(log10_abf)
   
}


make_Wg<-function(phi, omega, p, r){
 
   phi2 = phi^2
   omg2 = omega^2	
   W = matrix(ncol=r, rep(omg2,r*r)) + diag(rep(phi2,r))
   Wg = diag(rep(1,p))%x%W
   return(Wg)
}




#########################################################


#####################################################

attach(read.table("sim.R.dat",head=T))
Y = cbind(t1,t2,t3)
n = dim(Y)[1]

Xc = matrix(nrow=n,rep(1,n))
Xg = matrix(nrow=n,cbind(rs1, rs2, rs3))
Wg = make_Wg(phi=0.5,omega=0.5 ,p=dim(Xg)[2],r=dim(Y)[2])
compute_log10_ABF(Y,Xg,Xc, Wg)



