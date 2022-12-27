

// Author: Derya Uysal
// Version: 12/27/2022
// 
//Function for the asmptotic variance of  the regression estimator in Uysal (2022)
//Input: 
// Data & specifications
//     struct DS ds_struct;
//     ds_struct = dsCreate;
//     ds_struct = reshape(ds_struct,7,1);
//     ds_struct[1].dataMatrix = x; //n x k matrix of covariates including a constant
//     ds_struct[2].dataMatrix = y; //n x 1 vector of observed outcome 
//     ds_struct[3].dataMatrix = d; //n x 1 vector of treatment variable 
//     ds_struct[4].type = 1; //mf=1 or 2 3 (identiy, logit, poisson link)) only required for reg, dr1, dr2
//     ds_struct[5].type = 1; //pf=1 or 2 (logit or probit)) required for ipw1-ipw3, dr1, dr2
//     ds_struct[6].type = 1; //dr_wt = 1, 2, 3 for the weights for the dr1 corresponding to ipw1, ipw2, and ipw3 )
//     ds_struct[7].type = n; //sample size 
// Parameter Vector teta

// Output: 
// Asymptotic variance derived in the text
proc(1)= asyvar_reg(teta, struct DS ds_struct); 
    local AVb1, AVb11, AVb12, AVb0,
    AVb01, AVb02, psi1, psi2, Eetab1, Eetab0, AVtau;

    local k, b1,b0,mu1,mu0, xb1, xb0, my1, my0;
    local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
       N = ds_struct[7].Type;
    
        k=cols(indep);
b1=teta[1:k,.];
b0=teta[k+1:2*k,.];
mu1=teta[2*k+1,.];
mu0=teta[2*k+2,.];
    
//linear index
   xb1=indep*b1; 
   xb0=indep*b0;  
 if mf==1;
     my1=xb1;
     my0=xb0;
     AVb11= (1/rows(indep))*(treat.*indep)'indep;
     AVb01= (1/rows(indep))*((1-treat).*indep)'indep;
     Eetab1=meanc(indep);
     Eetab0=meanc(indep);
  elseif mf==2;
     my1=exp(xb1)./(1+exp(xb1));
     my0=exp(xb0)./(1+exp(xb0));  
     AVb11= (1/rows(indep))*(treat.*my1.*(1-my1).*indep)'*indep;
     AVb01= (1/rows(indep))*((1-treat).*my0.*(1-my0).*indep)'*indep;
     Eetab1=meanc(my1.*(1-my1).*indep);
     Eetab0=meanc(my0.*(1-my0).*indep);
  elseif mf==3;
     my1=exp(xb1);
     my0=exp(xb0);
     AVb11= (1/rows(indep))*(treat.*my1.*indep)'indep;
     AVb01= (1/rows(indep))*((1-treat).*my0.*indep)'indep;
     Eetab1=meanc(my1.*indep);
     Eetab0=meanc(my0.*indep);
 endif;      
 

    psi1=(indep.*(treat.*(dep-my1)));
    AVb12= (1/N)*psi1'*psi1;
    psi2=(indep.*((1-treat).*(dep-my0)));
    AVb02= (1/N)*psi2'*psi2;
 AVb1 =  inv(AVb11)*AVb12*inv(AVb11)';
 AVb0 =  inv(AVb01)*AVb02*inv(AVb01)';
 
 AVtau=meanc((my1-my0-(mu1-mu0))^2);
 retp(sqrt((1/N)*(AVtau+Eetab1'*AVb1*Eetab1+Eetab0'*AVb0*Eetab0)));
endp;
