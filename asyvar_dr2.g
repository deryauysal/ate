
// Author: Derya Uysal
// Version: 12/27/2022
// 
//Function for the asmptotic variance of  the second DR estimator in Uysal (2022)
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

proc(1)= asyvar_dr2(teta, struct DS ds_struct); 
        local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
       N = ds_struct[7].Type;
    local k, alp, b1, b0, mu1, mu0;
    local psest, psi1,dpsest;
    local xb1, xb0, my1, my0;
    local Emu1b1, Emu0b0, Em1b1, Em0b0;
    local psi2, psi3, psi4, psi5, psi11, psi22, psi33,psi21, psi31;
    local V, V44, V55, V45;
k=cols(indep);
alp=teta[1:k,.];
b1=teta[k+1:2*k,.];
b0=teta[2*k+1:3*k,.];
mu1=teta[3*k+1,.];
mu0=teta[3*k+2,.];

//Propensity score specification
if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(treat-psest)); //Prob. function logit.
elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = (indep.*dpsest.*(treat-psest)./(psest.*(1-psest))); 
endif;

//Conditional mean of potential outcome
xb1=indep*b1; 
xb0=indep*b0;
 if mf==1;
     my1=xb1;
     my0=xb0;
     Emu1b1=meanc(indep)';
     Emu0b0=meanc(indep)';
     Em1b1=-(1/N)*indep'indep;
     Em0b0=-(1/N)*indep'indep;
  elseif mf==2;
     my1=exp(xb1)./(1+exp(xb1));
     my0=exp(xb0)./(1+exp(xb0));  
     Emu1b1=meanc(my1.*(1-my1).*indep)';
     Emu0b0=meanc(my0.*(1-my0).*indep)';
     Em1b1=-(1/N)*((my0.*(1-my0).*indep)'indep);
     Em0b0=-(1/N)*((my0.*(1-my0).*indep)'indep);
  elseif mf==3;
     my1=exp(xb1);
     my0=exp(xb0);
     Emu1b1=meanc(my1.*indep)';
     Emu0b0=meanc(my0.*indep)';
     Em1b1=-(1/N)*(my1.*indep)'indep;
     Em0b0=-(1/N)*(my1.*indep)'indep;
 endif;      
 
psi2=(indep.*((treat./psest).*(dep-my1)));
psi3=(indep.*(((1-treat)./(1-psest)).*(dep-my0)));
psi4=(my1)-mu1;
psi5=(my0)-mu0;
psi22= (1/N)*psi2'*psi2;
psi33= (1/N)*psi3'*psi3;
psi21= (1/N)*psi2'*psi1;
psi31= (1/N)*psi3'*psi1;
psi11= (1/N)*psi1'*psi1;

 
 V44=Emu1b1*inv(Em1b1)*psi22*inv(Em1b1)'*Emu1b1'-Emu1b1*inv(Em1b1)*psi21'*inv(psi11)*psi21*inv(Em1b1)'*Emu1b1'+meanc(psi4.^2);
 V55=Emu0b0*inv(Em0b0)*psi33*inv(Em0b0)'*Emu0b0'-Emu0b0*inv(Em0b0)*psi31'*inv(psi11)*psi31*inv(Em0b0)'*Emu0b0' + meanc(psi5.^2);
 V45=Emu1b1*inv(Em1b1)*(-psi21)*inv(psi11)*(-psi31)*inv(Em0b0)'*Emu0b0'+ meanc(psi4.*psi5);
 
 V=V44+V55-2*V45;
 
// Emu1b1*inv(Em1b1)*psi22*inv(Em1b1)'*Emu1b1'+meanc(psi4.^2)+Emu0b0*inv(Em0b0)*psi33*inv(Em0b0)'*Emu0b0'+meanc(psi5.^2);

//  (Emu1b1*inv(Em1b1)*psi21' + Emu0b0*inv(Em0b0)*psi31')*inv(psi11)*(Emu1b1*inv(Em1b1)*psi21' + Emu0b0*inv(Em0b0)*psi31')'
 retp(sqrt((1/N)*V));
endp;
