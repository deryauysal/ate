
// Author: Derya Uysal
// Version: 12/27/2022
// 
//Moment functions for the first DR estimator in Uysal (2022)
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
// Parameter Vector
//     struct PV p
//      par= pvPack(p,alp_start, "alpha");         
//      par= pvPack(p,b1_start, "b1");
//      par= pvPack(p, b0_start, "b0");     
//      par= pvPack(p, mu1reg, "mu1");
//      par= pvPack(p, mu0reg, "mu0");   
// Output: 
// Moment functions


proc dr1_mom(struct PV p, struct DS ds_struct); 
local indep, dep, treat, mf, pf, dr_wt;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
//with three different weighting functions. dr_wt=1, 2,3
//two differen ps spec
//three differen mean y

local alp,b1,b0, omega1,omega0, mu1, mu0,  psi1,psi2,psi3,psi4, psi5,psest, dpsest, A,B,
    c0, c1, xb1, xb0, my1, my0;
     alp = pvUnpack(p, "alpha");
     b1 = pvUnpack(p, "b1");
     b0 = pvUnpack(p, "b0");
     mu1 = pvUnpack(p, "mu1");
     mu0 = pvUnpack(p, "mu0");
//Propensity score specification
if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = meanc(indep.*(treat-psest)); //Prob. function logit.
elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = meanc(indep.*dpsest.*(treat-psest)./(psest.*(1-psest))); 
endif;

//Weighting type
if dr_wt==1;
omega1=treat./psest;
omega0=(1-treat)./(1-psest);
elseif dr_wt==2;
omega1=(treat./psest)./meanc(treat./psest);
omega0=((1-treat)./(1-psest))./meanc((1-treat)./(1-psest));
elseif dr_wt==3;
A=(1-treat)./(1-psest);
B=(treat./psest);
c0=((1/(1-psest)).*meanc(A.*psest-treat))./meanc((A.*psest-treat)^2);
c1=((1/(psest)).*meanc(B.*(1-psest)-(1-treat)))./meanc((B.*(1-psest)-(1-treat))^2);
omega1=((treat./psest).*(1-c1))./(meanc((treat./psest).*(1-c1)));
omega0=(((1-treat)./(1-psest)).*(1-c0))./(meanc(((1-treat)./(1-psest)).*(1-c0)));    
endif;

//Conditional mean of potential outcome
xb1=indep*b1; 
xb0=indep*b0;
 if mf==1;
     my1=xb1;
     my0=xb0;
  elseif mf==2;
     my1=exp(xb1)./(1+exp(xb1));
     my0=exp(xb0)./(1+exp(xb0));   
  elseif mf==3;
     my1=exp(xb1);
     my0=exp(xb0);
 endif;      
psi2=meanc(indep.*(treat.*(dep-my1)));
psi3=meanc(indep.*((1-treat).*(dep-my0)));
psi4=meanc(omega1.*(dep-my1)+my1-mu1);
psi5=meanc(omega0.*(dep-my0)+my0-mu0);


retp (psi1|psi2|psi3|psi4|psi5);

endp;
