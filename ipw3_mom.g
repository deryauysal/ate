
// Author: Derya Uysal
// Version: 12/27/2022
// 
//Moment functions for the ipw3 estimator in Uysal (2022)
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
//     par= pvPack(par,alp_start, "alpha");
//     par= pvPack(par, eta1s, "eta1");
//     par= pvPack(par, eta0s, "eta0");
//     par= pvPack(par, mu1ipw3, "mu1");
//     par= pvPack(par, mu0ipw3, "mu0");     
// Output: 
// Moment functions


proc ipw3_mom(struct PV p, struct DS ds_struct); 
local indep, dep, treat, mf, pf, dr_wt;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
local alp, mu1, mu0,  psi1,psi2,psi3,psi4,psi5,eta1,eta0,psest,dpsest;
     alp = pvUnpack(p, "alpha");
eta1 = pvUnpack(p, "eta1");
eta0 = pvUnpack(p, "eta0");
     mu1 = pvUnpack(p, "mu1");
     mu0 = pvUnpack(p, "mu0");
if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = meanc(indep.*(treat-psest)); //Prob. function logit.
elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = meanc(indep.*dpsest.*(treat-psest)./(psest.*(1-psest))); 
endif;

psi2=meanc(treat.*(dep-mu1)./(psest.^2)+eta1.*((treat-psest)./psest).^2);
psi3=meanc((1-treat).*(dep-mu0)./((1-psest).^2)+eta0.*((treat-psest)./(1-psest)).^2);
psi4=meanc(treat.*(dep-mu1)./psest+eta1*((treat-psest)./psest));
psi5=meanc((1-treat).*(dep-mu0)./(1-psest)-eta0*((treat-psest)./(1-psest)));

retp (psi1|psi2|psi3|psi4|psi5);

endp;
