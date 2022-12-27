
// Author: Derya Uysal
// Version: 12/27/2022
// 
//Function for the asmptotic variance of  the ipw3 estimator in Uysal (2022)
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
proc(1)= asyvar_ipw3(teta, struct DS ds_struct); 
        local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
       N = ds_struct[7].Type;
    local k, alp, mu1, mu0,eta1,eta0, psest, psi1, dpsest, Emu1alp, Emu0alp, EH, Vp,
     psi5, psi4, V44, V55, V45;
k=cols(indep);
alp=teta[1:k,.];
eta1=teta[k+1,.];
eta0=teta[k+2,.];
mu1=teta[k+3,.];
mu0=teta[k+4,.];
if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(treat-psest)); //Prob. function logit.
    Emu1alp=-meanc(((treat.*(dep-mu1+eta1))./(psest.^2)).*psest.*(1-psest).*indep);
    Emu0alp=meanc((((1-treat).*(dep-mu0+eta0))./((1-psest).^2)).*psest.*(1-psest).*indep);

elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = (indep.*dpsest.*(treat-psest)./(psest.*(1-psest))); 
    Emu1alp=-meanc(((treat.*(dep-mu1+eta1))./(psest.^2)).*dpsest.*indep);
    Emu0alp=meanc((((1-treat).*(dep-mu0+eta0))./((1-psest).^2)).*dpsest.*indep);

endif;
EH = -(1/N)*psi1'*psi1;

psi4=(treat.*(dep-mu1)./psest+eta1*((treat-psest)./psest));
psi5=((1-treat).*(dep-mu0)./(1-psest)-eta0*((treat-psest)./(1-psest)));

v44=(1/N)*psi4'*psi4;
v45=(1/N)*psi5'*psi4;
v55=(1/N)*psi5'*psi5;

/*Vp1 = meanc((treat.*(dep-mu1).^2)./psest.^2);
Vp2 = eta1.*meanc((treat.*(dep-mu1))./psest.^2);
Vp3 = meanc(((1-treat).*(dep-mu0).^2)./((1-psest).^2));
Vp4= eta0.*meanc(((1-treat).*(dep-mu0))./((1-psest).^2));
Vp = Vp1+Vp2+Vp3+Vp4+2*eta1*eta0;
*/
Vp=V44+V55-2*V45;
retp(sqrt((1/N)*(-(Emu1alp-Emu0alp)'*inv(-EH)*(Emu1alp-Emu0alp)+Vp)));
endp;
