// Author: Derya Uysal
// Version: 12/27/2022
// 
//Procedure to estimate the ATE using eqSolvemt with user written moment functions
//Input: f and corresponding f2:
//       ipw1_mom & ipw1_crossmom (requires pf=1 or 2 (logit or probit))
//       ipw2_mom & ipw2_crossmom (requires pf=1 or 2 (logit or probit)) 
//       ipw3_mom & ipw3_crossmom(requires pf=1 or 2 (logit or probit)) 
//       reg_mom  & reg_crossmom (requires mf=1 or 2 3 (identiy, logit, poisson link)) 
//       reg_dr1  & dr1_crossmom (requires mf=1 or 2 3 (idetiy, logit, poisson link), 
//               (requires pf=1 or 2 (logit or probit)) and
//                 dr_wt = 1, 2, 3 for the weights corresponding to ipw1, ipw2, and ipw3 )
//       reg_dr2  & dr2_crossmom (requires mf=1 or 2 3 (identiy, logit, poisson link), 
//               (requires pf=1 or 2 (logit or probit)) 
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

//Output: 
//        outate.ate  = estimated treatment effect
//        outate.atese  = se of the treatment effect
//        outate.pars=   estimated coefficients (changes depending on the method)
//        outate.parsse = standard errors of the coefficients
//        outate.retcode = return code same as the retcode of eqSolvemt



proc(1) =  mm_ate(&f, &f2, struct DS ds_struct); 
        library dc, maxlik;
    
#include maxlik.ext
local indep, dep, treat, mf, pf, dr_wt, N, v_mest, se_mest;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
   N = ds_struct[7].Type;
 //Data, mean spec. linear: mf==1, logit: 2 poisson: 3
    local alp_start, estcoef,ate,y1,y0,x1,x0,b1_start,b0_start,mmi1,mmi0,mmi;
    local pshat, eta1s, eta0s;
    local f:proc; //one of the 8: ipw1, ipw2, ipw3, dr1a, dr1b,dr1c, dr2, reg
    local f2:proc;
    y1=selif(dep,treat.==1);
    y0=selif(dep,treat.==0);
    x1=selif(indep,treat.==1);
    x0=selif(indep,treat.==0);    
    trap 1;

// Attempt to compute the inverse of the moment matrix
mmi1 = inv(x1'x1);

// Check to see if 'mmi', contains a scalar error code
if scalmiss(mmi1);
    // Compute the pseudo-inverse of the moment matrix
    mmi1 = pinv(x1'x1);
endif;

// Solve the linear equations
    b1_start=mmi1*x1'y1;
    trap 1;

// Attempt to compute the inverse of the moment matrix
mmi0 = inv(x0'x0);

// Check to see if 'mmi', contains a scalar error code
if scalmiss(mmi0);
    // Compute the pseudo-inverse of the moment matrix
    mmi0 = pinv(x0'x0);
endif;

    b0_start=mmi0*x0'y0;
    trap 1;

// Attempt to compute the inverse of the moment matrix
mmi = inv(indep'indep);

// Check to see if 'mmi', contains a scalar error code
if scalmiss(mmi);
    // Compute the pseudo-inverse of the moment matrix
    mmi = pinv(indep'indep);
endif;
maxset;
        _max_CovPar = 0;
    _max_GradMethod = 1;
    __output = 0;
_max_GradProc = &llogitgrad;
_max_MaxIters = 1000;
    
    //starting values
    local alp_starts, nix1,nix2,nix3, retips;
     alp_starts=mmi*indep'treat;

 
// Step four: Call binary logit procedure
    if pf==1;

    {alp_start,nix1,nix2,nix3, retips} = maxlik(treat~indep~ones(rows(treat),1),0,&llogit,alp_starts);             


   /* elseif pf==2;
        dcout1 = binaryProbit(dcCt);*/
    endif;
//    alp_start=pvGetParVector(dcout1.par);
if pf==1;//if logit
pshat=exp((indep)*alp_start)./(1+exp((indep)*alp_start)); //Prob. function logit.
elseif pf==2;//if probit
pshat=cdfn((indep)*alp_start); //Prob. function probit.
endif;    


struct PV par;
par=pvCreate();

//Starting values for each method:
//REG;
local y1hat, y0hat, mu1reg, mu0reg, mu1ipw1, mu0ipw1, mu1ipw2, mu0ipw2,mu1ipw3,mu0ipw3,omega1,omega0,c1ipw3,c0ipw3;
y1hat = indep*b1_start;
y0hat = indep*b0_start;

mu1reg= meanc(y1hat);
mu0reg= meanc(y0hat);
//IPW1

mu1ipw1= meanc(dep.*treat./pshat);
mu0ipw1= meanc(dep.*(1-treat)./(1-pshat));

//IPW2
omega1=(treat./pshat)./meanc(treat./pshat);
omega0=((1-treat)./(1-pshat))./meanc((1-treat)./(1-pshat));
mu1ipw2= meanc(dep.*omega1);
mu0ipw2= meanc(dep.*omega0);
//IPW3

c1ipw3 = meanc((treat-pshat)./pshat)/meanc(((treat-pshat)./pshat).^2);
c0ipw3 = -meanc((treat-pshat)./(1-pshat))/meanc(((treat-pshat)./(1-pshat)).^2);
mu1ipw3 = meanc((dep.*treat./pshat).*(1-c1ipw3./pshat))/meanc((treat./pshat).*(1-c1ipw3./pshat));
mu0ipw3 = meanc((dep.*(1-treat)./(1-pshat)).*(1-c0ipw3./(1-pshat)))/meanc(((1-treat)./(1-pshat)).*(1-c0ipw3./(1-pshat)));

eta1s=-meanc(treat.*(dep-mu1ipw3)./(pshat.^2))./meanc(((treat-pshat)./pshat).^2);
eta0s=-meanc((1-treat).*(dep-mu0ipw3)./((1-pshat).^2))./meanc(((treat-pshat)./(1-pshat)).^2);

    //Regression
    if &f==&reg_mom;
par= pvPack(par,b1_start, "b1");
par= pvPack(par, b0_start, "b0");
par= pvPack(par, mu1reg, "mu1");
par= pvPack(par, mu0reg, "mu0");    
     //IPW 1   
     elseif &f==&ipw1_mom;
par= pvPack(par, alp_start, "alpha");
par= pvPack(par, mu1ipw1, "mu1");
par= pvPack(par, mu0ipw1, "mu0");            
      //IPW 2   
     elseif &f==&ipw2_mom;
par= pvPack(par,alp_start, "alpha");
par= pvPack(par, mu1ipw2, "mu1");
par= pvPack(par, mu0ipw2, "mu0");     
          
          
     //IPW 3  
     elseif &f==&ipw3_mom;
par= pvPack(par,alp_start, "alpha");
par= pvPack(par, eta1s, "eta1");
par= pvPack(par, eta0s, "eta0");
par= pvPack(par, mu1ipw3, "mu1");
par= pvPack(par, mu0ipw3, "mu0");     


             
         //DR1  
     elseif &f==&dr1_mom;  
par= pvPack(par,alp_start, "alpha"); 
par= pvPack(par,b1_start, "b1");
par= pvPack(par, b0_start, "b0");     
par= pvPack(par, mu1reg, "mu1");
par= pvPack(par, mu0reg, "mu0");           
            //DR2
     elseif &f==&dr2_mom; 
par= pvPack(par,alp_start, "alpha");         
par= pvPack(par,b1_start, "b1");
par= pvPack(par, b0_start, "b0");     
par= pvPack(par, mu1reg, "mu1");
par= pvPack(par, mu0reg, "mu0");       
 endif; 

// Declare control structure and fill with defaults
struct eqSolvemtControl c;
c = eqSolvemtControlCreate();
c.output = 0;
 c.maxIters = 1000;
// Declare output structure to hold results
struct eqSolvemtOut out;
// Solve the system of equations
out =  eqSolvemt(&f, par, ds_struct, c);

    ate=pvUnpack(out.par, "mu1")-pvUnpack(out.par, "mu0");
    estcoef=pvGetParVector(out.par);

    v_mest=(1/N)*inv(gradMT(&f,out.par, ds_struct))*f2(out.par, ds_struct)*inv(gradMT(&f,out.par, ds_struct))';
    se_mest=sqrt(v_mest[rows(estcoef)-1,rows(estcoef)-1]+v_mest[rows(estcoef),rows(estcoef)]-2*v_mest[rows(estcoef)-1,rows(estcoef)]);
    //Define structure

   struct OutputsATE outate;
   outate.ate  = ate;
   outate.atese  = se_mest;
   outate.pars=  estcoef;
   outate.parsse = sqrt(diag(v_mest));
   outate.retcode = out.retcode;

    retp(outate);
 
 
    //retp(ate, se_mest, estcoef, sqrt(diag(v_mest)), out.par,out.retcode);


endp;
    



proc reg_crossmom(struct PV p, struct DS ds_struct);
      local xb1, xb0, my1, my0,psi1, beta1, beta0, mu1, mu0, psi2, psi3, psi4;
     beta1 = pvUnpack(p, "b1");
     beta0 = pvUnpack(p, "b0");
     mu1 = pvUnpack(p, "mu1");
     mu0 = pvUnpack(p, "mu0");
    local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
   N = ds_struct[7].Type;
local v11,v12,v13,v14,v21,v22,v23,v24;
local v31,v32,v33,v34,v41,v42,v43,v44;


//linear index
  xb1=indep*beta1; 
   xb0=indep*beta0;  
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
psi1=(indep.*(treat.*(dep-my1)));
psi2=(indep.*((1-treat).*(dep-my0)));
psi3=(my1)-mu1;
psi4=(my0)-mu0;

v11=(1/N)*psi1'*psi1;
v21=(1/N)*psi2'*psi1;
v31=(1/N)*psi3'*psi1;
v41=(1/N)*psi4'*psi1;

v12=(1/N)*psi1'*psi2;
v22=(1/N)*psi2'*psi2;
v32=(1/N)*psi3'*psi2;
v42=(1/N)*psi4'*psi2;

v13=(1/N)*psi1'*psi3;
v23=(1/N)*psi2'*psi3;
v33=(1/N)*psi3'*psi3;
v43=(1/N)*psi4'*psi3;
v14=(1/N)*psi1'*psi4;
 
v24=(1/N)*psi2'*psi4;
v34=(1/N)*psi3'*psi4;
v44=(1/N)*psi4'*psi4;


retp ((v11~v12~v13~v14)|(v21~v22~v23~v24)|(v31~v32~v33~v34)|(v41~v42~v43~v44));

endp;



proc ipw1_crossmom(struct PV p, struct DS ds_struct);
local psi1, alp, mu1, mu0, psi2, psi3;
     alp = pvUnpack(p, "alpha");
     mu1 = pvUnpack(p, "mu1");
     mu0 = pvUnpack(p, "mu0");

    local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
   N = ds_struct[7].Type;
local   psest, dpsest,
v11,v12,v13,v21,v22,v23,v31,v32,v33;

if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(treat-psest)); //Prob. function logit.
elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = (indep.*dpsest.*(treat-psest)./(psest.*(1-psest))); 
endif;

psi2=((treat./psest).*dep)-mu1;
psi3=(((1-treat)./(1-psest)).*dep)-mu0;

v11=(1/N)*psi1'*psi1;
v21=(1/N)*psi2'*psi1;
v31=(1/N)*psi3'*psi1;

v12=(1/N)*psi1'*psi2;
v22=(1/N)*psi2'*psi2;
v32=(1/N)*psi3'*psi2;

v13=(1/N)*psi1'*psi3;
v23=(1/N)*psi2'*psi3;
v33=(1/N)*psi3'*psi3;


retp ((v11~v12~v13)|(v21~v22~v23)|(v31~v32~v33));

endp;


proc ipw2_crossmom(struct PV p, struct DS ds_struct);
local    psi1, alp, mu1, mu0, psi2, psi3;
     alp = pvUnpack(p, "alpha");
     mu1 = pvUnpack(p, "mu1");
     mu0 = pvUnpack(p, "mu0");

    local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
   N = ds_struct[7].Type;
local   psest, dpsest,
v11,v12,v13,v21,v22,v23,v31,v32,v33;



if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(treat-psest)); //Prob. function logit.
elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = (indep.*dpsest.*(treat-psest)./(psest.*(1-psest))); 
endif;

psi2=(treat.*(dep-mu1)./psest);
psi3=((1-treat).*(dep-mu0)./(1-psest));

v11=(1/N)*psi1'*psi1;
v21=(1/N)*psi2'*psi1;
v31=(1/N)*psi3'*psi1;

v12=(1/N)*psi1'*psi2;
v22=(1/N)*psi2'*psi2;
v32=(1/N)*psi3'*psi2;

v13=(1/N)*psi1'*psi3;
v23=(1/N)*psi2'*psi3;
v33=(1/N)*psi3'*psi3;


retp ((v11~v12~v13)|(v21~v22~v23)|(v31~v32~v33));

endp;


proc ipw3_crossmom(struct PV p, struct DS ds_struct);
local    psi1, alp, mu1, mu0, eta1, eta0, psi2, psi3,psi4,psi5;
     alp = pvUnpack(p, "alpha");
eta1 = pvUnpack(p, "eta1");
eta0 = pvUnpack(p, "eta0"); 
    mu1 = pvUnpack(p, "mu1");
     mu0 = pvUnpack(p, "mu0");

    local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
   N = ds_struct[7].Type;


local psest,dpsest,
v11,v12,v13,v14,v15,v21,v22,v23,v24,v25,v31,v32,v33,v34,v35,v41,v42,v43,v44,v45,v51,v52,v53,v54,v55;




if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(treat-psest)); //Prob. function logit.
elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = (indep.*dpsest.*(treat-psest)./(psest.*(1-psest))); 
endif;

psi2=(treat.*(dep-mu1)./(psest.^2)+eta1.*((treat-psest)./psest).^2);
psi3=((1-treat).*(dep-mu0)./((1-psest).^2)+eta0.*((treat-psest)./(1-psest)).^2);
psi4=(treat.*(dep-mu1)./psest+eta1*((treat-psest)./psest));
psi5=((1-treat).*(dep-mu0)./(1-psest)-eta0*((treat-psest)./(1-psest)));

v11=(1/N)*psi1'*psi1;
v21=(1/N)*psi2'*psi1;
v31=(1/N)*psi3'*psi1;
v41=(1/N)*psi4'*psi1;
v51=(1/N)*psi5'*psi1;

v12=(1/N)*psi1'*psi2;
v22=(1/N)*psi2'*psi2;
v32=(1/N)*psi3'*psi2;
v42=(1/N)*psi4'*psi2;
v52=(1/N)*psi5'*psi2;


v13=(1/N)*psi1'*psi3;
v23=(1/N)*psi2'*psi3;
v33=(1/N)*psi3'*psi3;
v43=(1/N)*psi4'*psi3;
v53=(1/N)*psi5'*psi3;

v14=(1/N)*psi1'*psi4;
v24=(1/N)*psi2'*psi4;
v34=(1/N)*psi3'*psi4;
v44=(1/N)*psi4'*psi4;
v54=(1/N)*psi5'*psi4;

v15=(1/N)*psi1'*psi5;
v25=(1/N)*psi2'*psi5;
v35=(1/N)*psi3'*psi5;
v45=(1/N)*psi4'*psi5;
v55=(1/N)*psi5'*psi5;

retp ((v11~v12~v13~v14~v15)|(v21~v22~v23~v24~v25)|(v31~v32~v33~v34~v35)|(v41~v42~v43~v44~v45)|(v51~v52~v53~v54~v55));

endp;


proc dr1_crossmom(struct PV p, struct DS ds_struct);
local alp,b1,b0, omega1,omega0, mu1, mu0,  psi1,psi2,psi3,psi4, psi5,psest, dpsest, A,B,
    c0, c1, xb1, xb0, my1, my0;
     alp = pvUnpack(p, "alpha");
     mu1 = pvUnpack(p, "mu1");
     mu0 = pvUnpack(p, "mu0");
     b1 = pvUnpack(p, "b1");
     b0 = pvUnpack(p, "b0");
    local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
   N = ds_struct[7].Type;
local v11,v12,v13,v14,v15,v21,v22,v23,v24,v25,v31,v32,v33,v34,v35,v41,v42,v43,v44,v45,v51,v52,v53,v54,v55;


//Propensity score specification
if pf==1;
psest=exp((indep)*alp)./(1+exp((indep)*alp)); //Prob. function logit.
psi1 = (indep.*(treat-psest)); //Prob. function logit.
elseif pf==2;
psest=cdfn((indep)*alp); //Prob. function probit.
dpsest=pdfn((indep)*alp); //density. function probit.
psi1 = (indep.*dpsest.*(treat-psest)./(psest.*(1-psest))); 
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
 psi2=(indep.*(treat.*(dep-my1)));
psi3=(indep.*((1-treat).*(dep-my0)));
psi4=(omega1.*(dep-my1)+my1-mu1);
psi5=(omega0.*(dep-my0)+my0-mu0);



v11=(1/N)*psi1'*psi1;
v21=(1/N)*psi2'*psi1;
v31=(1/N)*psi3'*psi1;
v41=(1/N)*psi4'*psi1;
v51=(1/N)*psi5'*psi1;

v12=(1/N)*psi1'*psi2;
v22=(1/N)*psi2'*psi2;
v32=(1/N)*psi3'*psi2;
v42=(1/N)*psi4'*psi2;
v52=(1/N)*psi5'*psi2;


v13=(1/N)*psi1'*psi3;
v23=(1/N)*psi2'*psi3;
v33=(1/N)*psi3'*psi3;
v43=(1/N)*psi4'*psi3;
v53=(1/N)*psi5'*psi3;

v14=(1/N)*psi1'*psi4;
v24=(1/N)*psi2'*psi4;
v34=(1/N)*psi3'*psi4;
v44=(1/N)*psi4'*psi4;
v54=(1/N)*psi5'*psi4;

v15=(1/N)*psi1'*psi5;
v25=(1/N)*psi2'*psi5;
v35=(1/N)*psi3'*psi5;
v45=(1/N)*psi4'*psi5;
v55=(1/N)*psi5'*psi5;

retp ((v11~v12~v13~v14~v15)|(v21~v22~v23~v24~v25)|(v31~v32~v33~v34~v35)|(v41~v42~v43~v44~v45)|(v51~v52~v53~v54~v55));


endp;

proc dr2_crossmom(struct PV p, struct DS ds_struct);
    local alp,b1,b0, mu1, mu0;
     alp = pvUnpack(p, "alpha");
     mu1 = pvUnpack(p, "mu1");
     mu0 = pvUnpack(p, "mu0");
     b1 = pvUnpack(p, "b1");
     b0 = pvUnpack(p, "b0");
    local indep, dep, treat, mf, pf, dr_wt, N;
   indep = ds_struct[1].dataMatrix;
   dep = ds_struct[2].dataMatrix;
   treat = ds_struct[3].dataMatrix;
   mf = ds_struct[4].Type;
   pf = ds_struct[5].Type;
   dr_wt = ds_struct[6].Type;
   N = ds_struct[7].Type;
local psi1,psi2,psi3,psi4, psi5,psest,dpsest, xb1, xb0, my1,my0,
v11,v12,v13,v14,v15,v21,v22,v23,v24,v25,v31,v32,v33,v34,v35,v41,v42,v43,v44,v45,v51,v52,v53,v54,v55;


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
  elseif mf==2;
     my1=exp(xb1)./(1+exp(xb1));
     my0=exp(xb0)./(1+exp(xb0));   
  elseif mf==3;
     my1=exp(xb1);
     my0=exp(xb0);
 endif; 


psi2=(indep.*((treat./psest).*(dep-my1)));
psi3=(indep.*(((1-treat)./(1-psest)).*(dep-my0)));
psi4=(my1)-mu1;
psi5=(my0)-mu0;
v11=(1/N)*psi1'*psi1;
v21=(1/N)*psi2'*psi1;
v31=(1/N)*psi3'*psi1;
v41=(1/N)*psi4'*psi1;
v51=(1/N)*psi5'*psi1;

v12=(1/N)*psi1'*psi2;
v22=(1/N)*psi2'*psi2;
v32=(1/N)*psi3'*psi2;
v42=(1/N)*psi4'*psi2;
v52=(1/N)*psi5'*psi2;


v13=(1/N)*psi1'*psi3;
v23=(1/N)*psi2'*psi3;
v33=(1/N)*psi3'*psi3;
v43=(1/N)*psi4'*psi3;
v53=(1/N)*psi5'*psi3;

v14=(1/N)*psi1'*psi4;
v24=(1/N)*psi2'*psi4;
v34=(1/N)*psi3'*psi4;
v44=(1/N)*psi4'*psi4;
v54=(1/N)*psi5'*psi4;

v15=(1/N)*psi1'*psi5;
v25=(1/N)*psi2'*psi5;
v35=(1/N)*psi3'*psi5;
v45=(1/N)*psi4'*psi5;
v55=(1/N)*psi5'*psi5;

retp ((v11~v12~v13~v14~v15)|(v21~v22~v23~v24~v25)|(v31~v32~v33~v34~v35)|(v41~v42~v43~v44~v45)|(v51~v52~v53~v54~v55));



endp;





/*Logit likelihood*/
proc llogit(b,z);//b coeff, z -> data z[.,1]: dependent, z[2:cols(z)-1] independent without a constant, z[.,cols(z)]: weights
    local m, w, cdfl, ll;
    if meanc(z[.,2]) .== 1;
    m = (z[.,2:cols(z)-1])*b;
    else;
    m = (ones(rows(z),1)~z[.,2:cols(z)-1])*b;
    endif;
    w = z[., cols(z)];
    cdfl= exp(m)./(1+exp(m));
    ll= w.*(z[.,1].*ln(cdfl)+(1-z[.,1]).*ln(1-cdfl));
    retp(ll);
endp;

proc llogitgrad(b,z);//b coeff, z -> data z[.,1]: dependent, z[2:cols(z)-1] independent without a constant, z[.,cols(z)]: weights
    local m, w, cdfl, gll;
    w = z[., cols(z)];
    if meanc(z[.,2]) .== 1;
    m = (z[.,2:cols(z)-1])*b;   
    cdfl= exp(m)./(1+exp(m));
    gll= w.*(z[.,1]-cdfl).*(z[.,2:cols(z)-1]);
    else;
    m = (ones(rows(z),1)~z[.,2:cols(z)-1])*b;
    cdfl= exp(m)./(1+exp(m));
    gll= w.*(z[.,1]-cdfl).*(ones(rows(z),1)~z[.,2:cols(z)-1]);        
    endif;
    retp(gll);
endp;



