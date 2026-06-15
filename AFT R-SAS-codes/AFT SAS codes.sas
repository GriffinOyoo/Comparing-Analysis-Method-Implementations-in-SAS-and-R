/*THIS FILE CONTAINS SAS CODES FOR SIX AFT MODELS AND THEIR MODEL DIAGNOSTIC PLOTS*/



/*********************PROC LIFEREG PROCEDURE*************************/
/*1. LOG-LOGISTIC MODEL*/
proc lifereg data=AFT1;
model LENFOL*FSTAT(0) =  AGE GENDER HR BMI AFB/ dist=llogistic;
run;

/*2. LOG-NORMAL MODEL*/
proc lifereg data=AFT1;
model LENFOL*FSTAT(0) =  AGE GENDER HR BMI AFB/ dist=lognormal;
run;

/*3. WEIBULL MODEL*/
proc lifereg data=AFT1;
model LENFOL*FSTAT(0) = AGE GENDER HR BMI AFB/ dist=weibull;
run;

/*4. EXPONENTIAL MODEL*/
proc lifereg data=AFT1;
model LENFOL*FSTAT(0) = AGE GENDER HR BMI AFB/ dist=exponential;
run;


/*5. GENERALIZED GAMMA MODEL*/
proc lifereg data=AFT1;
model LENFOL*FSTAT(0) = AGE GENDER HR BMI AFB/ dist=gamma;
run;
/*6. GENERALIZED F MODEL*/ 
/*Not enabled



/*********************PARMSURV MACRO PROCEDURE*************************/
/*1. LOG-LOGISTIC MODEL*/
%paramsurv(data=AFT1, t1=LENFOL, censor=FSTAT, covars=AGE GENDER HR BMI AFB, dist=LLOGIS)

/*2. LOG-NORMAL MODEL*/
%paramsurv(data=AFT1, t1=LENFOL, censor=FSTAT, covars=AGE GENDER HR BMI AFB, dist=LNORM)

/*3. WEIBULL MODEL*/
%paramsurv(data=AFT1,  t1=LENFOL, censor=FSTAT, covars=AGE GENDER HR BMI AFB, dist=WEIBULL)

/*4. EXPONENTIAL MODEL*/
%paramsurv(data=AFT1, t1=LENFOL, censor=FSTAT, covars=AGE GENDER HR BMI AFB, dist=EXP)

/*5. GENERALIZED GAMMA MODEL*/
/*User-defined initial values*/
%paramsurv(data=AFT1, t1=LENFOL, censor=FSTAT, covars=AGE GENDER HR BMI AFB, dist=GENGAMMA, init={0 0 0 0 0 0 1 0.5}) 

/*6. GENERALIZED F MODEL*/
%paramsurv(data=AFT1, t1=LENFOL, censor=FSTAT, covars=AGE GENDER HR BMI AFB, dist=GENF)







/***********************MODEL DIAGNOSTIC PLOT****************************************/
/*Cox-Snell Residuals*/

/*LOG-LOGISTIC MODEL */
proc lifereg data=AFT1;
    class AFB GENDER;
    model LENFOL*FSTAT(0)= AGE GENDER HR BMI AFB / dist=loglogistic;
    output out=res_loglogistic cres=coxsnell;
run;
proc lifetest data=res_loglogistic outsurv=km_loglogistic plots=none;
    time coxsnell*FSTAT(0);
run;
data km_loglogistic;
    set km_loglogistic;
    observed   = -log(survival);
    lowerlimit = -log(sdf_ucl);
    upperlimit = -log(sdf_lcl);
    expected   = coxsnell;
    model="Log-Logistic";
run;


/*LOG-NORMAL MODEL*/
proc lifereg data=AFT1;
    class AFB GENDER;
    model LENFOL*FSTAT(0)= AGE GENDER HR BMI AFB / dist=lnormal;
    output out=res_lognormal cres=coxsnell;
run;

proc lifetest data=res_lognormal outsurv=km_lognormal plots=none;
    time coxsnell*FSTAT(0);
run;

data km_lognormal;
    set km_lognormal;
    observed   = -log(survival);
    lowerlimit = -log(sdf_ucl);
    upperlimit = -log(sdf_lcl);
    expected   = coxsnell;
    model="Lognormal";
run;


/*  WEIBULL MODEL */
proc lifereg data=AFT1;
    class AFB GENDER;
    model LENFOL*FSTAT(0)= AGE GENDER HR BMI AFB / dist=weibull;
    output out=res_weibull cres=coxsnell;
run;

proc lifetest data=res_weibull outsurv=km_weibull plots=none;
    time coxsnell*FSTAT(0);
run;

data km_weibull;
    set km_weibull;
    observed   = -log(survival);
    lowerlimit = -log(sdf_ucl);
    upperlimit = -log(sdf_lcl);
    expected   = coxsnell;
    model="Weibull";
run;


/*EXPONENTIAL MODEL*/
proc lifereg data=AFT1;
    class AFB GENDER;
    model LENFOL*FSTAT(0)= AGE GENDER HR BMI AFB / dist=exponential;
    output out=res_exponential cres=coxsnell;
run;

proc lifetest data=res_exponential outsurv=km_exponential plots=none;
    time coxsnell*FSTAT(0);
run;

data km_exponential;
    set km_exponential;
    observed   = -log(survival);
    lowerlimit = -log(sdf_ucl);
    upperlimit = -log(sdf_lcl);
    expected   = coxsnell;
    model="Exponential";
run;


/* GAMMA MODEL */
proc lifereg data=AFT1;
    class AFB GENDER;
    model LENFOL*FSTAT(0)= AGE GENDER HR BMI AFB / dist=gamma;
    output out=res_gamma cres=coxsnell;
run;

proc lifetest data=res_gamma outsurv=km_gamma plots=none;
    time coxsnell*FSTAT(0);
run;

data km_gamma;
    set km_gamma;
    observed   = -log(survival);
    lowerlimit = -log(sdf_ucl);
    upperlimit = -log(sdf_lcl);
    expected   = coxsnell;
    model="Gamma";
run;


/*COMBINE ALL MODELS*/
data km_all;
    set km_weibull
        km_loglogistic
        km_lognormal  
        km_gamma  
        km_exponential;
run;


/*2×3 PANEL PLOT*/
proc sgpanel data=km_all;
    panelby model / columns=3 rows=2;
    /* Step cumulative hazard */
    series x=coxsnell y=observed /
        lineattrs=(color=black thickness=2);
    /* Confidence region */
    band x=coxsnell lower=lowerlimit upper=upperlimit /
        fillattrs=(color=gray transparency=0.7);
    /* Dashed confidence limits */
    series x=coxsnell y=lowerlimit /
        lineattrs=(color=gray pattern=dash);
    series x=coxsnell y=upperlimit /
        lineattrs=(color=gray pattern=dash);
    /* Exponential(1) reference line */
    series x=coxsnell y=expected /
        lineattrs=(color=red thickness=2);
    colaxis label="Cox–Snell residual";
    rowaxis label="Cumulative hazard of residuals";
run;
