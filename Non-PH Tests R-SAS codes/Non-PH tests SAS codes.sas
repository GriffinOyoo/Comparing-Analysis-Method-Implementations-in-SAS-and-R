

/*1. Weighted log-rank test.*/

/* Peto-Peto.*/
proc lifetest data=WHAS; 
time LENFOL*FSTAT(0); 
strata AFB / test=peto; 
run;

/*Modified Peto Peto*/
proc lifetest data=WHAS; 
time LENFOL*FSTAT(0); 
strata AFB / test=modpeto; 
run;

/*Tarone-Ware*/
proc lifetest data=WHAS; 
time LENFOL*FSTAT(0); 
strata AFB / test=taroneware; 
run;

/*Wilcoxon*/
proc lifetest data=wlr; 
time LENFOL*FSTAT(0); 
strata AFB / test= wilcoxon; 
run;

/*Fleming-Harrington*/
proc lifetest data=WHAS; 
time LENFOL*FSTAT(0); 
strata AFB / test=FH(0,1); 
run;


/*2. Max-Combo test*/

/*Max-Combo(Unstratified)*/

/*combo_wlr macro was used*/
%combo_wlr(
   data  = WHAS,
   group = AFB,
   time  = LENFOL,
   event = FSTAT
);

/*Stratified Max-Combo*/

%macro combo_wlr(data=,group=,time=,event=,weights=, strata=);

data one; 
 set &data;
 keep &group &time &event &strata; 
run;


proc lifetest data=one;
  time &time*&event(0);
  strata &strata/ group= &group test=FH(&r,&g);
 run;


%combo_wlr(
   data  = new,
   group = AFB_new,
   time  = LENFOL,
   event = FSTAT,
   strata= GENDER
);


/*3. Modestly Weighted Log-rank test.*/

/*Not SAS enabled for now.
