/*PROC FREQ*/


/*H0: P=0.5*/

/*1. Exact Binomial Test*/
proc freq data=coin_flips;
    tables result / binomial(p=0.5 );
    exact binomial;
run;


/*2. Wald Binomial Test(Asymptotic)*/
proc freq data=coin_flips;
    tables result / binomial(p=0.5);
run;

/*3. Mid-p adjusted Test*/
proc freq data=coin_flips;     
    tables result / binomial(p=0.5 level='H' cl=midp );  
    exact binomial / midp;
run;

/*4. Wilson Score Test*/
proc freq data=coin_flips;
    tables result / binomial(level='H' p=0.5 cl=score);
	run;





/*H0: P=0.19*/


	/*1. Exact Binomial Test*/
proc freq data=lung_cancer;
    tables death_flag / binomial(p=0.19 level='1');
    exact binomial;
    title "Exact Binomial Test for Death Proportion";
run;

/*2. Wald Binomial Test (Asymptotic)*/
proc freq data=lung_cancer;
    tables death_flag / binomial(p=0.19 level='1');
    title "Asymptotic Binomial Test for Death Proportion"; 
run;

/*3. Mid-p adjusted Exact Binomial Test*/
proc freq data=lung_cancer;     
    tables death_flag / binomial(p=0.19 level='1');  
    exact binomial / midp;
    title "Exact mid-P Binomial Test for Death Proportion"; 
run;

/*4. Wilson Score Test.*/
proc freq data=lung_cancer;
    tables death_flag / binomial(level='1' p=0.19 cl=score);
    title "Wilson Score Test for Death Proportion";
run;
