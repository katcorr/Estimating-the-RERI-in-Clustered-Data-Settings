/************************************************************************************************************
PROJECT:	simulation study for estimating RERI in clustered data setting
PROGRAM:	run_sim_functions_v01.sas
PROGRAMMER:     kat
DATE:		2/7/2017	
RE:		run macros for the various procedures we'll be comparing 

NOTE:		From SAS Technical Support 8/27/2015:
		"GLIMMIX does not produce a covariance matrix involving the fixed effect parameters and the elements of the
		 variance/covariance matrix. The MMEQ computations could be resource intensive computation in terms of memory and time so the
		 developers have have been hesitant to add an option to the GLIMMIX procedure.  I will let them know of your request for such 
		an option.

		For now, only PROC NLMIXED can produce the cross-covariances between the fixed and random effects.  If your GLIMMIX model is not
		too complicated, then it is not too difficult to write the code to estimate that model in NLMIXED.  If the covariance structure
		is complicated (r-side correlations or multiple random effects), then NLMIXED may not be very useful."		
OUTPUT:			
HISTORY:    	03.08.2017 Updating to look at interaction between A and D (a site-specific characteristic/D does not vary within site)
************************************************************************************************************/

options ls=80 ps=55 nodate nonumber nofmterr; *mprint mlogic symbolgen MSGLevel = I;

%let dir=/home/kc171/paper3;

%include "&dir./sim_functions_v01.sas";

%macro setup(var2=, clustnum=); 
	libname in "&dir.";
	
	%if %upcase(&var2)=B %then %do;
		 libname out "&dir./results_datasets/&clustnum.clust";
	%end;
	%else %if %upcase(&var2)=D %then %do;
		 libname out "&dir./results_datasets_D/&clustnum.clust";
	%end;

	data simdata0;
		set in.sim_dat_&clustnum.clust_v02;
	run;

	proc sort data=simdata0; 
		by Sim Cluster ID;
	run;

	proc print data=simdata0 (firstobs=5500 obs=5600);
	run;

	data simdata (keep=sim cluster id outcome outcomevalue A &var2.);
		set simdata0;
		by Sim Cluster ID;
	
		%if %upcase(&var2)=B %then %do;
			array outc (16) o_A1 o_A2 o_A3 o_A4 o_A5 o_B1 o_B2 o_B3 o_B4 o_C1 o_C2 o_C3 o_D1 o_D2 o_D3 o_E1;
		%end;
		%else %if %upcase(&var2)=D %then %do;
			array outc (16) o_A1_2 o_A2_2 o_A3_2 o_A4_2 o_A5_2 o_B1_2 o_B2_2 o_B3_2 o_B4_2 o_C1_2 o_C2_2 o_C3_2 o_D1_2 o_D2_2 o_D3_2 o_E1_2;
		%end;

		do i = 1 to 16;
			outcome = vname(outc(i));
			outcomevalue = outc(i);
			output;
		end;
	run;

	proc sort data=simdata; by Outcome Sim Cluster ID; run;

	proc print data=simdata (firstobs=5500 obs=5600);
	run;

%mend setup;

%setup(var2=B, clustnum=275);
	


%macro run_all(y=,x=);

	/* NAIVE MODELS*/
	
	/* 1. logistic not accounting for clustering */
	*%naive_logistic();

	/* 2. log binomial not accounting for clustering */
	*%naive_logbin();


	/* MARGINAL MODELS -- ACCOUNTING FOR CLUSTERING */

	/* 3. marginal logistic */
	*%marg_logistic();

	/* 4. marginal log binomial */
	*%marg_logbin();

	/* 5. marginal Poisson approximation to log binomial */
	*%marg_pois();
	
	/* 6. COPY method */
	*%marg_COPY(c=1000);

	/* MIXED EFECT MODELS */

	/* 7. logistic with random intercept*/
	*%mixed_logistic();

	/* 8. log binomial with random intercept */
	%mixed_logbin();

	/* 9. Poisson approximation to log binomial (default method is pseudo-likelihood / linearization) with random intercept*/
	*%mixed_pois();

	/* 10. COPY method with random intercept */
	*%mixed_COPY(c=1000);

	/* 9b. Poisson approximation to log binomial (method=Laplace) with random intercept / SAME AS R */
	*%mixed_pois_laplace();

	/* 10b. COPY method with random intercept */
	*%mixed_COPY_laplace(c=1000);

	/* 8b. log binomial with random intercept (method=laplace) */
	*%mixed_logbin_laplace();

	/* 8c. log binomial in NLMIXED (using glm estimates as starting values) */
	*%nlmixed_logbin();

 	/* 8d. log binomial with random intercept and empirical (robust) covariance for fixed effects */
	%mixed_logbin_empirical();

 	/* 9c. poisson with random intercept and empirical (robust) covariance for fixed effects */
	*%mixed_pois_empirical();

	/* 9d. poisson with random intercept (method=laplace) and empirical (robust) covariance for fixed effects */
	*%mixed_pois_laplace_empirical();

	/* 9e. poisson with random intercept and scale */
	*%mixed_pois_scale();

	/* 9f. poisson with random intercept (method=laplace) and scale */
	/*ERROR: R-side random effects are not supported for METHOD=LAPLACE.*/
	*%mixed_pois_laplace_scale();

	/* 9g. poisson with random intercept and scale and empirical (robust) covariance*/
	*%mixed_pois_scale_emp();

	/* 9h. poisson with random intercept (method=laplace) and scale and empirical (robust) covariance */
	/* ERROR: R-side random effects are not supported for METHOD=LAPLACE.*/
	*%mixed_pois_laplace_scale_emp();

%mend run_all;
 
ods exclude all;

%run_all(y=outcomevalue, x=A B A*B);
*%run_all(y=outcomevalue, x=A D A*D);




