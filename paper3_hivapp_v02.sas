/*****************************************************************************************************************************
PROD PROGRAM:  paper3_hivapp_v01.sas																											 										
WORK PROGRAM:  /home/kc171/paper3_application/paper3_hivapp_v01.sas																											 
																															 
PURPOSE:       consider additive interaction between low first cd4 count and NVP use at conception on preterm birth
																												 
SOURCE PRGM:   																										 
INPUT:         																											 
OUTPUT:       
MACROS USED:   																										 
EXEMPTIONS:    
							
AUTHOR:        Kat Correia
CREATION DATE: 7/19/2017

NOTES:         
MODIFICATIONS: 
*****************************************************************************************************************************/

options ls=80 ps=55 nodate nonumber nofmterr mprint mlogic symbolgen MergeNoBy=warn MSGLevel = I;

libname in "/home/kc171/paper3_application";

/**********************/
/* NAIVE LOG BINOMIAL */
/**********************/

%macro naive_logbin(var1,var2,covars);

	data temp;
		set arvs;
		X1=&var1.;
		X2=&var2.;
	run;

	proc genmod data=temp descending;
		model preterm = X1 X2 X1*X2 &covars. / dist=bin link=log covb;
		ods output ParameterEstimates = naive_logbin_est;
		ods output covb = naive_logbin_covb;
	run;
	
	proc print data=naive_logbin_est; run;
	proc print data=naive_logbin_covb; run;

	proc transpose data=naive_logbin_est
			out=reri00;
		id parameter;
		var estimate;
	run;

	proc transpose data=naive_logbin_covb
		out=cov0_var1;
		id rowname;
		var prm2;
	run;

	proc transpose data=naive_logbin_covb
		out=cov0_var2;
		id rowname;
		var prm3;
	run;

	proc transpose data=naive_logbin_covb
		out=cov0_i;
		id rowname;
		var prm4;
	run;

	data reri0;
		length one two $20 covars $100;
		merge reri00
		      cov0_var1 (keep=prm2 prm3 prm4
			 	rename=(prm2=Var_X1 prm3=Cov_X1X2 prm4=Cov_X1I))
	      	      cov0_var2 (keep=prm3 prm4
				rename=(prm3=Var_X2 prm4=COV_X2I))
	              cov0_i (keep=prm4
			      rename=(prm4=Var_I));

		one = "&var1.";
		two = "&var2.";
		covars = "&covars.";

		reri = exp(X1_X2 + X1 + X2) - exp(X1) - exp(X2) + 1;

		/* Delta method CI */
		theta1 = exp(X1 + X2 + X1_X2);
		theta2 = exp(X1);
		theta3 = exp(X2);

		RERI_delta_variance = ((theta1 - theta2)**2)*Var_X1 + ((theta1 - theta3)**2)*Var_X2 + (theta1**2)*Var_I 
				+ 2*( ((theta1-theta2)*(theta1-theta3)*Cov_X1X2) + ((theta1-theta2)*theta1*Cov_X1I) 
						+ ((theta1-theta3)*theta1*Cov_X2I) );

		Lower_delta = RERI - 1.96*sqrt(RERI_delta_variance);
		Upper_delta = RERI + 1.96*sqrt(RERI_delta_variance);

		/* Zou MOVER CI */
		Var_123 = Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I);
	
		l1 = exp(X1 + X2 + X1_X2 - 1.96*sqrt(Var_123));
		u1 = exp(X1 + X2 + X1_X2 + 1.96*sqrt(Var_123));
		l2 = exp(X1 - 1.96*sqrt(Var_X1));
		u2 = exp(X1 + 1.96*sqrt(Var_X1));
		l3 =  exp(X2 - 1.96*sqrt(Var_X2));
		u3 =  exp(X2 + 1.96*sqrt(Var_X2));

		r12 = (Var_X1 + Cov_X1X2 + Cov_X1I) / sqrt(Var_X1*(Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I)));
		r13 = (Cov_X1X2 + Var_X2 + Cov_X2I) / sqrt(Var_X2*(Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I)));
		r23 = Cov_X1X2 / sqrt(Var_X1*Var_X2);


		Lower_Zou = 1 + theta1 - theta2 - theta3 - ( (theta1 - l1)**2 + (u2 - theta2)**2 + (u3 - theta3)**2 - 2*r12*(theta1 - l1)*(u2 - theta2)
				- 2*r13*(theta1-l1)*(u3 - theta3) + 2*r23*(u2 - theta2)*(u3 - theta3) )**0.5;

		Upper_Zou = 1 + theta1 - theta2 - theta3 + ( (u1 - theta1)**2 + (theta2 - l2)**2 + (theta3 - l3)**2 - 2*r12*(u1 - theta1)*(theta2 - l2)
				- 2*r13*(u1 - theta1)*(theta3-l3) + 2*r23*(theta2 - l2)*(theta3 - l3) )**0.5;

	run;

	data naive_reri;
		set naive_reri
		    reri0 (in=new);
		if new then num=_n_;
	run;

	proc datasets memtype=data library=work;
		delete temp naive_logbin_est naive_logbin_covb reri00 cov0_var1 cov0_var2 cov0_i reri0;
	run; quit;

%mend naive_logbin;


/*****************************/
/* MARGINAL/GEE LOG BINOMIAL */
/*****************************/

%macro marg_logbin(var1,var2,covars);

	proc freq data=arvs;
		table &var1.*&var2.*preterm / nocol nopercent;
	run;

	data temp;
		set arvs;
		X1=&var1.;
		X2=&var2.;
	run;

	proc genmod data=temp descending;
		class instn1;
		model preterm = X1 X2 X1*X2 &covars. / dist=bin link=log;
		repeated subject=instn1 / type=CS ecovb;
		ods output GEEEmpPEst = marg_logbin_est;
		ods output GEERCov = marg_logbin_ecovb;
	run;
	
	proc transpose data=marg_logbin_est
			out=reri00;
		id parm;
		var estimate;
	run;

	proc transpose data=marg_logbin_ecovb
		out=cov0_var1;
		id rowname;
		var prm2;
	run;

	proc transpose data=marg_logbin_ecovb
		out=cov0_var2;
		id rowname;
		var prm3;
	run;

	proc transpose data=marg_logbin_ecovb
		out=cov0_i;
		id rowname;
		var prm4;
	run;

	data reri0;
		length one two $20 covars $100;
		merge reri00
		      cov0_var1 (keep=prm2 prm3 prm4
			 	rename=(prm2=Var_X1 prm3=Cov_X1X2 prm4=Cov_X1I))
	      	      cov0_var2 (keep=prm3 prm4
				rename=(prm3=Var_X2 prm4=COV_X2I))
	              cov0_i (keep=prm4
			      rename=(prm4=Var_I));

		one = "&var1.";
		two = "&var2.";
		covars = "&covars.";

		reri = exp(X1_X2 + X1 + X2) - exp(X1) - exp(X2) + 1;

		/* Delta method CI */
		theta1 = exp(X1 + X2 + X1_X2);
		theta2 = exp(X1);
		theta3 = exp(X2);

		RERI_delta_variance = ((theta1 - theta2)**2)*Var_X1 + ((theta1 - theta3)**2)*Var_X2 + (theta1**2)*Var_I 
				+ 2*( ((theta1-theta2)*(theta1-theta3)*Cov_X1X2) + ((theta1-theta2)*theta1*Cov_X1I) 
						+ ((theta1-theta3)*theta1*Cov_X2I) );

		Lower_delta = RERI - 1.96*sqrt(RERI_delta_variance);
		Upper_delta = RERI + 1.96*sqrt(RERI_delta_variance);

		/* Zou MOVER CI */
		Var_123 = Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I);
	
		l1 = exp(X1 + X2 + X1_X2 - 1.96*sqrt(Var_123));
		u1 = exp(X1 + X2 + X1_X2 + 1.96*sqrt(Var_123));
		l2 = exp(X1 - 1.96*sqrt(Var_X1));
		u2 = exp(X1 + 1.96*sqrt(Var_X1));
		l3 =  exp(X2 - 1.96*sqrt(Var_X2));
		u3 =  exp(X2 + 1.96*sqrt(Var_X2));

		r12 = (Var_X1 + Cov_X1X2 + Cov_X1I) / sqrt(Var_X1*(Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I)));
		r13 = (Cov_X1X2 + Var_X2 + Cov_X2I) / sqrt(Var_X2*(Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I)));
		r23 = Cov_X1X2 / sqrt(Var_X1*Var_X2);


		Lower_Zou = 1 + theta1 - theta2 - theta3 - ( (theta1 - l1)**2 + (u2 - theta2)**2 + (u3 - theta3)**2 - 2*r12*(theta1 - l1)*(u2 - theta2)
				- 2*r13*(theta1-l1)*(u3 - theta3) + 2*r23*(u2 - theta2)*(u3 - theta3) )**0.5;

		Upper_Zou = 1 + theta1 - theta2 - theta3 + ( (u1 - theta1)**2 + (theta2 - l2)**2 + (theta3 - l3)**2 - 2*r12*(u1 - theta1)*(theta2 - l2)
				- 2*r13*(u1 - theta1)*(theta3-l3) + 2*r23*(theta2 - l2)*(theta3 - l3) )**0.5;

	run;

	data reri;
		set reri
		    reri0 (in=new);
		if new then num=_n_;
	run;

	proc datasets memtype=data library=work;
		delete temp marg_logbin_est marg_logbin_ecovb reri00 cov0_var1 cov0_var2 cov0_i reri0;
	run; quit;

%mend marg_logbin;

/***************************/
/* RANDOM INTERCEPT MODELS */
/***************************/

%macro randint(var1,var2,covars,dist,meth);

	data temp;
		set arvs;
		X1=&var1.;
		X2=&var2.;
	run;

	proc glimmix data=temp method=&meth.;
		class instn1;
		model preterm (descending) = X1 X2 X1*X2 &covars. / dist=&dist. link=log s covb; 
		random int / subject=instn1 g;
		ods output ParameterEstimates = mixed_logbin_est;
		ods output covb = mixed_logbin_covb;
	run;
	
	proc transpose data=mixed_logbin_est
			out=reri00;
		id effect;
		var estimate;
	run;

	proc transpose data=mixed_logbin_covb
		out=cov0_var1;
		id effect;
		var col2;
	run;

	proc transpose data=mixed_logbin_covb
		out=cov0_var2;
		id effect;
		var col3;
	run;

	proc transpose data=mixed_logbin_covb
		out=cov0_i;
		id effect;
		var col4;
	run;

	data reri0;
		length one two $20 covars $100 model $50;
		merge reri00
		      cov0_var1 (keep=x1 x2 x1_x2
			 	rename=(x1=Var_X1 x2=Cov_X1X2 x1_x2=Cov_X1I))
	      	      cov0_var2 (keep=x2 x1_x2
				rename=(x2=Var_X2 x1_x2=COV_X2I))
	              cov0_i (keep=x1_x2
			      rename=(x1_x2=Var_I));

		one = "&var1.";
		two = "&var2.";
		covars = "&covars.";
		model = "&dist. &meth.";

		reri = exp(X1_X2 + X1 + X2) - exp(X1) - exp(X2) + 1;

		/* Delta method CI */
		theta1 = exp(X1 + X2 + X1_X2);
		theta2 = exp(X1);
		theta3 = exp(X2);

		RERI_delta_variance = ((theta1 - theta2)**2)*Var_X1 + ((theta1 - theta3)**2)*Var_X2 + (theta1**2)*Var_I 
				+ 2*( ((theta1-theta2)*(theta1-theta3)*Cov_X1X2) + ((theta1-theta2)*theta1*Cov_X1I) 
						+ ((theta1-theta3)*theta1*Cov_X2I) );

		Lower_delta = RERI - 1.96*sqrt(RERI_delta_variance);
		Upper_delta = RERI + 1.96*sqrt(RERI_delta_variance);

		/* Zou MOVER CI */
		Var_123 = Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I);
	
		l1 = exp(X1 + X2 + X1_X2 - 1.96*sqrt(Var_123));
		u1 = exp(X1 + X2 + X1_X2 + 1.96*sqrt(Var_123));
		l2 = exp(X1 - 1.96*sqrt(Var_X1));
		u2 = exp(X1 + 1.96*sqrt(Var_X1));
		l3 =  exp(X2 - 1.96*sqrt(Var_X2));
		u3 =  exp(X2 + 1.96*sqrt(Var_X2));

		r12 = (Var_X1 + Cov_X1X2 + Cov_X1I) / sqrt(Var_X1*(Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I)));
		r13 = (Cov_X1X2 + Var_X2 + Cov_X2I) / sqrt(Var_X2*(Var_X1 + Var_X2 + Var_I + 2*(Cov_X1X2 + Cov_X1I + Cov_X2I)));
		r23 = Cov_X1X2 / sqrt(Var_X1*Var_X2);


		Lower_Zou = 1 + theta1 - theta2 - theta3 - ( (theta1 - l1)**2 + (u2 - theta2)**2 + (u3 - theta3)**2 - 2*r12*(theta1 - l1)*(u2 - theta2)
				- 2*r13*(theta1-l1)*(u3 - theta3) + 2*r23*(u2 - theta2)*(u3 - theta3) )**0.5;

		Upper_Zou = 1 + theta1 - theta2 - theta3 + ( (u1 - theta1)**2 + (theta2 - l2)**2 + (theta3 - l3)**2 - 2*r12*(u1 - theta1)*(theta2 - l2)
				- 2*r13*(u1 - theta1)*(theta3-l3) + 2*r23*(theta2 - l2)*(theta3 - l3) )**0.5;

	run;

	data cond_reri;
		set cond_reri
		    reri0 (in=new);
		if new then num=_n_;
	run;

	proc datasets memtype=data library=work;
		delete temp mixed_logbin_est mixed_logbin_covb reri00 cov0_var1 cov0_var2 cov0_i reri0;
	run; quit;

%mend randint;



/***************************************************************************************/
/***************** 	RUN MODELS & REPORT RESULTS 	   *****************************/
/***************************************************************************************/

data arvs;
	set in.ptd_int (where=(preterm ne . and lowfcd4 ne . and onNVP_conception ne . and race_black ne .));

	/* three age groups */
	if (. < momage_conception < 30) or (momage_conception=.) then agegrp3_1 = 1;
	else agegrp3_1 = 0;
	
	if 30 <= momage_conception < 40 then agegrp3_2 = 1;
	else agegrp3_2 = 0;

	if momage_conception >= 40 then agegrp3_3 = 1;
	else agegrp3_3 = 0;
run;

proc sort data=arvs; by instn1; run;

/*** NAIVE MODEL ***/
data naive_reri;
run;

%naive_logbin(var1=onNVP_conception, var2=lowfcd4, covars=);
%naive_logbin(var1=onNVP_conception, var2=lowfcd4, covars=race_black);
%naive_logbin(var1=onNVP_conception, var2=lowfcd4, covars=agegrp3_2 agegrp3_3);
%naive_logbin(var1=onNVP_conception, var2=lowfcd4, covars=agegrp3_2 agegrp3_3 race_black);

/*** MARGINAL GEE ***/
data reri;
run;

%marg_logbin(var1=onNVP_conception, var2=lowfcd4, covars=);
%marg_logbin(var1=onNVP_conception, var2=lowfcd4, covars=race_black);
%marg_logbin(var1=onNVP_conception, var2=lowfcd4, covars=agegrp3_2 agegrp3_3);
%marg_logbin(var1=onNVP_conception, var2=lowfcd4, covars=agegrp3_2 agegrp3_3 race_black);


/*** RANDOM INTERCEPTS ***/
data cond_reri;
run;

%randint(var1=onNVP_conception,var2=lowfcd4,covars=,dist=binary,meth=RSPL initglm);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=,dist=binary,meth=Laplace);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=,dist=Poisson,meth=RSPL);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=,dist=Poisson,meth=Laplace);

%randint(var1=onNVP_conception,var2=lowfcd4,covars=race_black,dist=binary,meth=RSPL initglm);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=race_black,dist=binary,meth=Laplace);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=race_black,dist=binary,meth=Laplace noinitglm);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=race_black,dist=Poisson,meth=RSPL);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=race_black,dist=Poisson,meth=Laplace);

*%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3,dist=binary,meth=RSPL initglm);	*does not converge;
*%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3,dist=binary,meth=Laplace);	*does not converge;
*%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3,dist=binary,meth=Laplace noinitglm); *does not converge;
%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3,dist=Poisson,meth=RSPL);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3,dist=Poisson,meth=Laplace);

*%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3 race_black,dist=binary,meth=RSPL); *does not converge;
*%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3 race_black,dist=binary,meth=Laplace); *does not converge;
*%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3 race_black,dist=binary,meth=Laplace noinitglm); *does not converge;
%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3 race_black,dist=Poisson,meth=RSPL);
%randint(var1=onNVP_conception,var2=lowfcd4,covars=agegrp3_2 agegrp3_3 race_black,dist=Poisson,meth=Laplace);


/*** PRINT RESULTS ***/
ods rtf file="/home/kc171/paper3_application/paper3_hivapp_v02.doc";

	title "Naive model";
	proc print data=naive_reri;
		id covars;
		var reri lower_delta upper_delta lower_zou upper_zou;
	run;

	title "GEE marginal model";
	proc print data=reri;
		id covars;
		var reri lower_delta upper_delta lower_zou upper_zou;
	run;

	title "Random intercepts model";
	proc print data=cond_reri;
		id covars model;
		var reri lower_delta upper_delta lower_zou upper_zou;
	run;

ods rtf close;
