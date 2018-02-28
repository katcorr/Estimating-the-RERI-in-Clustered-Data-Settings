/************************************************************************************************************
PROJECT:	simulation study for estimating the RERI in clustered data settings
PROGRAM:	sim_functions_v02.sas
PROGRAMMER:     kat
DATE:		2/7/2017	
RE:		write macros for the various procedures we'll be comparing 
NOTE:		
OUTPUT:			
HISTORY:    
************************************************************************************************************/

/*************************************************/
/*************************************************/

%macro naive_logistic();

	proc logistic data=simdata descending;
		by Outcome Sim;
		model &y. = &x. / covb;
		ods output ConvergenceStatus = out.naive_logistic_convergence;
		ods output ParameterEstimates = out.naive_logistic_est;
		ods output covb = out.naive_logistic_covb;
	run;

%mend naive_logistic;

/*************************************************/
/*************************************************/

%macro naive_logbin();

	proc genmod data=simdata descending;
		by Outcome Sim;
		model &y. = &x. / dist=bin link=log covb;
		ods output ConvergenceStatus = out.naive_logbin_convergence;
		ods output ParameterEstimates = out.naive_logbin_est;
		ods output covb = out.naive_logbin_covb;
	run;

%mend naive_logbin;

/*************************************************/
/*************************************************/

%macro marg_logistic();

	proc genmod data=simdata descending;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=bin link=logit;
		repeated subject=Cluster / type=CS ecovb;
		ods output ConvergenceStatus = out.marg_logistic_convergence;
		ods output GEEEmpPEst = out.marg_logistic_est;
		ods output GEERCov = out.marg_logistic_ecovb;
	run;

%mend marg_logistic;

/*************************************************/
/*************************************************/

%macro marg_logbin();

	proc genmod data=simdata descending;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=bin link=log;
		repeated subject=Cluster / type=CS ecovb;
		ods output ConvergenceStatus = out.marg_logbin_convergence;
		ods output GEEEmpPEst = out.marg_logbin_est;
		ods output GEERCov = out.marg_logbin_ecovb;
	run;

%mend marg_logbin;

/*************************************************/
/*************************************************/

%macro marg_pois();

	proc genmod data=simdata;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log;
		repeated subject=Cluster / type=CS ecovb;
		ods output ConvergenceStatus = out.marg_pois_convergence;
		ods output GEEEmpPEst = out.marg_pois_est;
		ods output GEERCov = out.marg_pois_ecovb;
	run;

%mend marg_pois;


/*************************************************/
/*************************************************/

%macro marg_COPY(c=);

	data simdataCOPY;
		set simdata;
		by Outcome Sim Cluster ID;

		COPYoutcome = &y.;
		COPYweight = (&c.-1)/&c.;
		output;

		COPYoutcome = 1 - &y.;
		COPYweight = 1/&c.;
		output;
	run;

	proc genmod data=simdataCOPY descending;
		by Outcome Sim;
		class Cluster;
		weight COPYweight;
		model COPYoutcome = &x. / dist=bin link=log;
		repeated subject=Cluster / type=CS ecovb;
		ods output ConvergenceStatus = out.marg_COPY_convergence;
		ods output GEEEmpPEst = out.marg_COPY_est;
		ods output GEERCov = out.marg_COPY_ecovb;
	run;

%mend marg_COPY;


/*************************************************/
/*************************************************/

%macro mixed_logistic();

	proc glimmix data=simdata initglm;
		by Outcome Sim;
		class Cluster;
		model &y. (descending) = &x. / dist=binary link=logit s covb;
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_logistic_convergence;
		ods output ParameterEstimates = out.mixed_logistic_est;
		ods output g = out.mixed_logistic_g;
		ods output covb = out.mixed_logistic_covb;
	run;

%mend mixed_logistic;

/*************************************************/
/*************************************************/

%macro mixed_logbin();

	proc glimmix data=simdata initglm;
		by Outcome Sim;
		class Cluster;
		model &y. (descending) = &x. / dist=binary link=log s covb; 
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_logbin_convergence;
		ods output ParameterEstimates = out.mixed_logbin_est;
		ods output g = out.mixed_logbin_g;
		ods output covb = out.mixed_logbin_covb;
	run;
	

%mend mixed_logbin;

/*************************************************/
/*************************************************/

%macro mixed_pois();

	proc glimmix data=simdata initglm;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log s covb;
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_pois_convergence;
		ods output ParameterEstimates = out.mixed_pois_est;
		ods output g = out.mixed_pois_g;
		ods output covb = out.mixed_pois_covb;
	run;


%mend mixed_pois;


/*************************************************/
/*************************************************/

%macro mixed_COPY(c=);

	data simdataCOPY;
		set simdata;
		by Outcome Sim Cluster ID;

		COPYoutcome = &y.;
		COPYweight = (&c.-1)/&c.;
		output;

		COPYoutcome = 1 - &y.;
		COPYweight = 1/&c.;
		output;
	run;

	proc glimmix data=simdataCOPY initglm;
		by Outcome Sim;
		class Cluster;
		weight COPYweight;
		model COPYoutcome (descending) = &x. / dist=binary link=log s covb;
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_copy_convergence;
		ods output ParameterEstimates = out.mixed_copy_est;
		ods output g = out.mixed_copy_g;
		ods output covb = out.mixed_copy_covb;
	run;

%mend mixed_COPY;

/*************************************************/
/*************************************************/

%macro mixed_pois_laplace();

	proc glimmix data=simdata method=laplace;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log s covb;
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_pois_lap_convergence;
		ods output ParameterEstimates = out.mixed_pois_lap_est;
		ods output g = out.mixed_pois_lap_g;
		ods output covb = out.mixed_pois_lap_covb;
	run;


%mend mixed_pois_laplace;


/*************************************************/
/*************************************************/

%macro mixed_COPY_laplace(c=);

	data simdataCOPY;
		set simdata;
		by Outcome Sim Cluster ID;

		COPYoutcome = &y.;
		COPYweight = (&c.-1)/&c.;
		output;

		COPYoutcome = 1 - &y.;
		COPYweight = 1/&c.;
		output;
	run;

	proc glimmix data=simdataCOPY method=laplace;
		by Outcome Sim;
		class Cluster;
		weight COPYweight;
		model COPYoutcome (descending) = &x. / dist=binary link=log s covb;
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_copy_laplace_convergence;
		ods output ParameterEstimates = out.mixed_copy_laplace_est;
		ods output g = out.mixed_copy_laplace_g;
		ods output covb = out.mixed_copy_laplace_covb;
	run;

%mend mixed_COPY_laplace;


/*************************************************/
/*************************************************/

%macro mixed_logbin_laplace();

	proc glimmix data=simdata method=laplace;
		by Outcome Sim;
		class Cluster;
		model &y. (descending) = &x. / dist=binary link=log s covb; 
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_logbin_laplace_convergence;
		ods output ParameterEstimates = out.mixed_logbin_laplace_est;
		ods output g = out.mixed_logbin_laplace_g;
		ods output covb = out.mixed_logbin_laplace_covb;
	run;
	

%mend mixed_logbin_laplace;

/*************************************************/
/*************************************************/

%macro mixed_logbin_empirical();

	proc glimmix data=simdata initglm empirical=classical;
		by Outcome Sim;
		class Cluster;
		model &y. (descending) = &x. / dist=binary link=log s covb; 
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_logbin_emp_convergence;
		ods output ParameterEstimates = out.mixed_logbin_emp_est;
		ods output g = out.mixed_logbin_emp_g;
		ods output covb = out.mixed_logbin_emp_ecovb;
	run;
	

%mend mixed_logbin_empirical;

/*************************************************/
/*************************************************/

%macro mixed_pois_empirical();

	proc glimmix data=simdata initglm empirical=classical;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log s covb;
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_pois_emp_convergence;
		ods output ParameterEstimates = out.mixed_pois_emp_est;
		ods output g = out.mixed_pois_emp_g;
		ods output covb = out.mixed_pois_emp_ecovb;
	run;

%mend mixed_pois_empirical;

/*************************************************/
/*************************************************/

%macro mixed_pois_laplace_empirical();

	proc glimmix data=simdata method=laplace empirical=classical;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log s covb;
		random int / subject=Cluster g;
		ods output ConvergenceStatus = out.mixed_pois_lap_emp_conv;
		ods output ParameterEstimates = out.mixed_pois_lap_emp_est;
		ods output g = out.mixed_pois_lap_emp_g;
		ods output covb = out.mixed_pois_lap_emp_ecovb;
	run;


%mend mixed_pois_laplace_empirical;


/*************************************************/
/*************************************************/

%macro mixed_pois_scale();

	proc glimmix data=simdata initglm;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log s covb;
		random int / subject=Cluster g;
		random _residual_;
		ods output ConvergenceStatus = out.mixed_pois_scale_convergence;
		ods output ParameterEstimates = out.mixed_pois_scale_est;
		ods output g = out.mixed_pois_scale_g;
		ods output covb = out.mixed_pois_scale_covb;
	run;

%mend mixed_pois_scale;

/*************************************************/
/*************************************************/

%macro mixed_pois_laplace_scale();

	proc glimmix data=simdata method=laplace;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log s covb;
		random int / subject=Cluster g;
		random _residual_;
		ods output ConvergenceStatus = out.mixed_pois_lap_scale_conv;
		ods output ParameterEstimates = out.mixed_pois_lap_scale_est;
		ods output g = out.mixed_pois_lap_scale_g;
		ods output covb = out.mixed_pois_lap_scale_covb;
	run;


%mend mixed_pois_laplace_scale;


/*************************************************/
/*************************************************/

%macro mixed_pois_scale_emp();

	proc glimmix data=simdata initglm empirical=classical;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log s covb;
		random int / subject=Cluster g;
		random _residual_;
		ods output ConvergenceStatus = out.mixed_pois_scale_emp_conv;
		ods output ParameterEstimates = out.mixed_pois_scale_emp_est;
		ods output g = out.mixed_pois_scale_emp_g;
		ods output covb = out.mixed_pois_scale_emp_ecovb;
	run;

%mend mixed_pois_scale_emp;

/*************************************************/
/*************************************************/

%macro mixed_pois_laplace_scale_emp();

	proc glimmix data=simdata method=laplace empirical=classical;
		by Outcome Sim;
		class Cluster;
		model &y. = &x. / dist=Poisson link=log s covb;
		random int / subject=Cluster g;
		random _residual_;
		ods output ConvergenceStatus = out.mixed_pois_lap_scale_emp_conv;
		ods output ParameterEstimates = out.mixed_pois_lap_scale_emp_est;
		ods output g = out.mixed_pois_lap_scale_emp_g;
		ods output covb = out.mixed_pois_lap_scale_emp_ecovb;
	run;


%mend mixed_pois_laplace_scale_emp;