/************************************************************************
File pds_ns_supg_time_integration.c - time stepping procedures

Contains definition of routines:
  pdr_ns_supg_comp_sol_diff_norm - returns max difference between current and old
			  solutions vectors
  pdr_ns_supg_time_integration - time integration driver

------------------------------
History:
	2002    - Krzysztof Banas (pobanas@cyf-kr.edu.pl) for conv-diff
	2011    - Przemyslaw Plaszewski (pplaszew@agh.edu.pl)
	2011    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>

/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>		/* USES */
/* interface for all approximation modules */
#include <modfem/aph_intf.h>		/* USES */
/* interface for all solver modules */
#include <modfem/sih_intf.h>		/* USES */

/* problem dependent module interface */
#include <modfem/pdh_intf.h>		/* USES */
/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>	/* IMPLEMENTS */

#ifdef PARALLEL
#include <modfem/mmph_intf.h>		/* USES */
/* interface for all parallel approximation modules */
#include <modfem/apph_intf.h>		/* USES */
/* interface for parallel communication modules */
#include <modfem/pch_intf.h>		/* USES */
#endif

/* ns_supg problem module's types and functions */
#include <modfem/pd_ns_supg/pdh_ns_supg.h>	/* USES & IMPLEMENTS */
/* functions related to ns_supg weak formulation */
#include <modfem/pd_ns_supg/pdh_ns_supg_weakform.h>
/* types and functions related to ns_supg problem structures */
#include <modfem/pd_ns_supg/pdh_ns_supg_problem.h>
// bc and material header files are included in problem header files

// logs from execution
#include <modfem/uth_log.h>
#include <modfem/uth_io_results.h>

/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */


/*------------------------------------------------------------
pdr_ns_supg_sol_diff_norm - returns max difference in current and old
			  solutions vectors
------------------------------------------------------------*/
void pdr_ns_supg_sol_diff_norm(
  int Current, // in: current vector id
  int Old,     // in: old vector id
  double* sol_diff_norm_ns_supg_p // norm of difference current-old for ns_supg
  )
{

  double sol_dofs_current[APC_MAXELSD];	/* solution dofs */
  double sol_dofs_old[APC_MAXELSD];	/* solution dofs */
  int field_id, mesh_id;
  int node_id = 0;
  double temp, norm_ns_supg = 0.0;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  pdv_ns_supg_current_problem_id = PDC_NS_SUPG_ID;  // ns_supg problem
  i=3; field_id = pdr_ctrl_i_params(pdv_ns_supg_current_problem_id,i);
  mesh_id = apr_get_mesh_id(field_id);

  while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
	if (apr_get_ent_pdeg(field_id, APC_VERTEX, node_id) > 0) {
	  int numdofs = apr_get_ent_nrdofs(field_id, APC_VERTEX, node_id);
	  apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,
			Current, sol_dofs_current);
	  apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,
			Old, sol_dofs_old);
	  i = 0; // for ns_supg we take into account velocities ONLY!!!
	  for (i = 0; i < numdofs-1; ++i) {

/*kbw
  double norm=0.0;
  printf("node_id %d\ti %d\tcurrent %lf\told %lf\tnorm %lf\n",
	 node_id,i,sol_dofs_current[i],sol_dofs_old[i],norm);
  norm += (sol_dofs_current[i] - sol_dofs_old[i]) * (sol_dofs_current[i] - sol_dofs_old[i]);
/*kew*/

	temp = fabs(sol_dofs_current[i] - sol_dofs_old[i]);
	if (norm_ns_supg < temp)
	  norm_ns_supg = temp;
	  }
	}
  }



#ifdef PARALLEL
  pcr_allreduce_max_double(1, &norm_ns_supg, &temp);
  norm_ns_supg = temp;
#endif

  *sol_diff_norm_ns_supg_p = norm_ns_supg;

  return;
}

/*------------------------------------------------------------
pdr_ns_supg_time_integration - time integration driver
------------------------------------------------------------*/
void pdr_ns_supg_time_integration(
  char* Work_dir,
  FILE *Interactive_input,
  FILE *Interactive_output
)
{

  int mesh_id;
  int field_ns_supg_id;
  int  solver_ns_supg_id;
  double sol_norm_uk_ns_supg, sol_norm_un_ns_supg;
  char autodump_filename[300];
  int iadapt = 0;
  char solver_ns_supg_filename[300];
  int nr_iter, nonl_iter;
  double conv_meas, conv_rate;
  int i, iaux;
  double cfl_min=1000000, cfl_max=0, cfl_ave=0.0, daux;

  /*  Navier_Stokes problem parameters */
  pdt_ns_supg_problem *problem_ns_supg = &pdv_ns_supg_problem;
  pdt_ns_supg_ctrls *ctrl_ns_supg = &pdv_ns_supg_problem.ctrl;
  pdt_ns_supg_times *time_ns_supg = &pdv_ns_supg_problem.time;
  pdt_ns_supg_nonls *nonl_ns_supg = &pdv_ns_supg_problem.nonl;
  pdt_ns_supg_linss *lins_ns_supg = &pdv_ns_supg_problem.lins;
  pdt_ns_supg_adpts *adpt_ns_supg = &pdv_ns_supg_problem.adpt;


/*++++++++++++++++ executable statements ++++++++++++++++*/

  pdv_ns_supg_current_problem_id = PDC_NS_SUPG_ID;

  // get mesh and field parameters
  i = 2; mesh_id = pdr_ctrl_i_params(pdv_ns_supg_current_problem_id,i);

  field_ns_supg_id = pdr_ctrl_i_params(pdv_ns_supg_current_problem_id,3);

#ifdef TURBULENTFLOW
  int intStep=0;
  SetMesh(mesh_id);

  SetField(field_ns_supg_id);
  //SetTurbModel(pdv_settings->chlength);
  SetTurbModel(pdv_ns_supg_current_problem_id);

#endif //TURBULENTFLOW

  // when restarting after reaching one of time integration limits
  // new limits are read from interactive input
  if (time_ns_supg->cur_time >= time_ns_supg->final_time) {

	fprintf(Interactive_output,
		"\nCurrent time: %lf is bigger than final time: %lf.\n\n",
		time_ns_supg->cur_time, time_ns_supg->final_time);

	fprintf(Interactive_output,
		"\nCurrent time-step: %d, current time step length: %lf.\n\n",
		time_ns_supg->cur_step, time_ns_supg->cur_dtime);

	if(Interactive_input == stdin){

#ifdef PARALLEL
	  if (pcr_my_proc_id() == pcr_print_master()) {
#endif
	printf("How many steps to perform?\n");
	scanf("%d",&iaux);
#ifdef PARALLEL
	  }
	  pcr_bcast_int(pcr_print_master(), 1, &iaux);
#endif

#ifdef PARALLEL
	  if (pcr_my_proc_id() == pcr_print_master()) {
#endif
	printf("Set new time step length (or max CFL limit if < 0):\n");
	scanf("%lf",&daux);
#ifdef PARALLEL
	  }
	  pcr_bcast_double(pcr_print_master(), 1, &daux);
#endif
	  if(daux>0){
	time_ns_supg->CFL_control = 0;
	time_ns_supg->cur_dtime = daux;
	  }
	  else{
	double daux1, daux2, daux3;
	time_ns_supg->CFL_control = 1;
	time_ns_supg->CFL_limit_max = -daux;
#ifdef PARALLEL
	if (pcr_my_proc_id() == pcr_print_master()) {
#endif
	  printf("Set new time step length decrease multiplier:\n");
	  scanf("%lf",&daux1);
	  printf("Set new min CFL limit (possibly 0 if time step is never increased):\n");
	  scanf("%lf",&daux2);
	  if(daux2>0.0){
		printf("Set new time step length increase multiplier:\n");
		scanf("%lf",&daux3);
	  }
#ifdef PARALLEL
	}
	pcr_bcast_double(pcr_print_master(), 1, &daux1);
	pcr_bcast_double(pcr_print_master(), 1, &daux2);
	pcr_bcast_double(pcr_print_master(), 1, &daux3);
#endif
	time_ns_supg->CFL_time_step_length_decrease_mult = daux1;
	time_ns_supg->CFL_limit_min = daux2;
	time_ns_supg->CFL_time_step_length_increase_mult = daux3;
	  }

	} else{

	  fprintf(Interactive_output, "\nExiting!\n\n");
	  exit(0);;

	}

	time_ns_supg->final_step = time_ns_supg->cur_step + iaux;

	time_ns_supg->final_time = time_ns_supg->cur_time +
	  iaux*time_ns_supg->cur_dtime;
  }

  if (time_ns_supg->cur_step >= time_ns_supg->final_step) {

	fprintf(Interactive_output,
		"\nCurrent time-step: %d is bigger than final step: %d\n\n",
		time_ns_supg->cur_step, time_ns_supg->final_step);

	fprintf(Interactive_output,
		"\nCurrent time: %lf, current time step length: %lf.\n\n",
		time_ns_supg->cur_time, time_ns_supg->cur_dtime);

	if(Interactive_input == stdin){

#ifdef PARALLEL
	  if (pcr_my_proc_id() == pcr_print_master()) {
#endif
	printf("How many steps to perform?\n");
	scanf("%d",&iaux);
#ifdef PARALLEL
	  }
	  pcr_bcast_int(pcr_print_master(), 1, &iaux);
#endif

#ifdef PARALLEL
	  if (pcr_my_proc_id() == pcr_print_master()) {
#endif
	printf("Set new time step length (or CFL limit if < 0):\n");
	scanf("%lf",&daux);
#ifdef PARALLEL
	  }
	  pcr_bcast_double(pcr_print_master(), 1, &daux);
#endif
	  if(daux>0){
	time_ns_supg->CFL_control = 0;
	time_ns_supg->cur_dtime = daux;
	  }
	  else{

	double daux1, daux2, daux3;
	time_ns_supg->CFL_control = 1;
	time_ns_supg->CFL_limit_max = -daux;
#ifdef PARALLEL
	if (pcr_my_proc_id() == pcr_print_master()) {
#endif
	  printf("Set new time step length decrease multiplier:\n");
	  scanf("%lf",&daux1);
	  printf("Set new min CFL limit (possibly 0 if time step is never increased):\n");
	  scanf("%lf",&daux2);
	  if(daux2>0.0){
		printf("Set new time step length increase multiplier:\n");
		scanf("%lf",&daux3);
	  }
#ifdef PARALLEL
	}
	pcr_bcast_double(pcr_print_master(), 1, &daux1);
	pcr_bcast_double(pcr_print_master(), 1, &daux2);
	pcr_bcast_double(pcr_print_master(), 1, &daux3);
#endif
	time_ns_supg->CFL_time_step_length_decrease_mult = daux1;
	time_ns_supg->CFL_limit_min = daux2;
	time_ns_supg->CFL_time_step_length_increase_mult = daux3;
	  }


	} else{

	  fprintf(Interactive_output, "\nExiting!\n\n");
	  exit(0);;

	}

	time_ns_supg->final_step = time_ns_supg->cur_step + iaux;

	time_ns_supg->final_time = time_ns_supg->cur_time +
	  iaux*time_ns_supg->cur_dtime;

  }

  /* print some info */
  fprintf(Interactive_output, "\n\n**************** BEGINNING TIME INTEGRATION ***************\n");
  fprintf(Interactive_output, "\nTime integration will stop whichever comes first:\n");
  fprintf(Interactive_output, "final_time: %lf\n", time_ns_supg->final_time);
  fprintf(Interactive_output, "final_timestep: %d\n", time_ns_supg->final_step);
  fprintf(Interactive_output, "error less than time_integration_tolerance: %lf\n\n", time_ns_supg->conv_meas);

#ifndef PARALLEL
  if(Interactive_input == stdin) {
	printf("Type [Ctrl-C] to manually break time integration.\n");
  }
#endif


  if(pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0){ // mkb interface for solvers

	int problem_id = pdv_ns_supg_current_problem_id;
	int solver_type = pdr_lins_i_params(PDC_NS_SUPG_ID, 1);
	int max_iter = -1;
	int error_type = -1;
	double error_tolerance = -1;
	int monitoring_level = -1;

	// when no parameter file passed - take control parameters from problem input file
	if(0==strlen(ctrl_ns_supg->solver_filename)){

	  strcpy(solver_ns_supg_filename, ctrl_ns_supg->solver_filename);
	  max_iter = pdr_lins_i_params(problem_id, 2); // max_iter
	  error_type = pdr_lins_i_params(problem_id, 3); // error_type
	  error_tolerance = pdr_lins_d_params(problem_id, 4); // error_tolerance
	  monitoring_level = pdr_lins_i_params(problem_id, 5);  // monitoring level

	}
	else{

	  sprintf(solver_ns_supg_filename, "%s/%s",
		  ctrl_ns_supg->work_dir, ctrl_ns_supg->solver_filename);

	}

	int parallel = SIC_SEQUENTIAL;
#ifdef PARALLEL
	parallel = SIC_PARALLEL;
#endif

/*kbw*/
	fprintf(Interactive_output, "initializing linear solver (file %s)\n", solver_ns_supg_filename);
	fprintf(Interactive_output, "parameters: parallel %d, maxiter %d, error_type %d, error_meas %.15lf, monitor %d\n",
		parallel, max_iter,  error_type,  error_tolerance, monitoring_level);
/*kew*/
	int num_levels = pdr_get_max_num_grid_levels(pdv_ns_supg_current_problem_id);
	solver_ns_supg_id = sir_init(solver_type, parallel, num_levels,
				 solver_ns_supg_filename, max_iter,
				 error_type, error_tolerance, monitoring_level);

	ctrl_ns_supg->solver_id = solver_ns_supg_id;
	fprintf(Interactive_output, "Assigned solver ID: %d for NS_SUPG\n",
		solver_ns_supg_id);
  }

  /* set parameter indicating new mesh */
  iadapt = 1;

  /* start loop over time steps */
  while (utv_SIGINT_not_caught) {

	/* update time step and current time */
	++(time_ns_supg->cur_step);

	time_ns_supg->prev_dtime = time_ns_supg->cur_dtime;


	// here possible time step length adaptations !!!
	// compute maximal, minimal and average CFL number for ns_supg_problem
	if(time_ns_supg->CFL_control == 1){

	  /* time step length adaptation based on specified CFL number */
	  int cfl_control = time_ns_supg->CFL_control;
	  double cfl_limit_max = time_ns_supg->CFL_limit_max;
	  double cfl_limit_min = time_ns_supg->CFL_limit_min;
	  double cfl_increase_mult = time_ns_supg->CFL_time_step_length_increase_mult;
	  double cfl_decrease_mult = time_ns_supg->CFL_time_step_length_decrease_mult;


	  fprintf(Interactive_output,
		  "\nCFL numbers before time step length control:\n");
	  pdr_ns_supg_compute_CFL(PDC_NS_SUPG_ID, &cfl_min, &cfl_max, &cfl_ave);

#ifdef PARALLEL
	  {
	double cfl_max_new = cfl_max;
	pcr_allreduce_max_double(1,&cfl_max,&cfl_max_new);
	cfl_max = cfl_max_new;
	  }
#endif

	  fprintf(Interactive_output,
		  "CFL_min = %lf, CFL_max (global) = %lf, CFL_average = %lf\n",
		  cfl_min,cfl_max,cfl_ave);

	  // check whether CFL is not too big
	  if(cfl_max > cfl_limit_max){
	time_ns_supg->cur_dtime = time_ns_supg->cur_dtime*cfl_decrease_mult;
	  }
	  // check whether CFL is not too small
	  if(cfl_max < cfl_limit_min ){
	time_ns_supg->cur_dtime=time_ns_supg->cur_dtime*cfl_increase_mult;
	  }
	}


	fprintf(Interactive_output, "\n\n***************************************************************************\n");
	fprintf(Interactive_output, "\n\n***************************************************************************\n");
	fprintf(Interactive_output, "Solving time step %d (step_length: %lf, initial time: %lf)\n",
		time_ns_supg->cur_step, time_ns_supg->cur_dtime, time_ns_supg->cur_time);

	time_ns_supg->cur_time += time_ns_supg->cur_dtime;

	utr_io_result_add_value_int(RESULT_STEP, time_ns_supg->cur_step);
	utr_io_result_add_value_double(RESULT_CUR_TIME, time_ns_supg->cur_time);
	utr_io_result_add_value_double(RESULT_DT, time_ns_supg->cur_dtime);

	// compute maximal, minimal and average CFL number for ns_supg_problem
	pdr_ns_supg_compute_CFL(PDC_NS_SUPG_ID, &cfl_min, &cfl_max, &cfl_ave);
#ifdef PARALLEL
	{
	  double cfl_max_new = cfl_max;
	  pcr_allreduce_max_double(1,&cfl_max,&cfl_max_new);
	  cfl_max = cfl_max_new;
	}
#endif

	utr_io_result_add_value_double(RESULT_CFL_MIN, cfl_min);
	utr_io_result_add_value_double(RESULT_CFL_MAX, cfl_max);
	utr_io_result_add_value_double(RESULT_CFL_AVG, cfl_ave);

	fprintf(Interactive_output,
		"CFL_min = %lf, CFL_max = %lf, CFL_average = %lf\n",
		cfl_min,cfl_max,cfl_ave);
	// compute maximal, minimal and average CFL number

	// rewrite solution (COULD BE DONE AFTER FIRST NONLINEAR ITERATION!!!):
	// current from previous time step becomes previous for current time step
	// soldofs_1 are the most recent solution dofs
	// at the beginning of time step they are rewritten to soldofs_3
	// that holds the solution from the previous time step
	utr_rewr_sol(field_ns_supg_id,
			 Current_solution_ID, // 1,
			 Previous_time_step_sol_ID // 3
			 );	//rewrite current -> u_n (3)

	/* update time dependent boundary conditions for NS problem*/
	pdr_ns_supg_update_timedep_bc(&problem_ns_supg->bc,
				  time_ns_supg->cur_time);

	nonl_iter = 0; double initial_sol_norm_uk_ns_supg = 0.0;
	for(;;) {

		utr_io_result_add_value_int(RESULT_NON_LIN_STEP,nonl_iter+1);

	  fprintf(Interactive_output, "\n---------------------------------------------------------------------------\n");
	  fprintf(Interactive_output, "Solving nonlinear iteration %d\n", nonl_iter);

	  // when forming linear system - soldofs_1 are equal to soldofs_2
	  // after solving the system soldofs_1 are different than soldofs_2
	  // before starting new solution soldofs_1 are rewritten to soldofs_2
	  // simply: soldofs_1 are u_k+1, soldofs_2 are u_k and soldofs_3 are u_n
	  utr_rewr_sol(field_ns_supg_id,
		   Current_solution_ID, // 1,
		   Previous_iteration_sol_ID // 2
		   );	//rewrite current -> u_k (2)

	  if (iadapt == 1) {
#ifdef PARALLEL
	/* initiate exchange tables for DOFs - for one field, max_num_levels */
	// offset = 0, nreq = 4;
	int problem_id = pdv_ns_supg_current_problem_id;
	int nr_levels = pdr_get_max_num_grid_levels(problem_id);
	int field_id = pdr_ctrl_i_params(pdv_ns_supg_current_problem_id, 3);

	pdv_exchange_table_index = appr_init_exchange_tables(pcr_nr_proc(),
								 pcr_my_proc_id(),
								 field_id, 0, 4,
								 nr_levels);
#endif
	if(pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0){ // mkb interface for solvers
	  sir_create(solver_ns_supg_id, PDC_NS_SUPG_ID);
	}
	iadapt = 0;
	  }

	  if(pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0){ // mkb interface for solvers

	int problem_id = pdv_ns_supg_current_problem_id;
	int nr_iter = -1;
	double conv_meas = -1.0;
	double conv_rate = -1.0;
	int monitor = -1;

	// when no parameter file passed - take control parameters from problem input file
	if(0==strlen(ctrl_ns_supg->solver_filename)){
	  nr_iter = pdr_lins_i_params(problem_id, 2); // max_iter
	  conv_meas = pdr_lins_d_params(problem_id, 4); // error_tolerance
	  monitor = pdr_lins_i_params(problem_id, 5);  // monitoring level
	}

	int ini_guess = 1; // get initial guess from data structure

	/*---------- CALLING DIRECT OR ITERATIVE SOLVER THROUGH MKB INTERFACE ---------*/
	sir_solve(solver_ns_supg_id, SIC_SOLVE, ini_guess, monitor,
		  &nr_iter, &conv_meas, &conv_rate);
	if(pdr_lins_i_params(problem_id, 1) > 0){ // iterative solver
	  fprintf(Interactive_output,
		  "\nAfter %d iterations of linear solver for ns_supg problem\n",
		  nr_iter);
	  fprintf(Interactive_output,
		  "Convergence measure: %lf, convergence rate %lf\n",
		  conv_meas, conv_rate);
	}


	  } else {
#ifdef PARALLEL
	if(pcv_nr_proc > 1) {
	  mf_log_err("Shared memory PARDISO called for solving in distributed memory computations!");
	}
#endif
	sir_direct_solve_lin_sys(PDC_NS_SUPG_ID, SIC_SEQUENTIAL,
				 solver_ns_supg_filename);
	  }

	  if(nonl_iter==0){
	// rewrite solution - POSSIBLE ALTERNATIVE
	// zero iteration for current time step becomes the previous time step solution
	// i.e. we exchange U^(n-1) with U^(n) ( assuming we solve for U^(n+1)_(k+1) )
	// i.e. in the first iteration we have:
	//              previous time step ->  U^(n-1), previous iteration -> U^(n)
	// i.e. in any other iteration we have:
	//              previous time step ->  U^(n), previous iteration -> U^(n+1)_k
	// i.e. before solver we always have: current = previous =  U^(n+1)_k,
	// i.e. after solver we always have: current = U^(n+1)_(k+1)
	/* utr_rewr_sol(field_ns_supg_id,  */
	/* 	     Previous_iteration_sol_ID, // 2,  */
	/* 	     Previous_time_step_sol_ID // 3 */
	/* 	     );	//rewrite current -> u_n (3)	 */
	  }

	  /* post-process solution using slope limiter */
	  //if(slope) {
	  //  iaux=pdr_slope_limit(problem_id);
	  // pdr_slope_limit for DG is archived in pdd_conv_diff/approx_dg/..._util.c
	  //  if(iaux<0) { printf("\nError in slope!\n");getchar();}
	  //}

	  pdr_ns_supg_sol_diff_norm(
				Current_solution_ID, // 1,
				Previous_iteration_sol_ID, // 2,
				&sol_norm_uk_ns_supg);
	  fprintf(Interactive_output, "\nAfter linear solver in nonlinear iteration %d\n",
		  nonl_iter);
	  fprintf(Interactive_output, "Solution difference norm (u_k, prev. u_k): %lf (limit %lf)\n",
		  sol_norm_uk_ns_supg, nonl_ns_supg->conv_meas);

	  utr_io_result_add_value_double(RESULT_NORM_NS_SUPG_NONL,sol_norm_uk_ns_supg);

	  if(nonl_iter==0) initial_sol_norm_uk_ns_supg = sol_norm_uk_ns_supg;
	  ++nonl_iter;

	  int break_indi=0;
	  if( nonl_ns_supg->conv_type==0 ){

	if( sol_norm_uk_ns_supg < nonl_ns_supg->conv_meas ){

	  fprintf(Interactive_output,
		  "\nConvergence in non-linear iterations: %lf < %lf - breaking.\n",
		  sol_norm_uk_ns_supg, nonl_ns_supg->conv_meas);

	  //break;
	  break_indi = 1;
	}

	  }
	  else if( nonl_ns_supg->conv_type==1 ){

	if( sol_norm_uk_ns_supg < nonl_ns_supg->conv_meas*initial_sol_norm_uk_ns_supg ){
	  fprintf(Interactive_output,
		  "\nConvergence in non-linear iterations: %lf < %lf*%lf - breaking.\n",
		  sol_norm_uk_ns_supg, nonl_ns_supg->conv_meas, initial_sol_norm_uk_ns_supg);

	  //break;
	  break_indi = 1;
	}
	  }


	  if (nonl_iter >= nonl_ns_supg->max_iter) {
	fprintf(Interactive_output,
		"\nMax nonlinear iterations (%d) reached - breaking.\n"
		, nonl_ns_supg->max_iter);


#ifdef TURBULENTFLOW
	fprintf(Interactive_output, "test_tstart");

	if (intStep>0)
	  {
		//RevStep();
	  }
	intStep=intStep+1;
	int jturb;
	for (jturb=0;jturb<1;jturb++)
	  {
		Next(time_ns_supg->cur_time);
	  }
	fprintf(Interactive_output, "test_tstop");


	//		ConfStep();

#endif //TURBULENTFLOW

	//break;
	break_indi = 1;
	  }

	  // check whether we should change time step length
	  if(break_indi == 1){

	if(time_ns_supg->time_step_length_nonl_control == 1){

	  fprintf(Interactive_output,
		  "In time step control (iter_min %d <=? nonl_iter %d <=? iter_max %d)\n",
		  time_ns_supg->time_step_length_nonl_iter_min, nonl_iter,
		  time_ns_supg->time_step_length_nonl_iter_max );

	  /* time step length adaptation based on specified CFL number */
	  int time_step_length_control = time_ns_supg->time_step_length_nonl_control;
	  double time_step_length_limit_max = time_ns_supg->time_step_length_nonl_iter_max;
	  double time_step_length_limit_min = time_ns_supg->time_step_length_nonl_iter_min;
	  double time_step_length_increase_mult =
		time_ns_supg->time_step_length_nonl_iter_increase_mult;
	  double time_step_length_decrease_mult =
		time_ns_supg->time_step_length_nonl_iter_decrease_mult;

	  // if convergence was slow - we decrease time step length
	  if( (nonl_iter>=time_step_length_limit_max)
	  // or if norm is above required value
	  ||  (sol_norm_un_ns_supg > time_ns_supg->time_step_length_max_norm) ) {
		time_ns_supg->cur_dtime = time_ns_supg->cur_dtime*time_step_length_decrease_mult;
	  }
	  // if convergence is fast we increase time step length
	  else if( (nonl_iter<=time_step_length_limit_min)
	  // or if norm is below required value
		   || (sol_norm_un_ns_supg < time_ns_supg->time_step_length_min_norm) ){
		time_ns_supg->cur_dtime = time_ns_supg->cur_dtime*time_step_length_increase_mult;
	  }


	  if(time_ns_supg->cur_dtime < time_ns_supg->time_step_length_min) {
		time_ns_supg->cur_dtime = time_ns_supg->time_step_length_min;
		mf_log_info("Current time step(%lf) is below minimum, resetting to time_step_length_min(%lf)",
					time_ns_supg->cur_dtime, time_ns_supg->time_step_length_min);
	  }
	  else if(time_ns_supg->cur_dtime > time_ns_supg->time_step_length_max) {
		  time_ns_supg->cur_dtime = time_ns_supg->time_step_length_max;
		  mf_log_info("Current time step(%lf) exceeds maximum, resetting to time_step_length_max(%lf)",
					  time_ns_supg->cur_dtime, time_ns_supg->time_step_length_max);
	  }
	} //end of time step adapt.

	break;
	  }


	}  // the end of infinite loop over non-linear iterations


	pdr_ns_supg_sol_diff_norm(
				  Current_solution_ID, // 1,
				  Previous_time_step_sol_ID, // 3,
				  &sol_norm_un_ns_supg);

	fprintf(Interactive_output, "\n---------------------------------------------------------------------------\n");
	fprintf(Interactive_output, "After non-linear solver in time step %d\n",
		time_ns_supg->cur_step);
	fprintf(Interactive_output, "Solution difference norm (u_n, prev. u_n): %lf (limit %lf)\n",
		sol_norm_un_ns_supg, time_ns_supg->conv_meas);

	utr_io_result_add_value_double(RESULT_NORM_NS_SUPG,sol_norm_un_ns_supg);

	/* graphics data dumps */
	if (time_ns_supg->intv_graph > 0 &&
	time_ns_supg->cur_step % time_ns_supg->intv_graph == 0) {

//      sprintf(autodump_filename, "%s/%s_g%d.dmp", ctrl_ns_supg->work_dir, ctrl_ns_supg->field_dmp_filepattern, time_ns_supg->cur_step);
//      if (utr_export_field(field_id, PDC_NREQ, Current_solution_ID, time_ns_supg->intv_graph_accu, autodump_filename) < 0)
//        fprintf(Interactive_output, "Error in writing field data!\n");

#ifdef TURBULENTFLOW
	  DumpParaviewC(time_ns_supg->cur_step);
#endif

	  pdr_ns_supg_write_paraview(Work_dir,
				 Interactive_input, Interactive_output);

	}

	/* full data dumps (for restarting) */
	if (time_ns_supg->intv_dumpout > 0
	&& time_ns_supg->cur_step % time_ns_supg->intv_dumpout == 0) {

	  pdr_ns_supg_dump_data(Work_dir,
							Interactive_input, Interactive_output);

	}

	utr_ctrl_pts_add_values(problem_ns_supg->ctrl.field_id);
	utr_io_result_write_values_and_proceed();

	/* check for stop conditions */



	// stop if computational error occured
	if( pdv_ns_supg_problem.ctrl.error_indicator != NO_ERROR ) {
		mf_log_err("Breaking time integration because of error(%d).",
				   pdv_ns_supg_problem.ctrl.error_indicator);
		break;
	}

	/* stop if convergence reached (for stationary problems) */
	if (sol_norm_un_ns_supg < time_ns_supg->conv_meas) {
	  fprintf(Interactive_output,
		  "\nConvergence in time integration!\n");
	  break;
	}

	/* stop when final time reached */
	if (time_ns_supg->cur_time >= time_ns_supg->final_time) {
	  fprintf(Interactive_output,
		  "\nFinal time reached (%lf >= %lf). Stopping.\n",
		  time_ns_supg->cur_time, time_ns_supg->final_time);
	  //mf_log_info("Final time reached. Stopping.");
	  break;
	}

	/* cur_dtime could change through simulation (time step adaptation) */
	/* thus we need to check also if cur_step not bigger than final_step */
	if (time_ns_supg->cur_step >= time_ns_supg->final_step) {
	  fprintf(Interactive_output,
		  "\nFinal step reached. Stopping.\n");
	  //mf_log_info("Final step reached. Stopping.");
	  break;
	}


	// when time for adaptation
	if( adpt_ns_supg->type>0 && adpt_ns_supg->interval>0 &&
	(time_ns_supg->cur_step+1)%adpt_ns_supg->interval==0 ) {

#ifdef PARALLEL
	  /* free exchange tables for DOFs */
	  appr_free_exchange_tables(pdv_exchange_table_index);
	  pdv_exchange_table_index = -1;
#endif

	  if(pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0){ // mkb interface for solvers
	/* free solver data structures */
	sir_free(solver_ns_supg_id);
	  }
	  pdr_ns_supg_adapt(Work_dir,
			Interactive_input, Interactive_output);
	  /* indicate to recreate block structure */
	  iadapt=1;
	}


	//check if user wants to break time_integration
	/* if(Interactive_input == stdin) { */
	/*   int sec=0; */
	/*   char c=getc(stdin); */
	/*   fflush(stdout); */
	/*   if(c == 'q') { */
	/* 	printf("\nPress [q] again to finalize current step and exit to menu: "); */
	/* 	c=0; */
	/* 	c=getchar(); */
	/* 	if(c == 'q') { */
	/* 	  printf("\nBreaking time integration (user input)!"); */
	/* 	  break; */
	/* 	} */
	/*   } */
	/* } */


  }  //end loop over timesteps


  if (iadapt == 0) {
#ifdef PARALLEL
	/* free exchange tables for DOFs */
	appr_free_exchange_tables(pdv_exchange_table_index);
	pdv_exchange_table_index = -1;
#endif


	if(pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0){ // mkb interface for solvers
	  sir_free(solver_ns_supg_id);
	}

  }

  if(pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0){ // mkb interface for solvers
	sir_destroy(solver_ns_supg_id);
  }

  fprintf(Interactive_output, "\n**************** LEAVING TIME INTEGRATION ***************\n");

  return;
}

