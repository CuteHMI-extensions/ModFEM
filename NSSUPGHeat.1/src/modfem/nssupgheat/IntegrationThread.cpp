#include "../../../include/modfem/nssupgheat/IntegrationThread.hpp"

#include <cstdlib>
#include <iostream>

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>

#include <modfem/pd_ns_supg/pdh_ns_supg_weakform.h>

/* utilities - including simple time measurement library */
#include <modfem/uth_intf.h>		/* USES */
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
#include <modfem/uth_log.h>
#include <modfem/uth_bc.h>
#include <modfem/uth_mesh.h>
#include <modfem/uth_io_results.h>
#include <modfem/uth_system.h>
#include <modfem/uth_io_compression.h>
#include <modfem/uth_system.h>

#ifdef PARALLEL
/* interface of parallel mesh manipulation modules */
#include <modfem/mmph_intf.h>		/* USES */
/* interface for all parallel approximation modules */
#include <modfem/apph_intf.h>		/* USES */
/* interface for parallel communication modules */
#include <modfem/pch_intf.h>		/* USES */
#endif

/* visualization module */
//#include <modfem/mod_fem_viewer.h>	/* USES */

/* problem module's types and functions */
#include <modfem/pd_ns_supg_heat/pdh_ns_supg_heat.h>		/* USES */
/* types and functions related to problem structures */
#include <modfem/pd_ns_supg/pdh_ns_supg_problem.h>
#include <modfem/pd_heat/pdh_heat_problem.h>
// bc and material header files are included in problem header files
#include <modfem/pd_heat/pdh_heat_bc.h>

#include <QDir>
#include <QElapsedTimer>

namespace modfem {
namespace nssupgheat {

IntegrationThread::IntegrationThread(FILE * interactiveInput, FILE * interactiveOutput, const QString & workingDirectory, QObject * parent):
	QThread(parent),
	m(new Members{
	interactiveInput,
	interactiveOutput,
	workingDirectory})
{
}

void IntegrationThread::run()
{
	integrate();
}

void IntegrationThread::integrate()
{
	FILE * interactiveInput = m->interactiveInput;
	FILE * interactiveOutput = m->interactiveOutput;

	QByteArray workDirectory(QDir(m->workingDirectory).path().toLocal8Bit());

	utv_SIGINT_not_caught = 1;

	utr_io_result_set_filename(workDirectory.data(), "ns_supg_heat.csv");
	utr_io_result_clear_columns();
	utr_io_result_add_column(RESULT_STEP);
	utr_io_result_add_column(RESULT_CUR_TIME);
	utr_io_result_add_column(RESULT_DT);
	utr_io_result_add_column(RESULT_CFL_MIN);
	utr_io_result_add_column(RESULT_CFL_MAX);
	utr_io_result_add_column(RESULT_CFL_AVG);
	utr_io_result_add_column(RESULT_NON_LIN_STEP);
	utr_io_result_add_column(RESULT_NORM_NS_SUPG);
	utr_io_result_add_column(RESULT_NORM_NS_SUPG_NONL);
	utr_io_result_add_column(RESULT_NORM_HEAT);
	utr_io_result_add_column(RESULT_NORM_HEAT_NONL);
	utr_io_result_add_column(RESULT_N_DOFS);
	utr_io_result_add_column(RESULT_OUT_FILENAME);
	utr_ctrl_pts_init(workDirectory.data(), "control_points.dat",
			apr_get_nreq(pdv_ns_supg_problem.ctrl.field_id)
			+ apr_get_nreq(pdv_heat_problem.ctrl.field_id),
			"control_points.csv");
	utr_io_result_write_columns(pdv_ns_supg_problem.time.cur_step);
	char mm_name[255] = {0};
	mmr_module_introduce(mm_name);
	//printf("mesh module %s\n", mm_name);
	if ( 0 != strcmp("3D_remesh", mm_name) ) {
		//printf("mesh module different than 3D_remesh\n");
		utr_ctrl_pts_add_values(pdv_ns_supg_problem.ctrl.field_id);
		utr_ctrl_pts_add_values(pdv_heat_problem.ctrl.field_id);
	}
	else {
		mf_log_err("3D_remesh mesh module: adding control points switched off!!!!!\n");
	}
	utr_io_result_write_values_and_proceed();

	mf_log_info("Beginning solution of coupled ns_supg-heat problem");

	double timer_all = time_clock();

	/*---------- main time integration procedure ------------*/

	nsSupgHeatIntegrate(workDirectory.data(), interactiveInput, interactiveOutput);

	timer_all = time_clock() - timer_all;

	fprintf(m->interactiveOutput, "\nExecution time total: %lf\n", timer_all);

}

void IntegrationThread::nsSupgHeatIntegrate(char * Work_dir, FILE * Interactive_input, FILE * Interactive_output)
{

	/* one instance of mesh, three instances of fields, two instances of solvers */
	int mesh_id;
	int field_ns_supg_id, field_heat_id;
	//int field_heat_dtdt_id;
	int  solver_ns_supg_id, solver_heat_id;
	double sol_norm_uk_ns_supg, sol_norm_un_ns_supg;
	double sol_norm_uk_heat, sol_norm_un_heat;
	char autodump_filename[300];
	int iadapt = 0;
	char solver_ns_supg_filename[300], solver_heat_filename[300];;
	int nr_iter, nonl_iter, monitor;
	double conv_meas, conv_rate;
	int i, iaux;
	double cfl_min = 1000000, cfl_max = 0, cfl_ave = 0.0;

	// both problems use the same time integration parameters
	/*  Navier_Stokes problem parameters */
	pdt_ns_supg_problem * problem_ns_supg = &pdv_ns_supg_problem;
	pdt_ns_supg_ctrls * ctrl_ns_supg = &pdv_ns_supg_problem.ctrl;
	pdt_ns_supg_times * time_ns_supg = &pdv_ns_supg_problem.time;
	pdt_ns_supg_nonls * nonl_ns_supg = &pdv_ns_supg_problem.nonl;
	pdt_ns_supg_linss * lins_ns_supg = &pdv_ns_supg_problem.lins;
	pdt_ns_supg_adpts * adpt_ns_supg = &pdv_ns_supg_problem.adpt;

	/*  heat problem parameters */
	pdt_heat_problem * problem_heat = &pdv_heat_problem;
	pdt_heat_ctrls * ctrl_heat = &pdv_heat_problem.ctrl;
	pdt_heat_times * time_heat = &pdv_heat_problem.time;
	pdt_heat_nonls * nonl_heat = &pdv_heat_problem.nonl;
	pdt_heat_linss * lins_heat = &pdv_heat_problem.lins;
	pdt_heat_adpts * adpt_heat = &pdv_heat_problem.adpt;


	/* time adaptation */
	double initial_dtime = time_ns_supg->prev_dtime;
	double final_dtime = time_ns_supg->cur_dtime;
	double dtime_adapt, time_adapt;
	double daux;


	double initial_nonl_error_ns_supg = 0.0;
	double initial_nonl_error_heat = 0.0;
	double previous_un_error_ns_supg = 1.0e6;
	double previous_un_error_heat = 1.0e6;

	/*++++++++++++++++ executable statements ++++++++++++++++*/
	// get mesh and field parameters
	/* get mesh id for ns_supg problem (the same is for heat problem) */
	pdv_ns_supg_heat_current_problem_id = PDC_NS_SUPG_ID;
	i = 2;
	mesh_id = pdr_ctrl_i_params(pdv_ns_supg_heat_current_problem_id, i);

	pdv_ns_supg_heat_current_problem_id = PDC_NS_SUPG_ID;  // ns_supg problem
	field_ns_supg_id = pdr_ctrl_i_params(pdv_ns_supg_heat_current_problem_id, 3);
	pdv_ns_supg_heat_current_problem_id = PDC_HEAT_ID;  // heat problem
	field_heat_id = pdr_ctrl_i_params(pdv_ns_supg_heat_current_problem_id, 3);

	pdv_heat_problem.time.type = pdv_ns_supg_problem.time.type;
	pdv_heat_problem.time.alpha = pdv_ns_supg_problem.time.alpha;
	pdv_heat_problem.time.cur_step = pdv_ns_supg_problem.time.cur_step;
	pdv_heat_problem.time.cur_time = pdv_ns_supg_problem.time.cur_time;
	pdv_heat_problem.time.cur_dtime = pdv_ns_supg_problem.time.cur_dtime;
	pdv_heat_problem.time.prev_dtime = pdv_ns_supg_problem.time.prev_dtime;
	pdv_heat_problem.time.final_time = pdv_ns_supg_problem.time.final_time;
	pdv_heat_problem.time.final_step = pdv_ns_supg_problem.time.final_step;
	pdv_heat_problem.time.monitor = pdv_ns_supg_problem.time.monitor;
	pdv_heat_problem.time.intv_graph = pdv_ns_supg_problem.time.intv_graph;
	pdv_heat_problem.time.intv_dumpout = pdv_ns_supg_problem.time.intv_dumpout;
	pdv_heat_problem.time.time_step_length_nonl_control = pdv_ns_supg_problem.time.time_step_length_nonl_control;
	pdv_heat_problem.time.time_step_length_nonl_iter_max = pdv_ns_supg_problem.time.time_step_length_nonl_iter_max;
	pdv_heat_problem.time.time_step_length_nonl_iter_min = pdv_ns_supg_problem.time.time_step_length_nonl_iter_min;
	pdv_heat_problem.time.time_step_length_nonl_iter_increase_mult = pdv_ns_supg_problem.time.time_step_length_nonl_iter_increase_mult;
	pdv_heat_problem.time.time_step_length_nonl_iter_decrease_mult = pdv_ns_supg_problem.time.time_step_length_nonl_iter_decrease_mult;
	pdv_heat_problem.time.time_step_length_min = pdv_ns_supg_problem.time.time_step_length_min;
	pdv_heat_problem.time.time_step_length_max = pdv_ns_supg_problem.time.time_step_length_max;


	/* print some info */
	fprintf(Interactive_output, "\n\n**************** BEGINNING TIME INTEGRATION ***************\n");
	fprintf(Interactive_output, "\nTime integration will stop whichever comes first:\n");
	fprintf(Interactive_output, "final_time: %lf\n", time_ns_supg->final_time);
	fprintf(Interactive_output, "final_timestep: %d\n", time_ns_supg->final_step);
	fprintf(Interactive_output, "error less than time_integration_tolerance: %lf\n\n",
			time_ns_supg->conv_meas);

#ifndef PARALLEL
	if (Interactive_input == stdin) {
		printf("Type [Ctrl-C] to manually break time integration.\n");
	}
#endif


	if (pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0) { // mkb interface for solvers

		pdv_ns_supg_heat_current_problem_id = PDC_NS_SUPG_ID;  // ns_supg problem
		int problem_id = pdv_ns_supg_heat_current_problem_id;

		int solver_type = pdr_lins_i_params(PDC_NS_SUPG_ID, 1); // <=0 - direct, >0 - iterative
		// <0 - through direct solver interface, >=0 - through mkb interface

		int max_iter = -1;
		int error_type = -1;
		double error_tolerance = -1;
		int monitoring_level = -1;

		// when no parameter file passed - take control parameters from problem input file
		if (0 == strlen(ctrl_ns_supg->solver_filename)) {

			strcpy(solver_ns_supg_filename, ctrl_ns_supg->solver_filename);
			max_iter = pdr_lins_i_params(problem_id, 2); // max_iter
			error_type = pdr_lins_i_params(problem_id, 3); // error_type
			error_tolerance = pdr_lins_d_params(problem_id, 4); // error_tolerance
			monitoring_level = pdr_lins_i_params(problem_id, 5);  // monitoring level

		}
		else {

			sprintf(solver_ns_supg_filename, "%s/%s",
					ctrl_ns_supg->work_dir, ctrl_ns_supg->solver_filename);

		}

		int parallel = SIC_SEQUENTIAL;
#ifdef PARALLEL
		parallel = SIC_PARALLEL;
#endif

		/*kbw*/
		fprintf(Interactive_output, "initializing ns_supg solver (file %s)\n", solver_ns_supg_filename);
		fprintf(Interactive_output, "parameters: parallel %d, maxiter %d, error_type %d, error_meas %.15lf, monitor %d\n",
				parallel, max_iter,  error_type,  error_tolerance, monitoring_level);
		/*kew*/

		int max_num_grid_levels = 1; // single level solver for ns_supg_heat
		solver_ns_supg_id = sir_init(solver_type, parallel, max_num_grid_levels,
						solver_ns_supg_filename, max_iter,
						error_type, error_tolerance, monitoring_level);
		ctrl_ns_supg->solver_id = solver_ns_supg_id;
		fprintf(Interactive_output, "\nAssigned ns_supg solver ID:  %d\n", solver_ns_supg_id);

	}

	if (pdr_lins_i_params(PDC_HEAT_ID, 1) >= 0) { // mkb interface for solvers

		pdv_ns_supg_heat_current_problem_id = PDC_HEAT_ID;  // heat problem
		int problem_id = pdv_ns_supg_heat_current_problem_id;

		int solver_type = pdr_lins_i_params(PDC_HEAT_ID, 1);   // <=0 - direct, >0 - iterative
		// <0 - through direct solver interface, >=0 - through mkb interface

		int max_iter = -1;
		int error_type = -1;
		double error_tolerance = -1;
		int monitoring_level = -1;

		// when no parameter file passed - take control parameters from problem input file
		if (0 == strlen(ctrl_heat->solver_filename)) {

			strcpy(solver_heat_filename, ctrl_heat->solver_filename);
			max_iter = pdr_lins_i_params(problem_id, 2); // max_iter
			error_type = pdr_lins_i_params(problem_id, 3); // error_type
			error_tolerance = pdr_lins_d_params(problem_id, 4); // error_tolerance
			monitoring_level = pdr_lins_i_params(problem_id, 5);  // monitoring level

		}
		else
			sprintf(solver_heat_filename, "%s/%s", ctrl_heat->work_dir, ctrl_heat->solver_filename);

		int parallel = SIC_SEQUENTIAL;
#ifdef PARALLEL
		parallel = SIC_PARALLEL;
#endif

		/*kbw*/
		fprintf(Interactive_output, "initializing heat solver (file %s)\n", solver_heat_filename);
		fprintf(Interactive_output, "parameters: parallel %d, maxiter %d, error_type %d, error_meas %.15lf, monitor %d\n",
				parallel, max_iter,  error_type,  error_tolerance, monitoring_level);
		/*kew*/

		int max_num_grid_levels = 1; // single level solver for ns_supg_heat
		solver_heat_id = sir_init(solver_type, parallel, max_num_grid_levels,
						solver_heat_filename, max_iter,
						error_type, error_tolerance, monitoring_level);
		ctrl_heat->solver_id = solver_heat_id;

		fprintf(Interactive_output, "\nAssigned heat solver ID:  %d\n", solver_heat_id);
	}


	/* set parameter indicating new mesh */
	iadapt = 1;

	//mmr_init_all_change(1);
	//mmr_init_all_change(2);
	//mmr_copyMESH(1);
	//mmr_init_all_change(2);
	// Writing out step '0'.
	pdr_ns_supg_heat_write_paraview(Work_dir,
			Interactive_input, Interactive_output);

	QElapsedTimer realTime;
	double simulationTime = 0.0;
	double minCurDTime = time_ns_supg->cur_dtime;
	realTime.start();
	/* start loop over time steps */
	while (utv_SIGINT_not_caught && !isInterruptionRequested()) {
		/* update time step and current time */
		++(time_ns_supg->cur_step);
		++(time_heat->cur_step);

		time_ns_supg->prev_dtime = time_ns_supg->cur_dtime;
		time_heat->prev_dtime = time_ns_supg->prev_dtime;

		previous_un_error_heat = sol_norm_un_ns_supg;
		previous_un_error_ns_supg = sol_norm_un_heat;

		//(ileWarstw,obecny_krok,od_krok,ileKrok,minPoprawy,px0,py0,pz0,px1,py1,pz1,endX,endY,endZ);
		//	mmr_test_mesh_motion(1,time_ns_supg->cur_step,3,30,0.005,  0.5,0.5,0.0,   0.7,0.7,1.0,   0.6,0.6,0.0);


		// here possible time step length adaptations !!!

		// compute maximal, minimal and average CFL number for ns_supg_problem
		if (time_ns_supg->CFL_control == 1) {

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
				pcr_allreduce_max_double(1, &cfl_max, &cfl_max_new);
				cfl_max = cfl_max_new;
			}
#endif
			fprintf(Interactive_output,
					"CFL_min = %lf, CFL_max (global) = %lf, CFL_average = %lf\n",
					cfl_min, cfl_max, cfl_ave);

			// check whether CFL is not too big
			if (cfl_max > cfl_limit_max) {
				time_ns_supg->cur_dtime = time_ns_supg->cur_dtime * cfl_decrease_mult;
				time_heat->cur_dtime = time_ns_supg->cur_dtime;
			}
			// check whether CFL is not too small
			if (cfl_max < cfl_limit_min ) {
				time_ns_supg->cur_dtime = time_ns_supg->cur_dtime * cfl_increase_mult;
				time_heat->cur_dtime = time_ns_supg->cur_dtime;
			}
		}


		fprintf(Interactive_output, "\n\n***************************************************************************\n");
		fprintf(Interactive_output, "\n\n***************************************************************************\n");
		fprintf(Interactive_output, "Solving time step %d (step_length: %lf, initial time: %lf)\n",
				time_ns_supg->cur_step, time_ns_supg->cur_dtime, time_ns_supg->cur_time);

		time_ns_supg->cur_time += time_ns_supg->cur_dtime;
		time_heat->cur_time = time_ns_supg->cur_time;

		utr_io_result_add_value_int(RESULT_STEP, time_ns_supg->cur_step);
		utr_io_result_add_value_double(RESULT_CUR_TIME, time_ns_supg->cur_time);
		utr_io_result_add_value_double(RESULT_DT, time_ns_supg->cur_dtime);

		// compute maximal, minimal and average CFL number for ns_supg_problem
		pdr_ns_supg_compute_CFL(PDC_NS_SUPG_ID, &cfl_min, &cfl_max, &cfl_ave);
#ifdef PARALLEL
		double cfl_max_new = cfl_max;
		pcr_allreduce_max_double( 1, &cfl_max, &cfl_max_new);
		cfl_max = cfl_max_new;
#endif


		utr_io_result_add_value_double(RESULT_CFL_MIN, cfl_min);
		utr_io_result_add_value_double(RESULT_CFL_MAX, cfl_max);
		utr_io_result_add_value_double(RESULT_CFL_AVG, cfl_ave);

		fprintf(Interactive_output,
				"CFL_min = %lf, CFL_max = %lf, CFL_average = %lf\n",
				cfl_min, cfl_max, cfl_ave);

		// rewrite solution:
		// current from previous time step becomes previous for current time step
		// soldofs_1 are the most recent solution dofs
		// at the beginning of time step they are rewritten to soldofs_3
		// that holds the solution from the previous time step
		utr_rewr_sol(field_ns_supg_id,
				Current_solution_ID, // 1,
				Previous_time_step_sol_ID // 3
		);	//rewrite current -> u_n (3)
		utr_rewr_sol(field_heat_id,
				Current_solution_ID, // 1,
				Previous_time_step_sol_ID // 3
		);	//rewrite current -> u_n (3)
		//pdr_rewr_heat_dtdt_sol(field_heat_dtdt_id, 1, 2);	//rewrite current (1) to previous (2)
		/* update time dependent boundary conditions for NS problem*/
		pdv_ns_supg_heat_current_problem_id = PDC_NS_SUPG_ID;
		pdr_ns_supg_update_timedep_bc(&problem_ns_supg->bc,
				time_ns_supg->cur_time);
		/* update time dependent boundary conditions for thermal problem */
		pdv_ns_supg_heat_current_problem_id = PDC_HEAT_ID;
		pdr_heat_update_timedep_bc(&pdv_heat_problem.bc, time_heat->cur_time);

		nonl_iter = 0;
		for (;;) {

			utr_io_result_add_value_int(RESULT_NON_LIN_STEP, nonl_iter + 1);

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
			utr_rewr_sol(field_heat_id,
					Current_solution_ID, // 1,
					Previous_iteration_sol_ID // 2
			);	//rewrite current -> u_k (2)

			if (iadapt == 1) {
#ifdef PARALLEL
				int nr_levels = 1; // single level solver for ns_supg_heat
				// 4 DOFs for NS_SUPG
				pdv_ns_supg_exchange_table_index = appr_init_exchange_tables(pcr_nr_proc(),
								pcr_my_proc_id(),
								field_ns_supg_id,
								0, 4, nr_levels);
				// 1 DOF for heat
				pdv_heat_exchange_table_index = appr_init_exchange_tables(pcr_nr_proc(),
								pcr_my_proc_id(),
								field_heat_id,
								0, 1, nr_levels);
#endif
				if (pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0) { // mkb interface for solvers
					sir_create(solver_ns_supg_id, PDC_NS_SUPG_ID);
				}
				if (pdr_lins_i_params(PDC_HEAT_ID, 1) >= 0) { // mkb interface for solvers
					sir_create(solver_heat_id, PDC_HEAT_ID);
				}
				iadapt = 0;
			}

			if (pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0) { // mkb interface for solvers

				int problem_id = PDC_NS_SUPG_ID;
				int nr_iter = -1;
				double conv_meas = -1.0;
				double conv_rate = -1.0;
				int monitor = -1;

				// when no parameter file passed - take control parameters from problem input file
				if (0 == strlen(ctrl_ns_supg->solver_filename)) {
					nr_iter = pdr_lins_i_params(problem_id, 2); // max_iter
					conv_meas = pdr_lins_d_params(problem_id, 4); // error_tolerance
					monitor = pdr_lins_i_params(problem_id, 5);  // monitoring level
				}

				int ini_guess = 1; // get initial guess from data structure

				/*---------- CALLING DIRECT OR ITERATIVE SOLVER THROUGH MKB INTERFACE ---------*/
				sir_solve(solver_ns_supg_id, SIC_SOLVE, ini_guess, monitor,
						&nr_iter, &conv_meas, &conv_rate);
				if (pdr_lins_i_params(problem_id, 1) > 0) { // iterative solver
					fprintf(Interactive_output,
							"\nAfter %d iterations of linear solver for ns_supg problem\n",
							nr_iter);
					fprintf(Interactive_output,
							"Convergence measure: %lf, convergence rate %lf\n",
							conv_meas, conv_rate);
				}

			} else {
#ifdef PARALLEL
				if (pcv_nr_proc > 1) {
					mf_log_err("Shared memory direct linear solver called in distributed memory computations!");
				}
#endif

				sir_direct_solve_lin_sys(PDC_NS_SUPG_ID, SIC_SEQUENTIAL, solver_ns_supg_filename);
			}

			if (pdr_lins_i_params(PDC_HEAT_ID, 1) >= 0) { // mkb interface for solvers

				int problem_id = PDC_HEAT_ID;
				int nr_iter = -1;
				double conv_meas = -1.0;
				double conv_rate = -1.0;
				int monitor = -1;

				// when no parameter file passed - take control parameters from problem input file
				if (0 == strlen(ctrl_heat->solver_filename)) {
					nr_iter = pdr_lins_i_params(problem_id, 2); // max_iter
					conv_meas = pdr_lins_d_params(problem_id, 4); // error_tolerance
					monitor = pdr_lins_i_params(problem_id, 5);  // monitoring level
				}

				int ini_guess = 1; //  get initial guess from data structure

				/*---------- CALLING DIRECT OR ITERATIVE SOLVER THROUGH MKB INTERFACE ---------*/
				sir_solve(solver_heat_id, SIC_SOLVE, ini_guess, monitor,
						&nr_iter, &conv_meas, &conv_rate);
				if (pdr_lins_i_params(problem_id, 1) > 0) { // iterative solver
					fprintf(Interactive_output,
							"\nAfter %d iterations of linear solver for heat problem\n",
							nr_iter);
					fprintf(Interactive_output,
							"Convergence measure: %lf, convergence rate %lf\n",
							conv_meas, conv_rate);
				}

			} else {

#ifdef PARALLEL
				if (pcv_nr_proc > 1) {
					mf_log_err("Shared memory PARDISO called for solving in distributed memory computations!");
				}
#endif

				sir_direct_solve_lin_sys(PDC_HEAT_ID, SIC_SEQUENTIAL,
						solver_heat_filename);
			}


			nsSupgHeatSolDiffNorm(Current_solution_ID, Previous_iteration_sol_ID, &sol_norm_uk_ns_supg, &sol_norm_uk_heat);

			if (nonl_iter == 0) {
				initial_nonl_error_ns_supg = sol_norm_uk_ns_supg;
				initial_nonl_error_heat = sol_norm_uk_heat;
			}
			fprintf(Interactive_output, "\nAfter linear solver in nonlinear iteration %d\n",
					nonl_iter);
			fprintf(Interactive_output, "NS_SUPG solution difference norm (u_k, prev. u_k): %lf (limit %lf)\n",
					sol_norm_uk_ns_supg, nonl_ns_supg->conv_meas);
			fprintf(Interactive_output, "Heat solution difference norm (u_k, prev. u_k): %lf (limit %lf)\n",
					sol_norm_uk_heat, nonl_heat->conv_meas);

			utr_io_result_add_value_double(RESULT_NORM_NS_SUPG_NONL, sol_norm_uk_ns_supg);
			utr_io_result_add_value_double(RESULT_NORM_HEAT_NONL, sol_norm_uk_heat);

			++nonl_iter;

			int break_indi = 0;
			if ( nonl_ns_supg->conv_type == 0 ) { // && nonl_heat->conv_type==0){

				if ( (sol_norm_uk_ns_supg < nonl_ns_supg->conv_meas)
						&& (sol_norm_uk_heat < nonl_heat->conv_meas) ) {

					fprintf(Interactive_output,
							"\nConvergence in nonlinear iterations!\n");

					//break;
					break_indi = 1;
				}

			}
			else if ( nonl_ns_supg->conv_type == 1 ) { // && nonl_heat->conv_type==1){

				if ( (sol_norm_uk_ns_supg < nonl_ns_supg->conv_meas * initial_nonl_error_ns_supg)
						&& (sol_norm_uk_heat < nonl_heat->conv_meas * initial_nonl_error_heat) ) {

					fprintf(Interactive_output,
							"\nConvergence in nonlinear iterations!\n");

					//break;
					break_indi = 1;
				}

			}

			if (nonl_iter >= nonl_ns_supg->max_iter
					&& nonl_iter >= nonl_heat->max_iter) {
				fprintf(Interactive_output,
						"\nMax nonlinear iterations (ns_supg %d and heat %d) reached - breaking.\n",
						nonl_ns_supg->max_iter, nonl_heat->max_iter);

				//break;
				break_indi = 1;
			}


			if (break_indi == 1) {
				break;
			}


		} // the end of infinite loop over non-linear iterations

		nsSupgHeatSolDiffNorm(Current_solution_ID, Previous_time_step_sol_ID, &sol_norm_un_ns_supg, &sol_norm_un_heat);
		fprintf(Interactive_output, "\n---------------------------------------------------------------------------\n");
		fprintf(Interactive_output, "After non-linear solver in time step %d\n",
				time_ns_supg->cur_step);
		fprintf(Interactive_output, "NS_SUPG solution difference norm (u_n, prev. u_n): %lf (limit %lf)\n",
				sol_norm_un_ns_supg, time_ns_supg->conv_meas);
		fprintf(Interactive_output, "Heat solution difference norm (u_n, prev. u_n): %lf (limit %lf)\n",	    sol_norm_un_heat, time_heat->conv_meas);

		utr_io_result_add_value_double(RESULT_NORM_NS_SUPG, sol_norm_un_ns_supg);
		utr_io_result_add_value_double(RESULT_NORM_HEAT, sol_norm_un_heat);
		//    utr_io_result_add_value_double(RESULT_NORM_NS_SUPG,time_ns_supg->conv_meas);
		//    utr_io_result_add_value_double(RESULT_NORM_HEAT,time_heat->conv_meas);
		char mm_name[255] = {0};
		mmr_module_introduce(mm_name);
		//printf("mesh module %s\n", mm_name);
		if ( 0 != strcmp("3D_remesh", mm_name) ) {
			//printf("mesh module different than 3D_remesh\n");
			utr_ctrl_pts_add_values(pdv_ns_supg_problem.ctrl.field_id);
			utr_ctrl_pts_add_values(pdv_heat_problem.ctrl.field_id);
		}

		// check whether we should change time step length
		if (time_ns_supg->time_step_length_nonl_control == 1) {
			// && time_heat->time_step_length_nonl_control == 1){

			fprintf(Interactive_output,
					"In time step control (iter_min %d <=? nonl_iter %d <=? iter_max %d)\n",
					time_ns_supg->time_step_length_nonl_iter_min, nonl_iter,
					time_ns_supg->time_step_length_nonl_iter_max );

			fprintf(Interactive_output,
					"In time step control ns_supg (norm_min %lf <=? norm %lf <=? norm_max %lf)\n",
					time_ns_supg->time_step_length_min_norm, sol_norm_un_ns_supg,
					time_ns_supg->time_step_length_max_norm );

			fprintf(Interactive_output,
					"In time step control heat (norm_min %lf <=? norm %lf <=? norm_max %lf)\n",
					time_heat->time_step_length_min_norm, sol_norm_un_heat,
					time_heat->time_step_length_max_norm );

			/* time step length adaptation based on number of non-linear iterations */
			//int time_step_length_control = time_ns_supg->time_step_length_nonl_control;
			double time_step_length_limit_max = time_ns_supg->time_step_length_nonl_iter_max;
			double time_step_length_limit_min = time_ns_supg->time_step_length_nonl_iter_min;
			double time_step_length_increase_mult =
					time_ns_supg->time_step_length_nonl_iter_increase_mult;
			double time_step_length_decrease_mult =
					time_ns_supg->time_step_length_nonl_iter_decrease_mult;

			// if convergence was slow - we decrease time step length
			if ( (nonl_iter >= time_step_length_limit_max)
					|| (sol_norm_un_ns_supg > time_ns_supg->time_step_length_max_norm)
					|| (sol_norm_un_heat > time_heat->time_step_length_max_norm) ) {
				time_ns_supg->cur_dtime = time_ns_supg->cur_dtime * time_step_length_decrease_mult;
				time_heat->cur_dtime = time_ns_supg->cur_dtime;
			}
			// if convergence is fast we increase time step length
			else if ( (nonl_iter <= time_step_length_limit_min)
					|| (sol_norm_un_ns_supg < time_ns_supg->time_step_length_min_norm)
					|| (sol_norm_un_heat < time_heat->time_step_length_min_norm) ) {
				time_ns_supg->cur_dtime = time_ns_supg->cur_dtime * time_step_length_increase_mult;
				time_heat->cur_dtime = time_ns_supg->cur_dtime;
			}

			if (time_ns_supg->cur_dtime < time_ns_supg->time_step_length_min) {
				time_ns_supg->cur_dtime = time_ns_supg->time_step_length_min;
				mf_log_info("Current time step(%lf) is below minimum, resetting to time_step_length_min(%lf)",
						time_ns_supg->cur_dtime, time_ns_supg->time_step_length_min);
			}
			else if (time_ns_supg->cur_dtime > time_ns_supg->time_step_length_max) {
				time_ns_supg->cur_dtime = time_ns_supg->time_step_length_max;
				mf_log_info("Current time step(%lf) exceeds maximum, resetting to time_step_length_max(%lf)",
						time_ns_supg->cur_dtime, time_ns_supg->time_step_length_max);
			}

		} //end of time step adapt.


		// AS: heating-cooling problem dofs shift 2->3, 1->2
		//pdr_rewr_heat_dtdt_sol(field_heat_dtdt_id, 2, 3);	//rewrite previous (2) to early (3)
		//pdr_rewr_heat_dtdt_sol(field_heat_dtdt_id, 1, 2);	//rewrite current (1) to previous (2)
		//pdr_heat_dtdt(1,3);  // compute and write Current=1 heat_dtdt field

		/* graphics data dumps */
		if (time_ns_supg->intv_graph > 0 &&
				time_ns_supg->cur_step % time_ns_supg->intv_graph == 0) {

			mf_log_info("Writing field to disk (graphics)...");

			pdr_ns_supg_heat_write_paraview(Work_dir,
					Interactive_input, Interactive_output);
		}

		/* full data dumps (for restarting) */
		if (time_ns_supg->intv_dumpout > 0
				&& time_ns_supg->cur_step % time_ns_supg->intv_dumpout == 0) {

			mf_log_info("Writing field to disk (full precision)...");

			pdr_ns_supg_heat_dump_data(Work_dir,
					Interactive_input, Interactive_output);

		}

		utr_io_result_write_values_and_proceed();

		/* check for stop conditions */

		// stop if computational error occured
		if ( pdv_ns_supg_problem.ctrl.error_indicator != NO_ERROR ) {
			mf_log_err("Breaking time integration because of error(%d).",
					pdv_ns_supg_problem.ctrl.error_indicator);
			break;
		}

		/* stop if convergence reached (for stationary problems) */
		if (sol_norm_un_ns_supg < time_ns_supg->conv_meas
				&& sol_norm_un_heat < time_heat->conv_meas) {
			fprintf(Interactive_output,
					"\nConvergence in time integration!\n");
			break;
		}

		// when time for adaptation
		if ( (adpt_ns_supg->type > 0 && adpt_ns_supg->interval > 0 &&
						(time_ns_supg->cur_step + 1) % adpt_ns_supg->interval == 0)
				|| (adpt_heat->type > 0 && adpt_heat->interval > 0 &&
						(time_heat->cur_step + 1) % adpt_heat->interval == 0)) {

#ifdef PARALLEL
			/* free exchange tables for DOFs - for both fields: ns_supg and heat */
			// in reverse order
			appr_free_exchange_tables(pdv_heat_exchange_table_index);
			appr_free_exchange_tables(pdv_ns_supg_exchange_table_index);
#endif

			// to satisfy strange requirement of deallocating matrices in LIFO order
			// heat matrix is deallocated first
			if (pdr_lins_i_params(PDC_HEAT_ID, 1) >= 0) { // mkb interface for solvers
				sir_free(solver_heat_id);
			}
			if (pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0) { // mkb interface for solvers
				/* free solver data structures */
				sir_free(solver_ns_supg_id);
			}


			pdr_ns_supg_heat_adapt(Work_dir, Interactive_input, Interactive_output);
			/* indicate to recreate block structure */
			iadapt = 1;
		}

		simulationTime += time_ns_supg->cur_dtime;
		double timeDiff = static_cast<double>(realTime.elapsed() / 1000.0) - simulationTime;
		CUTEHMI_DEBUG("Difference between real time and simulation time: " << timeDiff);
		time_ns_supg->cur_dtime += timeDiff;

		time_ns_supg->cur_dtime = std::max(time_ns_supg->cur_dtime, minCurDTime);
		time_heat->cur_dtime = time_ns_supg->cur_dtime;
		CUTEHMI_DEBUG("Setting up new cur_dtime: " << time_ns_supg->cur_dtime);

		emit iterationFinished();
	}  //end loop over timesteps

	if (iadapt == 0) {
#ifdef PARALLEL
		/* free exchange tables for DOFs - for both fields: ns_supg and heat */
		// in reverse order
		appr_free_exchange_tables(pdv_heat_exchange_table_index);
		appr_free_exchange_tables(pdv_ns_supg_exchange_table_index);
#endif

		// to satisfy strange requirement of deallocating matrices in LIFO order
		// heat matrix is deallocated first
		if (pdr_lins_i_params(PDC_HEAT_ID, 1) >= 0) { // mkb interface for solvers
			sir_free(solver_heat_id);
		}
		if (pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0) { // mkb interface for solvers
			sir_free(solver_ns_supg_id);
		}
	}

	// to satisfy strange requirement of deallocating matrices in LIFO order
	// heatsolver is destroyed first
	if (pdr_lins_i_params(PDC_HEAT_ID, 1) >= 0) { // mkb interface for solvers
		sir_destroy(solver_heat_id);
	}
	if (pdr_lins_i_params(PDC_NS_SUPG_ID, 1) >= 0) { // mkb interface for solvers
		sir_destroy(solver_ns_supg_id);
	}

	fprintf(Interactive_output, "\n**************** LEAVING TIME INTEGRATION ***************\n");
}

void IntegrationThread::nsSupgHeatSolDiffNorm(int Current, int Old, double * sol_diff_norm_ns_supg_p, double * sol_diff_norm_heat_p)
{
	double sol_dofs_current[APC_MAXELSD];	/* solution dofs */
	double sol_dofs_old[APC_MAXELSD];	/* solution dofs */
	int field_id, mesh_id;
	int node_id = 0;
	double temp, norm_ns_supg = 0.0, norm_heat = 0.0;
	int i;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	pdv_ns_supg_heat_current_problem_id = PDC_NS_SUPG_ID;  // ns_supg problem
	i = 3;
	field_id = pdr_ctrl_i_params(pdv_ns_supg_heat_current_problem_id, i);
	mesh_id = apr_get_mesh_id(field_id);

	while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
		if (apr_get_ent_pdeg(field_id, APC_VERTEX, node_id) > 0) {
			int numdofs = apr_get_ent_nrdofs(field_id, APC_VERTEX, node_id);
			apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,
					Current, sol_dofs_current);
			apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,
					Old, sol_dofs_old);
			i = 0; // for ns_supg we take into account velocities ONLY!!!
			for (i = 0; i < numdofs - 1; ++i) {

				/*kbw
				  double norm=0.0;
				  printf("node_id %d\ti %d\tcurrent %lf\told %lf\tnorm %lf\n",
				  node_id,i,sol_dofs_current[i],sol_dofs_old[i],norm);
				  norm += (sol_dofs_current[i] - sol_dofs_old[i]) * (sol_dofs_current[i] - sol_dofs_old[i]);
				/*kew*/

				temp = fabs(sol_dofs_current[i] - sol_dofs_old[i]);
				if (norm_ns_supg < temp) norm_ns_supg = temp;
			}
		}
	}

	pdv_ns_supg_heat_current_problem_id = PDC_HEAT_ID;  // heat problem
	i = 3;
	field_id = pdr_ctrl_i_params(pdv_ns_supg_heat_current_problem_id, i);
	mesh_id = apr_get_mesh_id(field_id);

	while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
		if (apr_get_ent_pdeg(field_id, APC_VERTEX, node_id) > 0) {
			int numdofs = apr_get_ent_nrdofs(field_id, APC_VERTEX, node_id);
			apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,
					Current, sol_dofs_current);
			apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,
					Old, sol_dofs_old);
			i = 0;
			for (i = 0; i < numdofs; ++i) {
				temp = fabs(sol_dofs_current[i] - sol_dofs_old[i]);
				if (norm_heat < temp)
					norm_heat = temp;
			}
		}
	}

#ifdef PARALLEL
	pcr_allreduce_max_double(1, &norm_ns_supg, &temp);
	norm_ns_supg = temp;
	pcr_allreduce_max_double(1, &norm_heat, &temp);
	norm_heat = temp;
#endif

	*sol_diff_norm_ns_supg_p = norm_ns_supg;
	*sol_diff_norm_heat_p = norm_heat;
}

}
}
