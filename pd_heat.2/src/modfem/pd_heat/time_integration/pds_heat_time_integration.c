/************************************************************************
File pds_heat_time_integration.c - time stepping procedures

Contains definition of routines:
  pdr_heat_comp_sol_diff_norm - returns max difference between
										current and old solutions vectors
  pdr_heat_time_integration - time integration driver

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

/* USED DATA STRUCTURES AND INTERFACES FROM OTHER MODULES */

/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>		/* USES */
/* interface for all approximation modules */
#include <modfem/aph_intf.h>		/* USES */
/* interface for all solver modules */
#include <modfem/sih_intf.h>		/* USES */


#include <modfem/uth_log.h>
#include <modfem/uth_system.h>

#ifdef PARALLEL
/* interface of parallel mesh manipulation modules */
#include "mmph_intf.h"		/* USES */
/* interface for all parallel approximation modules */
#include "apph_intf.h"		/* USES */
/* interface for parallel communication modules */
#include "pch_intf.h"		/* USES */
#endif

/* utilities - including simple time measurement library */
#include <modfem/uth_intf.h>

/* interface for thread management modules */
#include <modfem/tmh_intf.h>


/* USED AND IMPLEMENTED PROBLEM DEPENDENT DATA STRUCTURES AND INTERFACES */

/* problem dependent module interface */
#include <modfem/pdh_intf.h>		/* USES */

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>	/* USES */

/* heat problem module's types and functions */
#include <modfem/pd_heat/pdh_heat.h> /* USES & IMPLEMENTS */
#include <modfem/pd_heat/pdh_heat_weakform.h>
#include <modfem/pd_heat/pdh_heat_problem.h>
// bc and material header files are included in problem header files


/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */

/*------------------------------------------------------------
pdr_heat_sol_diff_norm - returns max difference in current and old
			  solutions vectors
------------------------------------------------------------*/
void pdr_heat_sol_diff_norm(
		int Current, // in: current vector id
		int Old,     // in: old vector id
		double* sol_diff_norm_heat_p // norm of difference current-old for heat problem
		)
{

	double sol_dofs_current[APC_MAXELSD];	/* solution dofs */
	double sol_dofs_old[APC_MAXELSD];	/* solution dofs */
	int field_id, mesh_id;
	int node_id = 0;
	double temp, norm_heat = 0.0;
	int i;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	pdv_heat_current_problem_id = PDC_HEAT_ID;
	i=3; field_id = pdr_ctrl_i_params(pdv_heat_current_problem_id,i);
	mesh_id = apr_get_mesh_id(field_id);

	node_id = 0;
	while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
		if (apr_get_ent_pdeg(field_id, APC_VERTEX, node_id) > 0) {
			int numdofs = apr_get_ent_nrdofs(field_id, APC_VERTEX, node_id);
			apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,
							  Current, sol_dofs_current);
			apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,
							  Old, sol_dofs_old);
			i = 0;
			for (i = 0; i < numdofs; ++i) {

				/*kbw
  printf("el_id %d\ti %d\tcurrent %lf\told %lf\tnorm %lf\n",el_id,i,sol_dofs_current[i],sol_dofs_old[i],norm);
  norm += (sol_dofs_current[i] - sol_dofs_old[i]) * (sol_dofs_current[i] - sol_dofs_old[i]);
/*kew*/

				temp = fabs(sol_dofs_current[i] - sol_dofs_old[i]);
				if (norm_heat < temp){
					norm_heat = temp;
				}

			}
		}
	}


#ifdef PARALLEL
	pcr_allreduce_max_double(1, &norm_heat, &temp);
	norm_heat = temp;
#endif

	*sol_diff_norm_heat_p = norm_heat;

	return;
}

/* /\*------------------------------------------------------------ */
/*   pdr_rewr_heat_dtdt_sol - to rewrite solution from one vector to another */
/* ------------------------------------------------------------*\/ */
/* int pdr_rewr_heat_dtdt_sol( */
/*   int Field_id,      /\* in: data structure to be used  *\/ */
/*   int Sol_from,      /\* in: ID of vector to read solution from *\/ */
/*   int Sol_to         /\* in: ID of vector to write solution to *\/ */
/*   ) */
/* { */

/*   int nel, mesh_id, num_eq; */
/*   double dofs_loc[APC_MAXELSD]; */

/* /\*++++++++++++++++ executable statements ++++++++++++++++*\/ */

/*   mesh_id = apr_get_mesh_id(Field_id); */

/*   nel=0; */
/*   num_eq = pdv_heat_dtdt_problem.ctrl.nreq; */
/*   while((nel=mmr_get_next_act_elem(mesh_id,nel))!=0){ */
/*     printf("\n\tAS: nel = %d", nel); */
/*     apr_read_ent_dofs(Field_id, APC_ELEMENT, nel, num_eq, Sol_from, dofs_loc); */
/*     apr_write_ent_dofs(Field_id, APC_ELEMENT, nel, num_eq, Sol_to, dofs_loc); */

/*   } */

/*   return(0); */
/* } */

/* /\*------------------------------------------------------------ */
/* pdr_heat_dtdt - returns heating-cooling rate field between */
/* 		    current and old heat solution vector */
/* ------------------------------------------------------------*\/ */
/* void pdr_heat_dtdt( */
/*   int Current, // in: current vector id  */
/*   int Old     // in: old vector id */
/*   ) */
/* { */

/*   double sol_dofs_current[APC_MAXELSD];	/\* solution dofs *\/ */
/*   double sol_dofs_old[APC_MAXELSD];	/\* solution dofs *\/ */
/*   double sol_dofs_dtdt[APC_MAXELSD];	/\* solution dofs *\/ */
/*   int field_id, mesh_id; */
/*   int field_dtdt_id; */
/*   int node_id = 0; */
/*   double dtdt;; */
/*   int i, j; */
/*   int num_eq; */
/*   pdt_heat_times *time_heat = &pdv_heat_problem.time; */
/*   pdt_heat_material_query_params query_params; */
/*   pdt_heat_material_query_result query_result; */
/*   double VOF_t1, VOF_t2, latent_heat_of_fusion, specific_heat, specific_heat_t2, specific_heat_t1; */
/*   pdt_heat_problem *problem_heat = &pdv_heat_problem; */

/* /\*++++++++++++++++ executable statements ++++++++++++++++*\/ */

/*   mesh_id = apr_get_mesh_id(field_id); */
/*   pdv_heat_current_problem_id = PDC_HEAT_DTDT_ID;  // heating-cooling problem */
/*   i=3; field_dtdt_id = pdr_ctrl_i_params(pdv_heat_current_problem_id,i); */
/*   //printf("\n\tAS: field_heat_id = %d\tfield_dtdt_id = %d\n", field_id, field_dtdt_id); */

// CHECK WHETHER MATERIAL DATABASE USED AT ALL
/*   query_params.material_idx = 0;	//material by idx */
/*   query_params.name = ""; */
/*   query_params.temperature = pdv_heat_problem.ctrl.ref_temperature; */
/*   pdr_heat_material_query(&problem_heat->materials,  */
/* 			  &query_params, &query_result); */
/*   // latent_heat_of_fusion - temperature independent */
/*   latent_heat_of_fusion = query_result.latent_heat_of_fusion; */
/*   //printf("\n\tAS: ref_specific_heat = %lf\tlatent_heat_of_fusion = %lf\n", specific_heat, latent_heat_of_fusion); */
/*   node_id = 0; */
/*   while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) { */
/*     if (apr_get_ent_pdeg(field_id, APC_VERTEX, node_id) > 0) { */
/*       int numdofs = apr_get_ent_nrdofs(field_id, APC_VERTEX, node_id); */
/*       apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,  */
/* 			Current, sol_dofs_current); */
/*       apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs,  */
/* 			Old, sol_dofs_old); */

/*       //for (i = 0; i < numdofs; ++i) { */

/*       j = 0; // temperature numdofs = 1 */
/*       i = 0; // latent heat of fusion coefficient */
/*       // get material data at two instants */
/*       query_params.temperature = sol_dofs_old[j]; */
/*       pdr_heat_material_query(&problem_heat->materials,  */
/* 				  &query_params, &query_result); */
/*       VOF_t1 = query_result.VOF; */
/*       specific_heat_t1 = query_result.specific_heat; */

/*       query_params.temperature = sol_dofs_current[j]; */
/*       pdr_heat_material_query(&problem_heat->materials,  */
/* 				  &query_params, &query_result); */
/*       VOF_t2 = query_result.VOF; */
/*       specific_heat_t2 = query_result.specific_heat; */

/*       specific_heat = (specific_heat_t1 + specific_heat_t2) / 2.0; */
/*       if ( VOF_t1 != VOF_t2 ) { */
/* 	sol_dofs_dtdt[i] = -1.0*(latent_heat_of_fusion / specific_heat) * (VOF_t2 - VOF_t1) / time_heat->cur_dtime; */
/*       } else { */
/* 	sol_dofs_dtdt[i] = 0.0; */
/*       } */

/*       i = 1; // heating-cooling rate */
/*       sol_dofs_dtdt[i] = (sol_dofs_current[j] - sol_dofs_old[j]) / time_heat->cur_dtime; */

/*       //if ( VOF_t1 != VOF_t2 ) { */
/*       /\* if ( (sol_dofs_current[j]>=query_result.temp_solidus) || (sol_dofs_old[j]>=query_result.temp_solidus) ) { *\/ */
/*       /\* 	printf("\n\tAS: node_id = %d\tnumdofs(PDC_HEAT_ID) = %d\tsol_dofs_dtdt[0] = %lf\tsol_dofs_dtdt[1] = %lf",  *\/ */
/*       /\* 	       node_id, numdofs, sol_dofs_dtdt[0], sol_dofs_dtdt[1]); *\/ */
/*       /\* 	printf("\n\tAS: Lf * dfL/dt / Cp = %lf", -1.0*(latent_heat_of_fusion / specific_heat) * (VOF_t2 - VOF_t1) / time_heat->cur_dtime); *\/ */
/*       /\* 	printf("\n\tAS: VOF_t1 = %lf\tVOF_t2 = %lf", VOF_t1, VOF_t2); *\/ */
/*       /\* 	printf("\n\tAS: dT/dt = %lf", (sol_dofs_current[j] - sol_dofs_old[j]) / time_heat->cur_dtime); *\/ */
/*       /\* 	printf("\n\tAS: T_t1 = %lf\tT_t2 = %lf\n", sol_dofs_old[j], sol_dofs_current[j]); *\/ */
/*       /\* } *\/ */

/*       //} */

/*       //printf("\n\tAS: node_id = %d\tsol_dofs_dtdt[0] = %lf\tsol_dofs_dtdt[1] = %lf\n", node_id, sol_dofs_dtdt[0], sol_dofs_dtdt[1]); */

/*       num_eq = pdv_heat_dtdt_problem.ctrl.nreq; */
/*       //num_shap = apr_get_ent_numshap(field_id, APC_VERTEX, node_id); */
/*       //printf("\n\tAS:num_shap = %d", num_shap); */
/*       //apr_write_ent_dofs(field_dtdt_id, APC_VERTEX, node_id, num_shap, Current, sol_dofs_dtdt); */
/*       //AS: num_shap = i+1, nreq = 2 (1 - latent heat factor, 2 - heating-cooling rate */
/*       apr_write_ent_dofs(field_dtdt_id, APC_VERTEX, node_id, num_eq, Current, sol_dofs_dtdt); */

/*       //i = 2, 3, etc. for other properties */
/*     } */
/*   } */

/*   return; */
/* } */

/*------------------------------------------------------------
pdr_heat_time_integration - time integration driver
------------------------------------------------------------*/
void pdr_heat_time_integration(
		char* Work_dir,
		FILE *Interactive_input,
		FILE *Interactive_output
		)
{


	/* one instance of mesh, two instances of fields */
	int mesh_id;
	int field_heat_id;
	//int field_heat_dtdt_id;
	int  solver_heat_id;
	double sol_norm_uk_heat, sol_norm_un_heat;
	char autodump_filename[300];
	int iadapt = 0;
	char solver_heat_filename[300];;
	int nr_iter, nonl_iter, monitor;
	double conv_meas, conv_rate;
	int i, iaux;
	//double cfl_min=1000000, cfl_max=0, cfl_ave=0.0;


	/*  heat problem parameters */
	pdt_heat_problem *problem_heat = &pdv_heat_problem;
	pdt_heat_ctrls *ctrl_heat = &pdv_heat_problem.ctrl;
	pdt_heat_times *time_heat = &pdv_heat_problem.time;
	pdt_heat_nonls *nonl_heat = &pdv_heat_problem.nonl;
	pdt_heat_linss *lins_heat = &pdv_heat_problem.lins;
	pdt_heat_adpts *adpt_heat = &pdv_heat_problem.adpt;

	/*  heating-cooling problem parameters */
	//pdt_heat_dtdt_problem *problem_heat_dtdt = &pdv_heat_dtdt_problem;
	//pdt_heat_dtdt_ctrls *ctrl_heat_dtdt = &pdv_heat_dtdt_problem.ctrl;

	/* time adaptation */
	double initial_dtime = time_heat->prev_dtime;
	double final_dtime = time_heat->cur_dtime;
	double dtime_adapt, time_adapt;
	double daux;


	/*++++++++++++++++ executable statements ++++++++++++++++*/


	// get mesh and field parameters
	i = 2; mesh_id = pdr_ctrl_i_params(pdv_heat_current_problem_id,i);

	field_heat_id = pdr_ctrl_i_params(pdv_heat_current_problem_id,3);
	//printf("\n\tAS: field_heat_id = %d", field_heat_id);
	//field_heat_dtdt_id = pdr_ctrl_i_params(pdv_heat_current_problem_id,3);
	//printf("\n\tAS: field_heat_id = %d", field_heat_id);
	//printf("\n\tAS: field_heat_dtdt_id = %d\n", field_heat_dtdt_id);


	if (time_heat->cur_time >= time_heat->final_time) {

		fprintf(Interactive_output,
				"\nCurrent time: %lf is bigger than final time: %lf.\n\n",
				time_heat->cur_time, time_heat->final_time);

		fprintf(Interactive_output,
				"\nCurrent time-step: %d, current time step length: %lf.\n\n",
				time_heat->cur_step, time_heat->cur_dtime);

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
				printf("Set new time step length:\n");
				//printf("Set new time step length (or CFL limit if < 0):\n");
				scanf("%lf",&daux);
#ifdef PARALLEL
			}
			pcr_bcast_double(pcr_print_master(), 1, &daux);
#endif
			time_heat->cur_dtime = daux;
			/*       if(daux>0){ */
			/* 	time_heat->CFL_control = 0; */
			/* 	time_heat->cur_dtime = daux; */
			/*       } */
			/*       else{ */
			/* 	time_heat->CFL_control = 1; */
			/* 	time_heat->CFL_limit = -daux; */

			/* #ifdef PARALLEL */
			/* 	if (pcr_my_proc_id() == pcr_print_master()) { */
			/* #endif */
			/* 	  printf("Set new time step length multiplier:\n"); */
			/* 	  scanf("%lf",&daux); */
			/* #ifdef PARALLEL */
			/* 	} */
			/* 	pcr_bcast_double(pcr_print_master(), 1, &daux); */
			/* #endif */

			/* 	time_heat->time_step_length_mult = daux; */

			/*       } */

		} else{

			fprintf(Interactive_output, "\nExiting!\n\n");
			exit(0);;

		}

		time_heat->final_step = time_heat->cur_step + iaux;

		time_heat->final_time = time_heat->cur_time +
				iaux*time_heat->cur_dtime;


	}

	if (time_heat->cur_step >= time_heat->final_step) {

		fprintf(Interactive_output,
				"\nCurrent time-step: %d is bigger than final step: %d\n\n",
				time_heat->cur_step, time_heat->final_step);

		fprintf(Interactive_output,
				"\nCurrent time: %lf, current time step length: %lf.\n\n",
				time_heat->cur_time, time_heat->cur_dtime);

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
				//printf("Set new time step length (or CFL limit if < 0):\n");
				printf("Set new time step length:\n");
				scanf("%lf",&daux);
#ifdef PARALLEL
			}
			pcr_bcast_double(pcr_print_master(), 1, &daux);
#endif
			time_heat->cur_dtime = daux;
			/*       if(daux>0){ */
			/* 	time_heat->CFL_control = 0; */
			/* 	time_heat->cur_dtime = daux; */
			/*       } */
			/*       else{ */
			/* 	time_heat->CFL_control = 1; */
			/* 	time_heat->CFL_limit = -daux; */

			/* #ifdef PARALLEL */
			/* 	if (pcr_my_proc_id() == pcr_print_master()) { */
			/* #endif */

			/* 	  printf("Set new time step length multiplier:\n"); */
			/* 	  scanf("%lf",&daux); */

			/* #ifdef PARALLEL */
			/* 	} */
			/* 	pcr_bcast_double(pcr_print_master(), 1, &daux); */
			/* #endif */

			/* 	time_heat->time_step_length_mult = daux; */

			/*       } */

		} else{

			fprintf(Interactive_output, "\nExiting!\n\n");
			exit(0);;

		}

		time_heat->final_step = time_heat->cur_step + iaux;

		time_heat->final_time = time_heat->cur_time +
				iaux*time_heat->cur_dtime;

	}

	/* print some info */
	fprintf(Interactive_output, "\n\n**************** BEGINNING TIME INTEGRATION ***************\n");
	fprintf(Interactive_output, "\nTime integration will stop whichever comes first:\n");
	fprintf(Interactive_output, "final_time: %lf\n", time_heat->final_time);
	fprintf(Interactive_output, "final_timestep: %d\n", time_heat->final_step);
	fprintf(Interactive_output, "error less than time_integration_tolerance: %lf\n\n", time_heat->conv_meas);

#ifndef PARALLEL
	if(Interactive_input == stdin) {
		printf("Type [Ctrl-C] to manually break time integration.\n");
	}
#endif


	if(pdr_lins_i_params(pdv_heat_current_problem_id, 1) >= 0){ // mkb interface for solvers

		int problem_id = pdv_heat_current_problem_id;
		int solver_type = pdr_lins_i_params(pdv_heat_current_problem_id, 1);
		int max_iter = -1;
		int error_type = -1;
		double error_tolerance = -1;
		int monitoring_level = -1;

		// when no parameter file passed - take control parameters from problem input file
		if(0==strlen(ctrl_heat->solver_filename)){

			strcpy(solver_heat_filename, ctrl_heat->solver_filename);
			max_iter = pdr_lins_i_params(problem_id, 2); // max_iter
			error_type = pdr_lins_i_params(problem_id, 3); // error_type
			error_tolerance = pdr_lins_d_params(problem_id, 4); // error_tolerance
			monitoring_level = pdr_lins_i_params(problem_id, 5);  // monitoring level


		}
		else{

			sprintf(solver_heat_filename, "%s/%s",
					ctrl_heat->work_dir, ctrl_heat->solver_filename);

		}

		int parallel = SIC_SEQUENTIAL;
#ifdef PARALLEL
		parallel = SIC_PARALLEL;
#endif

		/*kbw*/
		fprintf(Interactive_output, "initializing heat solver (file %s)\n", solver_heat_filename);
		fprintf(Interactive_output, "parameters: parallel %d, maxiter %d, error_type %d, error_meas %.15lf, monitor %d\n",
				parallel, max_iter,  error_type,  error_tolerance, monitoring_level);
		/*kew*/

		solver_heat_id = sir_init(solver_type, parallel,
								  pdr_get_max_num_grid_levels(problem_id),
								  solver_heat_filename, max_iter,
								  error_type, error_tolerance, monitoring_level);
		ctrl_heat->solver_id = solver_heat_id;
		fprintf(Interactive_output, "\nAssigned solver IDs for heat %d\n",
				solver_heat_id);
	}

	/* set parameter indicating new mesh */
	iadapt = 1;

	/* if(ctrl_ns_supg->name == 110){ // for flux_radconv_isothermal problem  */
	/*   dtime_adapt = initial_dtime; */
	/*   time_ns_supg->cur_dtime = dtime_adapt; */
	/*   if ( dtime_adapt < final_dtime ) { */
	/*     time_adapt = log10(final_dtime / dtime_adapt); */
	/*   } else { */
	/*     time_adapt = -1.0; */
	/*   } */
	/*   fprintf(Interactive_output, */
	/* 	    "\n\tAS: Time adaptation (1) --> dtime_adapt = %1.12lf\tfinal_dtime = %1.12lf\ttime_adapt in %.0lf steps\n", dtime_adapt, final_dtime, time_adapt); */
	/* } */

	// #ifdef PHASE_TRANSFORMATION
	// 	fprintf(Interactive_output, "initializing phases field\n");
	//     pdr_phases_field_init(field_heat_id);
	// #endif

	/* start loop over time steps */
	while (utv_SIGINT_not_caught) {

		/* update time step and current time */
		++(time_heat->cur_step);

		// here possible time step length adaptations !!!

		// compute maximal, minimal and average CFL number for ns_supg_problem
		/*     if(time_heat->CFL_control == 1){ */

		/*       cfl_limit = time_heat->CFL_limit; */
		/*       time_step_length_mult = time_heat->time_step_length_mult; */

		/*       fprintf(Interactive_output, */
		/* 	      "\nCFL numbers before time step length control:\n"); */
		/*       pdr_heat_compute_CFL(pdv_heat_current_problem_id, &cfl_min, &cfl_max, &cfl_ave); */

		/* #ifdef PARALLEL */
		/*       double cfl_max_new = cfl_max; */
		/*       pcr_allreduce_max_double( 1, &cfl_max, &cfl_max_new); */
		/*       cfl_max = cfl_max_new;  */
		/* #endif */

		/*       fprintf(Interactive_output, */
		/* 	      "CFL_min = %lf, CFL_max = %lf, CFL_average = %lf\n", */
		/* 	      cfl_min,cfl_max,cfl_ave); */

		/*       // check whether CFL is not too big */
		/*       if(cfl_max > cfl_limit){ */
		/* 	time_heat->cur_dtime = time_heat->cur_dtime/(4*time_step_length_mult); */
		/*       } */
		/*       // check whether CFL is not too small */
		/*       if(cfl_max < cfl_limit/(time_step_length_mult) ){ */
		/* 	if(cfl_max > 1.e-3){ // to exclude pure heat diffusion problems */
		/* 	  time_heat->cur_dtime=time_heat->cur_dtime*time_step_length_mult; */
		/* 	} */
		/*       } */
		/*     } */

		fprintf(Interactive_output, "\n\n***************************************************************************\n");
		fprintf(Interactive_output, "\n\n***************************************************************************\n");
		fprintf(Interactive_output, "\n\nSolving time step %d (step_length: %lf, initial time: %lf)\n",
				time_heat->cur_step, time_heat->cur_dtime, time_heat->cur_time);

		time_heat->cur_time += time_heat->cur_dtime;

		// compute maximal, minimal and average CFL number for heat_problem
		/*     pdr_heat_compute_CFL(pdv_heat_current_problem_id, &cfl_min, &cfl_max, &cfl_ave); */
		/* #ifdef PARALLEL */
		/*     double cfl_max_new = cfl_max; */
		/*     pcr_allreduce_max_double( 1, &cfl_max, &cfl_max_new); //INSIDE ASSERT FAILS!!! */
		/*     cfl_max = cfl_max_new; */
		/* #endif */
		/*     fprintf(Interactive_output, */
		/* 	    "CFL_min = %lf, CFL_max = %lf, CFL_average = %lf\n", */
		/* 	    cfl_min,cfl_max,cfl_ave); */

		// rewrite solution:
		// current from previous time step becomes previous for current time step
		// soldofs_1 are the most recent solution dofs
		// at the beginning of time step they are rewritten to soldofs_3
		// that holds the solution from the previous time step
		utr_rewr_sol(field_heat_id, 1, 3);	//rewrite current -> u_n (3)
		//pdr_rewr_heat_dtdt_sol(field_heat_dtdt_id, 1, 2);	//rewrite current (1) to previous (2)
		/* update time dependent boundary conditions for NS problem*/
		/* update time dependent boundary conditions for thermal problem */
		pdr_heat_update_timedep_bc(&pdv_heat_problem.bc, time_heat->cur_time);

		nonl_iter = 0; double initial_sol_norm_uk_heat = 0.0;
		for(;;) {
			fprintf(Interactive_output, "\n---------------------------------------------------------------------------\n");
			fprintf(Interactive_output, "Solving nonlinear iteration %d\n", nonl_iter);

			// when forming linear system - soldofs_1 are equal to soldofs_2
			// after solving the system soldofs_1 are different than soldofs_2
			// before starting new solution soldofs_1 are rewritten to soldofs_2
			// simply: soldofs_1 are u_k+1, soldofs_2 are u_k and soldofs_3 are u_n
			utr_rewr_sol(field_heat_id,
						 Current_solution_ID, // 1,
						 Previous_iteration_sol_ID // 2
						 );	//rewrite current -> u_k (2)


			if (iadapt == 1) {
#ifdef PARALLEL
				/* initiate exchange tables for DOFs - for one field, max_num_levels */
				// one DOF - offset = 0, nreq = 1;
				int nr_levels = pdr_get_max_num_grid_levels(pdv_heat_current_problem_id);
				pdv_exchange_table_index = appr_init_exchange_tables(pcr_nr_proc(),
																	 pcr_my_proc_id(),
																	 field_heat_id, 0, 1,
																	 nr_levels);
#endif
				if(pdr_lins_i_params(pdv_heat_current_problem_id, 1) >= 0){ // mkb interface for solvers
					sir_create(solver_heat_id, pdv_heat_current_problem_id);
				}
				iadapt = 0;
			}


			if(pdr_lins_i_params(pdv_heat_current_problem_id, 1) >= 0){ // mkb interface for solvers

				int problem_id = pdv_heat_current_problem_id;
				int nr_iter = -1;
				double conv_meas = -1.0;
				double conv_rate = -1.0;
				int monitor = -1;

				// when no parameter file passed - take control parameters from problem input file
				if(0==strlen(ctrl_heat->solver_filename)){
					nr_iter = pdr_lins_i_params(problem_id, 2); // max_iter
					conv_meas = pdr_lins_d_params(problem_id, 4); // error_tolerance
					monitor = pdr_lins_i_params(problem_id, 5);  // monitoring level
				}

				int ini_guess = 1; //  get initial guess from data structure
				/*---------- CALLING DIRECT OR ITERATIVE SOLVER THROUGH MKB INTERFACE ---------*/
				sir_solve(solver_heat_id, SIC_SOLVE, ini_guess, monitor,
						  &nr_iter, &conv_meas, &conv_rate);
				if(pdr_lins_i_params(problem_id, 1) > 0){ // iterative solver
					fprintf(Interactive_output,
							"\nAfter %d iterations of linear solver for conv_diff problem\n",
							nr_iter);
					fprintf(Interactive_output,
							"Convergence measure: %lf, convergence rate %lf\n",
							conv_meas, conv_rate);
				}

			} else {

				sir_direct_solve_lin_sys(pdv_heat_current_problem_id, SIC_SEQUENTIAL,
										 solver_heat_filename);
			}


			/* post-process solution using slope limiter */
			//if(slope) {
			//  iaux=pdr_slope_limit(problem_id);
			// pdr_slope_limit for DG is archived in pdd_conv_diff/approx_dg/..._util.c
			//  if(iaux<0) { printf("\nError in slope!\n");getchar();}
			//}


			pdr_heat_sol_diff_norm(
						Current_solution_ID, // 1,
						Previous_iteration_sol_ID, // 2,
						&sol_norm_uk_heat);
			fprintf(Interactive_output, "\nAfter linear solver in nonlinear iteration %d\n",
					nonl_iter);
			fprintf(Interactive_output, "Solution difference norm (u_k, prev. u_k): %15.12f (limit %15.12f)\n",
					sol_norm_uk_heat, nonl_heat->conv_meas);

			if(nonl_iter==0) initial_sol_norm_uk_heat = sol_norm_uk_heat;

			++nonl_iter;

			int break_indi=0;
			if( nonl_heat->conv_type==0 ){

				if( sol_norm_uk_heat < nonl_heat->conv_meas ){

					fprintf(Interactive_output,
							"\nConvergence in non-linear iterations: %15.12f < %15.12f - breaking.\n",
							sol_norm_uk_heat, nonl_heat->conv_meas);

					//break;
					break_indi = 1;
				}
			}
			else if( nonl_heat->conv_type==1 ){

				if( sol_norm_uk_heat < nonl_heat->conv_meas*initial_sol_norm_uk_heat ){
					fprintf(Interactive_output,
							"\nConvergence in non-linear iterations: %15.12f < %15.12f * %15.12f - breaking.\n",
							sol_norm_uk_heat, nonl_heat->conv_meas, initial_sol_norm_uk_heat);
					//break;
					break_indi = 1;
				}
			}

			if (nonl_iter >= nonl_heat->max_iter) {
				fprintf(Interactive_output,
						"\nMax nonlinear iterations (%d) reached - breaking.\n",
						nonl_heat->max_iter);
				//break;
				break_indi = 1;
			}

			// check whether we should change time step length
			if(break_indi == 1){

				if(time_heat->time_step_length_nonl_control == 1){

					fprintf(Interactive_output,
							"In time step control (iter_min %d <=? nonl_iter %d <=? iter_max %d)\n",
							time_heat->time_step_length_nonl_iter_min, nonl_iter,
							time_heat->time_step_length_nonl_iter_max );

					/* time step length adaptation based on number non-linear iterations for convergence */
					int time_step_length_control = time_heat->time_step_length_nonl_control;
					double time_step_length_limit_max = time_heat->time_step_length_nonl_iter_max;
					double time_step_length_limit_min = time_heat->time_step_length_nonl_iter_min;
					double time_step_length_increase_mult =
							time_heat->time_step_length_nonl_iter_increase_mult;
					double time_step_length_decrease_mult =
							time_heat->time_step_length_nonl_iter_decrease_mult;

					// if convergence was slow - we decrease time step length
					if( (nonl_iter>=time_step_length_limit_max)
							// or if norm is above required value
							||  (sol_norm_un_heat > time_heat->time_step_length_max_norm) ) {
						time_heat->cur_dtime = time_heat->cur_dtime*time_step_length_decrease_mult;
					}
					// if convergence is fast we increase time step length
					else if( (nonl_iter<=time_step_length_limit_min)
							 // or if norm is below required value
							 || (sol_norm_un_heat < time_heat->time_step_length_min_norm) ){
						time_heat->cur_dtime = time_heat->cur_dtime*time_step_length_increase_mult;
					}

					if(time_heat->cur_dtime < time_heat->time_step_length_min) {
						time_heat->cur_dtime = time_heat->time_step_length_min;
						mf_log_info("Current time step(%lf) is below minimum, reseting to time_step_length_min(%lf)",
									time_heat->cur_dtime, time_heat->time_step_length_min);
					}
					else if(time_heat->cur_dtime > time_heat->time_step_length_max) {
						time_heat->cur_dtime = time_heat->time_step_length_max;
						mf_log_info("Current time step(%lf) exceeds maximum, reseting to time_step_length_max(%lf)",
									time_heat->cur_dtime, time_heat->time_step_length_max);
					}
				}


				break;
			}



		} // end of infinite loop over non-linear iterations


		pdr_heat_sol_diff_norm(1, 3, &sol_norm_un_heat);
		fprintf(Interactive_output, "\n---------------------------------------------------------------------------\n");
		fprintf(Interactive_output, "After non-linear solver in time step %d\n",
				time_heat->cur_step);
		fprintf(Interactive_output, "Solution difference norm (u_n, prev. u_n): %15.12f (limit %15.12f)\n",
				sol_norm_un_heat,  time_heat->conv_meas);

		// AS: heating-cooling problem dofs shift 2->3, 1->2
		//pdr_rewr_heat_dtdt_sol(field_heat_dtdt_id, 2, 3);	//rewrite previous (2) to early (3)
		//pdr_rewr_heat_dtdt_sol(field_heat_dtdt_id, 1, 2);	//rewrite current (1) to previous (2)
		//pdr_heat_dtdt(1,3);  // compute and write Current=1 heat_dtdt field

#ifdef PHASE_TRANSFORMATION
		pdr_phases_transformation(field_heat_id);
#endif

		/* graphics data dumps */
		if (time_heat->intv_graph > 0 &&
				time_heat->cur_step % time_heat->intv_graph == 0) {

			fprintf(Interactive_output, "(Writing field to disk (graphics)...)\n");

			pdr_heat_write_paraview(Work_dir,
									Interactive_input, Interactive_output);
		}

		/* full data dumps (for restarting) */
		if (time_heat->intv_dumpout > 0
				&& time_heat->cur_step % time_heat->intv_dumpout == 0) {

			fprintf(Interactive_output, "\n(Writing field to disk - not implemented)\n");

			//pdr_heat_dump_data(Work_dir,
			//			   Interactive_input, Interactive_output);

		}


		/* check for stop conditions */

		/* stop when final time reached */
		if (time_heat->cur_time >= time_heat->final_time) {
			fprintf(Interactive_output, "\nFinal time reached. Stopping.\n");
			break;
		}

		/* cur_dtime could change through simulation (time step adaptation) */
		/* thus we need to check also if cur_step not bigger than final_step */
		if (time_heat->cur_step >= time_heat->final_step) {
			fprintf(Interactive_output, "\nFinal step reached. Stopping.\n");
			break;
		}

		/* stop if convergence reached (for stationary problems) */
		if ( sol_norm_un_heat < time_heat->conv_meas) {
			fprintf(Interactive_output,
					"\nConvergence in time integration!\n");
			break;
		}

		// when time for adaptation
		if( (adpt_heat->type>0 && adpt_heat->interval>0 &&
			 (time_heat->cur_step+1)%adpt_heat->interval==0)){

#ifdef PARALLEL
			/* free exchange tables for DOFs */
			appr_free_exchange_tables(pdv_exchange_table_index);
			pdv_exchange_table_index = -1;
#endif

			if(pdr_lins_i_params(pdv_heat_current_problem_id, 1) >= 0){ // mkb interface for solvers
				sir_free(solver_heat_id);
			}


			pdr_heat_adapt(Work_dir, Interactive_input, Interactive_output);
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


		if(pdr_lins_i_params(pdv_heat_current_problem_id, 1) >= 0){ // mkb interface for solvers
			sir_free(solver_heat_id);
		}


	}

	if(pdr_lins_i_params(pdv_heat_current_problem_id, 1) >= 0){ // mkb interface for solvers
		sir_destroy(solver_heat_id);
	}

	fprintf(Interactive_output, "\n**************** LEAVING TIME INTEGRATION ***************\n");

	return;
}

