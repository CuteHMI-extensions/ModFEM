/************************************************************************
File pds_ns_supg_problem_io.c - problem data read

Contains definition of routines:
  pdr_ns_supg_problem_read - read problem data

------------------------------
History:
	2011    - Przemyslaw Plaszewski (pplaszew@agh.edu.pl)
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libconfig.h>

/* problem dependent module interface */
#include <modfem/pdh_intf.h>

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>

/* problem module's types and functions */
#include <modfem/pd_ns_supg/pdh_ns_supg.h>
/* problem module's types and functions */
#include <modfem/pd_ns_supg/pdh_ns_supg_problem.h> /* USES & IMPLEMENTS */
/* include approximation module header */
#include <modfem/aph_intf.h>

#include <modfem/uth_log.h>

/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */

/*------------------------------------------------------------
pdr_ns_supg_problem_read - read problem data
------------------------------------------------------------*/
int pdr_ns_supg_problem_read(
  char *Work_dir,
  char *Filename,
  FILE *Interactive_output,
  pdt_ns_supg_problem *Problem,
  int Nr_sol // nr_sol is time integration dependent - specified by main
)
{
  // ??????????????????????????????????????
  config_t cfg;
  config_setting_t *root_setting;
  //const char *str;
  char problem_file[300];

  Problem->ctrl.nreq = PDC_NS_SUPG_NREQ;
  Problem->ctrl.nr_sol = Nr_sol;

  sprintf(problem_file, "%s/%s", Work_dir, Filename);
  fprintf(Interactive_output, "\nOpening NS_SUPG problem file %s\n",
	  problem_file);

  config_init(&cfg);

  /* Read the file. If there is an error, report it and exit. */
  if(! config_read_file(&cfg, problem_file))
  {
//    fprintf(Interactive_output, "Problem file error: %s:%d - %s\n",
//	    config_error_file(&cfg),
//            config_error_line(&cfg), config_error_text(&cfg));
	config_destroy(&cfg);
	return(EXIT_FAILURE);
  }

  root_setting = config_lookup(&cfg, "problem");
  if(root_setting != NULL)
  {
	int count = config_setting_length(root_setting);

	int i, j;
	int setting_length;
	//const char *name;

	if(count <= 0)
	{
	  fprintf(Interactive_output, "Problem file error: Problem not defined\n");
	  config_destroy(&cfg);
	  return(EXIT_FAILURE);
	}

	if(count > 1)
	{
	  fprintf(Interactive_output,
		  "Problem file error: Currently only one problem supported\n");
	  config_destroy(&cfg);
	  return(EXIT_FAILURE);
	}


	for(i = 0; i < count; ++i)
	{
	  config_setting_t *problem_sett = config_setting_get_elem(root_setting, i);
	  mf_check_mem(problem_sett);

	  config_setting_t *setting=NULL;
	  const char* str=NULL;

/* CONTROL PARAMETERS - GENERIC: FIELDS REQUIRED FOR ALL PROBLEMS */

	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "name", &str) ) {
		  mf_fatal_err("Problem file error at 'name' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.name, str);

	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "mesh_type", &str)) {
		  mf_fatal_err("Problem file error at 'mesh_type' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.mesh_type, str);


	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "mesh_file_in", &str)) {
		  mf_fatal_err("Problem file error at 'mesh_file_in' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.mesh_filename, str);

	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "field_file_in", &str)) {
		  mf_fatal_err("Problem file error at 'field_file_in' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.field_filename, str);

	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "materials_file", &str)) {
		  mf_fatal_err("Problem file error at 'materials_file' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.material_filename, str);

	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "bc_file", &str)) {
		  mf_fatal_err("Problem file error at 'bc_file' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.bc_filename, str);

	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "solver_file", &str)) {
		  mf_fatal_err("Problem file error at 'solver_file' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.solver_filename, str);


	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "mesh_file_out", &str)) {
		  mf_fatal_err("Problem file error at 'mesh_file_out' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.mesh_dmp_filepattern, str);

	  if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "field_file_out", &str)) {
		  mf_fatal_err("Problem file error at 'field_file_out' string.");
	  }
	  mf_check_mem(str);
	  strcpy(Problem->ctrl.field_dmp_filepattern, str);

	  // get PDEG value
	  setting = config_setting_get_member(problem_sett, "pdeg");
	  if(setting != NULL) {
	Problem->ctrl.pdeg = (int)config_setting_get_int(setting);
	  } else {
	Problem->ctrl.pdeg = APC_P1_NONE_G1;
	  }

/* CONTROL PARAMETERS - SPECIFIC TO NS_SUPG PROBLEMS */

	  setting = config_setting_get_member(problem_sett, "penalty");
	  Problem->ctrl.penalty = (double)config_setting_get_float(setting);

	  for(j=0; j<3; ++j)
	  {
		setting = config_setting_get_member(problem_sett, "gravity");
	Problem->ctrl.gravity_field[j]=config_setting_get_float_elem(setting, j);
	  }

	  // if material file not specified material data are read from problem file
	  if(strlen(Problem->ctrl.material_filename)==0){


	// indicate material data are constant
	Problem->ctrl.ref_temperature = -1.0;

	setting = config_setting_get_member(problem_sett, "density");
	Problem->ctrl.density =
								 (double)config_setting_get_float(setting);


	setting = config_setting_get_member(problem_sett, "dynamic_viscosity");
	Problem->ctrl.dynamic_viscosity =
									 (double)config_setting_get_float(setting);

	  }
	  else{ // if material file specified


	// used to distinguish between constant and temperature dependent material data
	// if non-constant data reference_temperature MUST be specified (it is used to compute
	// reference_viscosity and reference_density as data at reference_tmperature)
	// otherwise it is automatically set to -1.0 to indicate constant data
	setting = config_setting_get_member(problem_sett, "reference_temperature");
	Problem->ctrl.ref_temperature = (double)config_setting_get_float(setting);

	Problem->ctrl.density = -1.0;
	Problem->ctrl.dynamic_viscosity = -1.0;

	  }

	  // used only to compute Reynolds number
	  setting = config_setting_get_member(problem_sett, "reference_length");
	  Problem->ctrl.ref_length =
								 (double)config_setting_get_float(setting);

	  // used only to compute Reynolds number
	  setting = config_setting_get_member(problem_sett, "reference_velocity");
	  Problem->ctrl.ref_velocity =
								 (double)config_setting_get_float(setting);

	  // OBSOLETE - Reynolds number always computed
	  /* Problem->ctrl.reynolds_number = Problem->ctrl.density *  */
	  /* 	Problem->ctrl.ref_velocity * Problem->ctrl.ref_length / */
	  /* 	Problem->ctrl.dynamic_viscosity; */


	  pdv_ns_supg_problem.ctrl.error_indicator = NO_ERROR;

	  setting = config_setting_get_member(problem_sett, "final_velocity_mag");
	  if(setting != NULL) {
		   Problem->ctrl.final_velocity_mag = (int)config_setting_get_float(setting);
	  }
	  else {
		  Problem->ctrl.final_velocity_mag = 0.0;
	  }


/* TIME INTEGRATION PARAMETERS - GENERIC: FOR ALL TYPES OF PROBLEMS */

	  setting = config_setting_get_member(problem_sett, "time_integration_type");
	  Problem->time.type = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett,"implicitness_parameter");
	  Problem->time.alpha =
								 (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "current_timestep");
	  Problem->time.cur_step = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "current_time");
	  Problem->time.cur_time = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "current_timestep_length");
	  Problem->time.cur_dtime = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "previous_timestep_length");
	  Problem->time.prev_dtime = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "final_timestep");
	  Problem->time.final_step = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "final_time");
	  Problem->time.final_time = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "time_error_type");
	  Problem->time.conv_type = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "time_error_tolerance");
	  Problem->time.conv_meas = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "time_monitoring_level");
	  Problem->time.monitor = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "full_dump_intv");
	  Problem->time.intv_dumpout = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "graph_dump_intv");
	  Problem->time.intv_graph = (int)config_setting_get_int(setting);

/* TIME INTEGRATION PARAMETERS - SPECIFIC FOR FLOW PROBLEMS */
	  setting = config_setting_get_member(problem_sett, "CFL_control");
	  if(setting == NULL) { //added by KMichalik 17.09.2012
		fprintf(Interactive_output,
				"Warning: CFL_control not found in problem file. Assuming defaults.\n");
		Problem->time.CFL_control = 0;
	  }
	  else {
		Problem->time.CFL_control = (int)config_setting_get_int(setting);
	  }

	  if(Problem->time.CFL_control == 1){
	setting = config_setting_get_member(problem_sett, "CFL_limit_max");
	Problem->time.CFL_limit_max = (double)config_setting_get_float(setting);
	setting = config_setting_get_member(problem_sett, "CFL_limit_min");
	Problem->time.CFL_limit_min = (double)config_setting_get_float(setting);
	setting = config_setting_get_member(problem_sett, "CFL_time_step_length_increase_mult");
	Problem->time.CFL_time_step_length_increase_mult = (double)config_setting_get_float(setting);
	setting = config_setting_get_member(problem_sett, "CFL_time_step_length_decrease_mult");
	Problem->time.CFL_time_step_length_decrease_mult = (double)config_setting_get_float(setting);
	  }

/* TIME STEP LENGTH CONTROL PARAMETERS BASED ON CONVERGENCE OF NONLINEAR ITERATIONS - GENERIC: FOR ALL PROBLEMS */
	  setting = config_setting_get_member(problem_sett, "time_step_length_nonl_control");
	  if(setting == NULL) {
	fprintf(Interactive_output,
		"Warning: time_step_length_nonl_control not found in problem file. Assuming defaults.\n");
	Problem->time.time_step_length_nonl_control = 0;
	  }
	  else {
	Problem->time.time_step_length_nonl_control = (int)config_setting_get_int(setting);
	  }

	  if(Problem->time.time_step_length_nonl_control == 1){
	setting = config_setting_get_member(problem_sett, "time_step_length_nonl_iter_max");
	Problem->time.time_step_length_nonl_iter_max = (int)config_setting_get_int(setting);
	setting = config_setting_get_member(problem_sett, "time_step_length_nonl_iter_min");
	Problem->time.time_step_length_nonl_iter_min = (int)config_setting_get_int(setting);
	setting = config_setting_get_member(problem_sett, "time_step_length_nonl_iter_increase_mult");
	Problem->time.time_step_length_nonl_iter_increase_mult = (double)config_setting_get_float(setting);
	setting = config_setting_get_member(problem_sett, "time_step_length_nonl_iter_decrease_mult");
	Problem->time.time_step_length_nonl_iter_decrease_mult = (double)config_setting_get_float(setting);

	// Min and Max allowed values for time step value.
	setting = config_setting_get_member(problem_sett, "time_step_length_min");
	if(setting != NULL) {
		Problem->time.time_step_length_min = (double)config_setting_get_float(setting);
	}
	else {
	  mf_log_err("time_step_length_min not found - assuming 0.0");
	  Problem->time.time_step_length_min = 0.0;
	}

	setting = config_setting_get_member(problem_sett, "time_step_length_max");
	if(setting != NULL) {
		Problem->time.time_step_length_max = (double)config_setting_get_float(setting);
	}
	else {
	  mf_log_err("time_step_length_max not found - assuming 1.0e10");
	  Problem->time.time_step_length_max = 1.0e10;
	}


	// MIN and MAX values for error Norm forcing time step adaptation, regardless form other conditions.
	setting = config_setting_get_member(problem_sett, "time_step_length_min_norm");
	if(setting != NULL) {
		Problem->time.time_step_length_min_norm = (double)config_setting_get_float(setting);
	}
	else {
		mf_log_err("time_step_length_min_norm not found - assuming %lf ",0.0 );
		Problem->time.time_step_length_min_norm = 0.0 ;
	}

	setting = config_setting_get_member(problem_sett, "time_step_length_max_norm");
	if(setting != NULL) {
		Problem->time.time_step_length_max_norm = (double)config_setting_get_float(setting);
	}
	else {
		mf_log_err("time_step_length_max_norm not found - assuming %lf",10e10);
		Problem->time.time_step_length_max_norm = 10e10;
	}
	  }

/* NONLINEAR SOLVER PARAMETERS - GENERIC: FOR ALL PROBLEMS */

	  setting = config_setting_get_member(problem_sett, "nonlinear_solver_type");
	  Problem->nonl.type = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "max_nonl_iter");
	  Problem->nonl.max_iter = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "nonl_error_type");
	  Problem->nonl.conv_type = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "nonl_error_tolerance");
	  Problem->nonl.conv_meas = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "nonl_monitoring_level");
	  Problem->nonl.monitor = (int)config_setting_get_int(setting);

/* LINEAR SOLVER PARAMETERS - GENERIC: FOR ALL PROBLEMS */

	  setting = config_setting_get_member(problem_sett, "linear_solver_type");
	  Problem->lins.type = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "max_linear_solver_iter");
	  Problem->lins.max_iter = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "linear_solver_error_type");
	  Problem->lins.conv_type = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "linear_solver_error_tolerance");
	  Problem->lins.conv_meas = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "lin_solv_monitoring_level");
	  Problem->lins.monitor = (int)config_setting_get_int(setting);

/* ADAPTATION PARAMETERS - GENERIC: FOR ALL PROBLEMS */

	  setting = config_setting_get_member(problem_sett, "adapt_type");
	  Problem->adpt.type = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "adapt_interval");
	  Problem->adpt.interval = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "adapt_maxgen");
	  Problem->adpt.maxgen = (int)config_setting_get_int(setting);

	  setting = config_setting_get_member(problem_sett, "adapt_tolerance");
	  Problem->adpt.eps = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett, "adapt_deref_ratio");
	  Problem->adpt.ratio = (double)config_setting_get_float(setting);

	  setting = config_setting_get_member(problem_sett,"adapt_monitoring_level");
	  Problem->adpt.monitor = (int)config_setting_get_int(setting);


	}

  }

  config_destroy(&cfg);

  return(EXIT_SUCCESS);
}

/*------------------------------------------------------------
pdr_ns_supg_problem_clear - clear problem data
------------------------------------------------------------*/
int pdr_ns_supg_problem_clear(pdt_ns_supg_problem *Problem)
{
  //Problem->ctrl.name = -1;
  Problem->ctrl.mesh_id = -1;
  Problem->ctrl.field_id = -1;
  Problem->ctrl.nr_sol = -1;
  Problem->ctrl.nreq = -1;
  Problem->ctrl.solver_id = -1;
  Problem->ctrl.pdeg = -1;

  // ns_supg specific
  Problem->ctrl.penalty = -1.0;
  Problem->ctrl.gravity_field[0] = -1.0;
  Problem->ctrl.gravity_field[1] = -1.0;
  Problem->ctrl.gravity_field[2] = -1.0;
  Problem->ctrl.ref_temperature = -1.0;
  Problem->ctrl.density = -1.0;
  Problem->ctrl.dynamic_viscosity = -1.0;
  Problem->ctrl.ref_length = -1.0;
  Problem->ctrl.ref_velocity = -1.0;
  Problem->ctrl.error_indicator = NO_ERROR;
  Problem->ctrl.final_velocity_mag = 0.0;

  Problem->time.type = -1;
  Problem->time.alpha = -1.0;
  Problem->time.cur_step = -1;
  Problem->time.final_step = -1;
  Problem->time.cur_time = -1.0;
  Problem->time.final_time = -1.0;
  Problem->time.cur_dtime = -1.0;
  Problem->time.prev_dtime = -1.0;
  Problem->time.conv_type = -1;
  Problem->time.conv_meas = -1.0;
  Problem->time.monitor = -1;
  Problem->time.intv_dumpout = -1;
  Problem->time.intv_graph = -1;
  Problem->time.graph_accu = -1;

  // ns_supg specific
  Problem->time.CFL_control = -1;
  Problem->time.CFL_limit_max = -1.0;
  Problem->time.CFL_limit_min = -1.0;
  Problem->time.CFL_time_step_length_decrease_mult = -1.0;
  Problem->time.CFL_time_step_length_increase_mult = -1.0;

  Problem->time.time_step_length_nonl_control = -1;
  Problem->time.time_step_length_nonl_iter_max = -1.0;
  Problem->time.time_step_length_nonl_iter_min = -1.0;
  Problem->time.time_step_length_nonl_iter_decrease_mult = -1.0;
  Problem->time.time_step_length_nonl_iter_increase_mult = -1.0;
  Problem->time.time_step_length_min = 0.0;
  Problem->time.time_step_length_max = 10e10;
  Problem->time.time_step_length_max_norm = 10e10;
  Problem->time.time_step_length_min_norm = 0.0;

  Problem->nonl.type = -1;
  Problem->nonl.max_iter = -1;
  Problem->nonl.conv_type = -1;
  Problem->nonl.conv_meas = -1.0;
  Problem->nonl.monitor = -1;

  Problem->lins.type = -1;
  Problem->lins.max_iter = -1;
  Problem->lins.conv_type = -1;
  Problem->lins.conv_meas = -1.0;
  Problem->lins.monitor = -1;

  Problem->adpt.type = -1;
  Problem->adpt.interval = -1;
  Problem->adpt.maxgen = -1;
  Problem->adpt.eps = -1.0;
  Problem->adpt.ratio = -1.0;
  Problem->adpt.monitor = -1;

  return 0;
}
