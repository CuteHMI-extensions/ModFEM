/************************************************************************
File pds_heat_problem_io.c - problem data read

Contains definition of routines:
  pdr_heat_problem_read - read problem data

------------------------------
History:
	2011    - Przemyslaw Plaszewski (pplaszew@agh.edu.pl)
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libconfig.h>
#include <float.h>

/* problem module's types and functions */
#include <modfem/pd_heat/pdh_heat.h>		/* USES */
/* problem module's types and functions */
#include <modfem/pd_heat/pdh_heat_problem.h> /* USES & IMPLEMENTS */
/* types and functions related to materials handling */
#include <modfem/pd_heat/pdh_heat_materials.h>
/* include approximation module header */
#include <modfem/aph_intf.h>

#include <modfem/uth_log.h>

/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_heat_ */
/* - argument names start uppercase */

/*------------------------------------------------------------
pdr_heat_problem_read - read problem data
------------------------------------------------------------*/
int pdr_heat_problem_read(
		char * Work_dir,
		char * Filename,
		FILE * Interactive_output,
		pdt_heat_problem * Problem,
		int Nr_sol // nr_sol is time integration dependent - specified by main
)
{
	config_t cfg;
	config_setting_t * root_setting;
	//const char *str;
	char problem_file[300];

	Problem->ctrl.nreq = PDC_HEAT_NREQ;
	Problem->ctrl.nr_sol = Nr_sol;

#ifdef PHASE_TRANSFORMATION
	Problem->ctrl.phase_transformation_nreq = PDC_HEAT_PT_MAXEQ;
	Problem->ctrl.phase_transformation_nr_sol = Nr_sol;
	Problem->ctrl.phases_nreq = PDC_HEAT_P_MAXEQ;
	Problem->ctrl.phases_nr_sol = Nr_sol;
#endif

	sprintf(problem_file, "%s/%s", Work_dir, Filename);
	fprintf(Interactive_output, "\nOpening HEAT problem file %s\n",
			problem_file);

	config_init(&cfg);

	/* Read the file. If there is an error, report it and exit. */
	if (! config_read_file(&cfg, problem_file))
	{
		printf("Error reading heat problem file\n");
		exit(-1);

//    fprintf(Interactive_output, "Problem file error: %s:%d - %s\n", config_error_file(&cfg),
//            config_error_line(&cfg), config_error_text(&cfg));
//    config_destroy(&cfg);
		return (EXIT_FAILURE);
	}

	root_setting = config_lookup(&cfg, "problem");
	if (root_setting != NULL)
	{
		int count = config_setting_length(root_setting);
		int i, j;
		int setting_length;
		//const char *name;

		if (count == 0)
		{
			fprintf(Interactive_output, "Problem file error: Problem not defined\n");
			config_destroy(&cfg);
			return (EXIT_FAILURE);
		}

		if (count > 1)
		{
			fprintf(Interactive_output, "Problem file error: Currently only one problem supported\n");
			config_destroy(&cfg);
			return (EXIT_FAILURE);
		}


		for (i = 0; i < count; ++i)
		{
			config_setting_t * problem_sett = config_setting_get_elem(root_setting, i);
			config_setting_t * setting;
			const char * str;

			/* CONTROL PARAMETERS - GENERIC: FIELDS REQUIRED FOR ALL PROBLEMS */

			if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "name", &str)) {
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

			// get PDEG value
			setting = config_setting_get_member(problem_sett, "pdeg");
			if (setting != NULL) {
				Problem->ctrl.pdeg = (int)config_setting_get_int(setting);
			} else {
				Problem->ctrl.pdeg = APC_P1_NONE_G1;
			}

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

#ifdef PHASE_TRANSFORMATION
			if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "phase_transformation_field_file_in", &str)) {
				mf_fatal_err("Problem file error at 'phase_transformation_field_file_in' string.");
			}
			mf_check_mem(str);
			strcpy(Problem->ctrl.phase_transformation_field_filename, str);
			if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "phase_transformation_field_file_out", &str)) {
				mf_fatal_err("Problem file error at 'phase_transformation_field_file_out' string.");
			}
			mf_check_mem(str);
			strcpy(Problem->ctrl.phase_transformation_field_dmp_filepattern, str);

			if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "phases_field_file_in", &str)) {
				mf_fatal_err("Problem file error at 'phases_field_file_in' string.");
			}
			mf_check_mem(str);
			strcpy(Problem->ctrl.phases_field_filename, str);
			if (CONFIG_FALSE == config_setting_lookup_string(problem_sett, "phases_field_file_out", &str)) {
				mf_fatal_err("Problem file error at 'phases_field_file_out' string.");
			}
			mf_check_mem(str);
			strcpy(Problem->ctrl.phases_field_dmp_filepattern, str);

			setting = config_setting_get_member(problem_sett, "austenitizing_temperature");
			Problem->ctrl.austenitizing_temperature = (double)config_setting_get_float(setting);

#endif

			/* CONTROL PARAMETERS - SPECIFIC TO HEAT PROBLEMS */

			setting = config_setting_get_member(problem_sett, "penalty");
			Problem->ctrl.penalty = (double)config_setting_get_float(setting);


			// for boundary conditions and initial condition as well
			setting = config_setting_get_member(problem_sett, "ambient_temperature");
			Problem->ctrl.ambient_temperature =
					(double)config_setting_get_float(setting);


			// if material file not specified material data are read from problem file
			if (strlen(Problem->ctrl.material_filename) == 0) {

				// indicate material data are constant
				Problem->ctrl.ref_temperature = -1.0;

				setting = config_setting_get_member(problem_sett, "thermal_conductivity");
				Problem->ctrl.thermal_conductivity = (double)config_setting_get_float(setting);

				setting = config_setting_get_member(problem_sett, "density");
				Problem->ctrl.density = (double)config_setting_get_float(setting);

				setting = config_setting_get_member(problem_sett, "specific_heat");
				Problem->ctrl.specific_heat = (double)config_setting_get_float(setting);

			}
			else { // if material file specified


				// used to distinguish between constant and temperature dependent material data
				// if non-constant data reference_temperature MUST be specified (it is used to compute
				// reference_viscosity and reference_density as data at reference_tmperature)
				// otherwise it is automatically set to -1.0 to indicate constant data
				setting = config_setting_get_member(problem_sett, "reference_temperature");
				Problem->ctrl.ref_temperature = (double)config_setting_get_float(setting);

				Problem->ctrl.thermal_conductivity = -1.0;
				Problem->ctrl.density = -1.0;
				Problem->ctrl.specific_heat = -1.0;

			}


			/* TIME INTEGRATION PARAMETERS - GENERIC: FOR ALL TYPES OF PROBLEMS */

			setting = config_setting_get_member(problem_sett, "time_integration_type");
			Problem->time.type = (int)config_setting_get_int(setting);

			setting = config_setting_get_member(problem_sett, "implicitness_parameter");
			Problem->time.alpha = (double)config_setting_get_float(setting);

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


			/* TIME STEP LENGTH CONTROL PARAMETERS BASED ON CONVERGENCE OF NONLINEAR ITERATIONS - GENERIC: FOR ALL PROBLEMS */
			setting = config_setting_get_member(problem_sett, "time_step_length_nonl_control");
			if (setting == NULL) {
				fprintf(Interactive_output,
						"Warning: time_step_length_nonl_control not found in problem file. Assuming defaults.\n");
				Problem->time.time_step_length_nonl_control = 0;
			}
			else {
				Problem->time.time_step_length_nonl_control = (int)config_setting_get_int(setting);
			}

			if (Problem->time.time_step_length_nonl_control == 1) {
				setting = config_setting_get_member(problem_sett, "time_step_length_nonl_iter_max");
				Problem->time.time_step_length_nonl_iter_max = (int)config_setting_get_int(setting);
				setting = config_setting_get_member(problem_sett, "time_step_length_nonl_iter_min");
				Problem->time.time_step_length_nonl_iter_min = (int)config_setting_get_int(setting);
				setting = config_setting_get_member(problem_sett, "time_step_length_nonl_iter_increase_mult");
				Problem->time.time_step_length_nonl_iter_increase_mult = (double)config_setting_get_float(setting);
				setting = config_setting_get_member(problem_sett, "time_step_length_nonl_iter_decrease_mult");
				Problem->time.time_step_length_nonl_iter_decrease_mult = (double)config_setting_get_float(setting);

				setting = config_setting_get_member(problem_sett, "time_step_length_min");
				if (setting != NULL) {
					Problem->time.time_step_length_min = (double)config_setting_get_float(setting);
				}
				else {
					mf_log_err("time_step_length_min not found - assuming 0.0");
					Problem->time.time_step_length_min = 0.0;
				}


				setting = config_setting_get_member(problem_sett, "time_step_length_max");
				if (setting != NULL) {
					Problem->time.time_step_length_max = (double)config_setting_get_float(setting);
				}
				else {
					mf_log_err("time_step_length_max not found - assuming 1.0e10");
					Problem->time.time_step_length_max = 1.0e10;
				}

				// MIN and MAX values for error Norm forcing time step adaptation, regardless form other conditions.
				setting = config_setting_get_member(problem_sett, "time_step_length_min_norm");
				if (setting != NULL) {
					Problem->time.time_step_length_min_norm = (double)config_setting_get_float(setting);
				}
				else {
					mf_log_err("time_step_length_min_norm not found - assuming  %lf ", 0.0 );
					Problem->time.time_step_length_min_norm = 0.0 ;
				}

				setting = config_setting_get_member(problem_sett, "time_step_length_max_norm");
				if (setting != NULL) {
					Problem->time.time_step_length_max_norm = (double)config_setting_get_float(setting);
				}
				else {
					mf_log_err("time_step_length_max_norm not found - assuming  %lf", 10e10);
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

			config_setting_lookup_string(problem_sett, "solver_file", &str);
			strcpy(Problem->ctrl.solver_filename, str);

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

			setting = config_setting_get_member(problem_sett, "adapt_monitoring_level");
			Problem->adpt.monitor = (int)config_setting_get_int(setting);


		}

	}

	config_destroy(&cfg);


	return (EXIT_SUCCESS);
}

/*------------------------------------------------------------
pdr_heat_problem_clear - clear problem data
------------------------------------------------------------*/
int pdr_heat_problem_clear(pdt_heat_problem * Problem)
{

	Problem->ctrl.mesh_id = -1;
	Problem->ctrl.field_id = -1;
	Problem->ctrl.nr_sol = -1;
	Problem->ctrl.nreq = -1;
	Problem->ctrl.solver_id = -1;
	Problem->ctrl.pdeg = -1;

	// heat specific
	Problem->ctrl.penalty = -1.0;
	Problem->ctrl.ref_temperature = -1.0;
	Problem->ctrl.ambient_temperature = -1.0;
	Problem->ctrl.thermal_conductivity = -1.0;
	Problem->ctrl.density = -1.0;
	Problem->ctrl.specific_heat = -1.0;

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

