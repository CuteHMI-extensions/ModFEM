#include "../../../include/modfem/nssupgheat/Problem.hpp"

#include <QSettings>
#include <QDir>
#include <QDataStream>
#include <QFuture>
#include <QtConcurrent>

#include <cstdlib>
#include <iostream>

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>

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
#include <modfem/mod_fem_viewer.h>	/* USES */

/* problem module's types and functions */
#include <modfem/pd_ns_supg_heat/pdh_ns_supg_heat.h>		/* USES */
/* types and functions related to problem structures */
#include <modfem/pd_ns_supg/pdh_ns_supg_problem.h>
#include <modfem/pd_heat/pdh_heat_problem.h>
// bc and material header files are included in problem header files
#include <modfem/pd_heat/pdh_heat_bc.h>

/***************************************/
/* DECLARATIONS OF INTERNAL PROCEDURES */
/***************************************/
/* Rules:
/* - name always begins with pdr_ns_supg_heat */
/* - argument names start uppercase */

/*------------------------------------------------------------
pdr_ns_supg_heat_post_process - simple post-processing
------------------------------------------------------------*/
double pdr_ns_supg_heat_post_process(
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output
);

/*------------------------------------------------------------
pdr_ns_supg_heat_profile - to dump a set of values along a line
------------------------------------------------------------*/
int pdr_ns_supg_heat_write_profile(
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output
);

/*------------------------------------------------------------
pdr_heat_initial_condition - procedure passed as argument
  to field initialization routine in order to provide problem
  dependent initial condition data
------------------------------------------------------------*/
double pdr_heat_initial_condition(
		int Field_id, // field_id - each problem should know its field id
		double * Coor,  // point coordinates
		int Sol_comp_id // solution component
);

/*------------------------------------------------------------
pdr_ns_supg_initial_condition - procedure passed as argument
  to field initialization routine in order to provide problem
  dependent initial condition data
------------------------------------------------------------*/
double pdr_ns_supg_initial_condition(
		int Field_id, // field_id - each problem should know its field id
		double * Coor,  // point coordinates
		int Sol_comp_id // solution component
);

/*------------------------------------------------------------
 TODO: not implemented yet
------------------------------------------------------------*/
int  pdr_ns_supg_heat_bc_free();

/*------------------------------------------------------------
 TODO: not implemented yet
------------------------------------------------------------*/
int  pdr_ns_supg_heat_material_free();

/***************************************/
/* DEFINITIONS OF PROCEDURES */
/***************************************/

namespace modfem {
namespace nssupgheat {

Problem::Problem(QObject * parent):
	QObject(parent),
	m(new Members{".",
					nullptr,
					nullptr,
					0,
					0,
					0,
					0,
					0,
					new ElementData(this)})
{
	QSettings settings;
	setDirectory(settings.value("ModFEM/NSSUPGHeat.0/directory", ".").toString());
}

Problem::~Problem()
{
	QSettings settings;
	settings.setValue("ModFEM/NSSUPGHeat.0/directory", directory());
}

QString Problem::directory() const
{
	return m->directory;
}

void Problem::setDirectory(const QString & directory)
{
	if (m->directory != directory) {
		m->directory = directory;
		emit directoryChanged();
	}
}

int Problem::problemId() const
{
	return m->problemId;
}

int Problem::meshId() const
{
	return m->meshId;
}

int Problem::fieldId() const
{
	return m->fieldId;
}

int Problem::solutionCount() const
{
	return m->solutionCount;
}

int Problem::equationCount() const
{
	return m->equationCount;
}

ElementData * Problem::elementData() const
{
	return m->elementData;
}

void Problem::setDirectoryFromURL(const QUrl & url)
{
	setDirectory(url.toLocalFile());
}

void Problem::init()
{
	QByteArray problemDir(QDir(directory()).path().toLocal8Bit());

	char * argv[] = {};
	int argc = 0;
	CUTEHMI_DEBUG("Initializing Navier-Stokes (SUPG) heat problem...");
	utr_set_interactive(problemDir.data(), argc, argv, & m->interactiveInput, & m->interactiveOutput);
	std::cout.flush();

	// utr_set_interactive resets all time counters
	// we can reset any particular time counter at any moment
	utv_time.total = 0.0;
	double time_begin = time_clock();

#ifdef DEBUG
	fprintf(interactive_output, "Starting program in debug mode.\n");
#endif
#ifdef DEBUG_MMM
	fprintf(interactive_output, "Starting mesh module in debug mode.\n");
#endif
#ifdef DEBUG_APM
	fprintf(interactive_output, "Starting approximation module in debug mode.\n");
#endif
#ifdef DEBUG_SIM
	fprintf(interactive_output, "Starting solver interface in debug mode.\n");
#endif
#ifdef DEBUG_LSM
	fprintf(interactive_output, "Starting linear solver (adapter) module in debug mode.\n");
#endif
	std::cout.flush();

	// Initialization of problem data (including mesh and two fields).

	// initialize  heat structures.
	int iaux = pdr_ns_supg_heat_init(problemDir.data(), m->interactiveInput, m->interactiveOutput);
	std::cout.flush();
	if (iaux == EXIT_FAILURE)
		CUTEHMI_CRITICAL("Failed to initialize NS SUPG Heat problem.");
	else {
#ifdef PARALLEL
		if (pcr_my_proc_id() == pcr_print_master()) {
#endif
			// to force writing out current waiting messages
			fflush(m->interactiveOutput);
#ifdef PARALLEL
		}
#endif

		utv_time.problem_initialization += time_clock() - time_begin;

		setProblemId(pdv_ns_supg_heat_current_problem_id);
		setMeshId(pdr_ctrl_i_params(problemId(), 2));
		setFieldId(pdr_ctrl_i_params(problemId(), 3));
		setSolutionCount(pdr_ctrl_i_params(problemId(), 4));
		setEquationCount(pdr_ctrl_i_params(problemId(), 5));

		m->elementData->init(meshId());
	}

	pdv_ns_supg_problem.time.final_step = 1;	// TEMP force single step
	solve();	//temp
}

void Problem::solve()
{
	QByteArray workDirectory(QDir(directory()).path().toLocal8Bit());
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

	pdr_ns_supg_heat_time_integration(workDirectory.data(), m->interactiveInput, m->interactiveOutput);


	timer_all = time_clock() - timer_all;

	fprintf(m->interactiveOutput, "\nExecution time total: %lf\n", timer_all);

	pdv_ns_supg_problem.time.final_step++;	// TEMP force single step
}

void Problem::integrate()
{
	FILE * interactiveInput = m->interactiveInput;
	FILE * interactiveOutput = m->interactiveOutput;

	QtConcurrent::run([ = ]() {
		QByteArray workDirectory(QDir(directory()).path().toLocal8Bit());

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

		pdr_ns_supg_heat_time_integration(workDirectory.data(), interactiveInput, interactiveOutput);


		timer_all = time_clock() - timer_all;

		fprintf(m->interactiveOutput, "\nExecution time total: %lf\n", timer_all);
	});
}

void Problem::writeParaview()
{
//	QDir problemDir(directory());
//	FILE * interactiveInput = stdin;
//	FILE * interactiveOutput = stdout;

//	pdr_heat_write_paraview(problemDir.path().toLocal8Bit().data(), interactiveInput, interactiveOutput);
}

void Problem::setProblemId(int problemId)
{
	if (m->problemId != problemId) {
		m->problemId = problemId;
		emit problemIdChanged();
	}
}

void Problem::setMeshId(int meshId)
{
	if (m->meshId != meshId) {
		m->meshId = meshId;
		emit meshIdChanged();
	}
}

void Problem::setFieldId(int fieldId)
{
	if (m->fieldId != fieldId) {
		m->fieldId = fieldId;
		emit fieldIdChanged();
	}
}

void Problem::setSolutionCount(int solutionCount)
{
	if (m->solutionCount != solutionCount) {
		m->solutionCount = solutionCount;
		emit solutionCountChanged();
	}
}

void Problem::setEquationCount(int equationCount)
{
	if (m->equationCount != equationCount) {
		m->equationCount = equationCount;
		emit equationCountChanged();
	}
}

}
}

/*------------------------------------------------------------
pdr_ns_supg_heat_post_process
------------------------------------------------------------*/
double pdr_ns_supg_heat_post_process(
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output)
{
	double x[3], xg[3];
	int pdeg;			/* degree of polynomial */
	int base;		/* type of basis functions for quadrilaterals */
	int num_shap;			/* number of element shape functions */
	int ndofs;			/* local dimension of the problem */
	double xcoor[3];		/* global coord of gauss point */
	double u_val[PDC_MAXEQ];	/* computed solution */
	double u_x[PDC_MAXEQ];		/* gradient of computed solution */
	double u_y[PDC_MAXEQ];		/* gradient of computed solution */
	double u_z[PDC_MAXEQ];		/* gradient of computed solution */
	double base_phi[APC_MAXELVD];	/* basis functions */
	double base_dphix[APC_MAXELVD];	/* x-derivatives of basis function */
	double base_dphiy[APC_MAXELVD];	/* y-derivatives of basis function */
	double base_dphiz[APC_MAXELVD];	/* y-derivatives of basis function */
	int el_nodes[MMC_MAXELVNO + 1];	/* list of nodes of El */
	double node_coor[3 * MMC_MAXELVNO];	/* coord of nodes of El */
	double dofs_loc[APC_MAXELSD];	/* element solution dofs */
	double dofs_loc2[APC_MAXELSD];	/* element solution dofs */
	int i, j, iel, ki, iaux, mat_num, nel, sol_vec_id, nreq;
	int list_el[20];
	int problem_id, field_id, mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	fprintf(Interactive_output, "Select problem: 1 - ns_supg, 2 - heat:\n");
	fscanf(Interactive_input, "%d", &problem_id);

	i = 3;
	field_id = pdr_ctrl_i_params(problem_id, i);
	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(field_id);
	i = 5;
	nreq = pdr_ctrl_i_params(problem_id, i);

	fprintf(Interactive_output, "Give global coordinates of a point (x,y,z):\n");
	fscanf(Interactive_input, "%lf", &x[0]);
	fscanf(Interactive_input, "%lf", &x[1]);
	fscanf(Interactive_input, "%lf", &x[2]);
	fprintf(Interactive_output, "x=%lf,y=%lf,z=%lf\n", x[0], x[1], x[2]);

	iaux = apr_sol_xglob(field_id, x, 1, list_el, xg, u_val, NULL, NULL, NULL, APC_CLOSE,
					APE_SOL_XGLOB_DEFAULT
					| APE_SOL_XGLOB_MATCH_ALL_ELEMENTS
					| APE_SOL_XGLOB_MATCH_WITH_ADAPTATION);
	if (iaux == 1) {
		fprintf(Interactive_output, "\nSolution at point %.2lf %.2lf %.2lf in element %d:\n\n", xg[0], xg[1], xg[2], list_el[1]);
		for (j = 0; j < nreq; j++) {
			fprintf(Interactive_output, "u_val[%d]=%lf\n", j, u_val[j]);
		}
	}
	else {
		printf("Local coordinates not found within a family in apr_sol_xglob\n");
	}

	return (0);

	fprintf(Interactive_output, "Give element number:\n");
	fscanf(Interactive_input, "%d", &nel);
	fprintf(Interactive_output, "Give local coordinates of a point (x,y,z):\n");
	fscanf(Interactive_input, "%lf", &x[0]);
	fscanf(Interactive_input, "%lf", &x[1]);
	fscanf(Interactive_input, "%lf", &x[2]);
	//fprintf(Interactive_output, "nel=%d, x=%lf, y=%lf, z=%lf\n",
	//        nel,x[0],x[1],x[2]);
	base = apr_get_base_type(field_id, nel);
	pdeg = apr_get_el_pdeg(field_id, nel, &pdeg);
	num_shap = apr_get_el_pdeg_numshap(field_id, nel, &pdeg);
	i = 5;
	nreq = pdr_ctrl_i_params(problem_id, i);
	ndofs = nreq * num_shap;
	/* get the coordinates of the nodes of El in the right order */
	mmr_el_node_coor(mesh_id, nel, el_nodes, node_coor);
	/* get the most recent solution degrees of freedom */
	sol_vec_id = 1;
	for (j = 0; j < nreq; j++) {
		for (i = 0; i < el_nodes[0]; i++) {
			apr_read_ent_dofs(field_id, APC_VERTEX, el_nodes[i + 1], nreq, sol_vec_id, dofs_loc);
			dofs_loc2[j * num_shap + i] = dofs_loc[j];
		}
	}
	/* calculations with jacobian but not on the boundary */
	iaux = 2;
	apr_elem_calc_3D(iaux, nreq, &pdeg, base, x, node_coor, dofs_loc2,
			base_phi, base_dphix, base_dphiy, base_dphiz,
			xcoor, u_val, u_x, u_y, u_z, NULL);
	fprintf(Interactive_output,
			"\nSolution at point %.2lf %.2lf %.2lf in element %d:\n\n",
			x[0], x[1], x[2], nel);
	for (j = 0; j < nreq; j++) {
		fprintf(Interactive_output, "u_val[%d]=%lf\n", j, u_val[j]);
	}
	fprintf(Interactive_output,
			"\nGlobal coordinates of the point:  %.2lf %.2lf %.2lf:\n\n",
			xcoor[0], xcoor[1], xcoor[2]);

	return (1);
}

/*------------------------------------------------------------
pdr_ns_supg_heat_write_profile
------------------------------------------------------------*/
int pdr_ns_supg_heat_write_profile(
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output
)
{
	double x1[3], x2[3];
	int solNr; // solution component ID
	int nSol; // number of solution components
	int nPoints; // number of points along the line
	int problem_id, field_id, mesh_id, i;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	fprintf(Interactive_output, "Select problem: 1 - ns_supg, 2 - heat:\n");
	fscanf(Interactive_input, "%d", &problem_id);

	pdv_ns_supg_heat_current_problem_id = problem_id;
	i = 3;
	field_id = pdr_ctrl_i_params(problem_id, i);
	/* select the corresponding mesh */
	mesh_id = apr_get_mesh_id(field_id);
	//i=5; nSol = pdr_ctrl_i_params(problem_id, i);
	nSol = apr_get_nreq(field_id);

	fprintf(Interactive_output, "Give solution component number (>=0):\n");
	fscanf(Interactive_input, "%d", &solNr);
	fprintf(Interactive_output, "Give number of points (>0):\n");
	fscanf(Interactive_input, "%d", &nPoints);
	fprintf(Interactive_output, "Give global coordinates of a point1 (x,y,z):\n");
	fscanf(Interactive_input, "%lf", &x1[0]);
	fscanf(Interactive_input, "%lf", &x1[1]);
	fscanf(Interactive_input, "%lf", &x1[2]);
	fprintf(Interactive_output, "Give global coordinates of a point2 (x,y,z):\n");
	fscanf(Interactive_input, "%lf", &x2[0]);
	fscanf(Interactive_input, "%lf", &x2[1]);
	fscanf(Interactive_input, "%lf", &x2[2]);
	apr_get_profile(Interactive_output, field_id, solNr, nSol, x1, x2, nPoints);

	return 1;
}



/*---------------------------------------------------------
pdr_change_data - to change some of control data
---------------------------------------------------------*/
void pdr_change_data(
		int Problem_id	/* in: data structure to be used  */
)
{

	/******************************* WARNING ******************************/
	/* in its current version the procedure changes ns_supg and heat data */
	/* separately - consistency must be ensured by the user, e.g.         */
	/* ns_supg->dtime must be equal heat->dtime                           */

	/* local variables */
	pdt_ns_supg_adpts * adpts_ns_supg = &pdv_ns_supg_problem.adpt;
	pdt_ns_supg_times * times_ns_supg = &pdv_ns_supg_problem.time;
	pdt_ns_supg_linss * linss_ns_supg = &pdv_ns_supg_problem.lins;
	pdt_ns_supg_ctrls * ctrls_ns_supg = &pdv_ns_supg_problem.ctrl;
	pdt_heat_adpts * adpts_heat = &pdv_heat_problem.adpt;
	pdt_heat_times * times_heat = &pdv_heat_problem.time;
	pdt_heat_linss * linss_heat = &pdv_heat_problem.lins;
	pdt_heat_ctrls * ctrls_heat = &pdv_heat_problem.ctrl;
	char c, d, pans[100]; /* string variable to read menu */

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper problem data structure */
	if (Problem_id == PDC_NS_SUPG_ID) {

		do {

			do {
				/* define a menu */
				printf("\nChoose a group of data:\n");
				printf("\tc - general control data \n");
				printf("\tt - time integration parameters \n");
				printf("\ta - adaptation parameters \n");
				printf("\tl - linear solver parameters \n");
				printf("\tq - quit changing data for problem %d\n", Problem_id);

				scanf("%s", pans);
				getchar();
			} while ( *pans != 'c' && *pans != 't' && *pans != 'a'
					&& *pans != 'l' && *pans != 'q' && *pans != 'q' );

			c = *pans;

			if (c == 'c') {

				do {

					do {
						/* define a menu */
						printf("\nChoose variable to change:\n");
						printf("\tv - dynamic viscosity\n");
						printf("\tq - quit changing general control data\n");

						scanf("%s", pans);
						getchar();
					} while ( *pans != 's' && *pans != 's' && *pans != 's'
							&& *pans != 'l' && *pans != 'c' && *pans != 'q' );

					d = *pans;

					if (d == 'v') {

						printf("Old value: %lf, new value: ",
								ctrls_ns_supg->dynamic_viscosity);
						scanf("%lg", &ctrls_ns_supg->dynamic_viscosity);
						getchar();

					}

				} while (d != 'q');

			}
			else if (c == 't') {

				do {

					do {
						/* define a menu */
						printf("\nChoose variable to change:\n");
						printf("\ta - method identifier\n");
						printf("\tc - current time-step number\n");
						printf("\td - final time-step number\n");
						printf("\te - current time-step length\n");
						printf("\tf - previous time-step length\n");
						printf("\tg - current time\n");
						printf("\th - final time\n");
						printf("\ti - convergence in time criterion number\n");
						printf("\tj - convergence in time treshold value\n");
						printf("\tk - implicitness parameter alpha (theta)\n");
						printf("\tq - quit changing time integration data\n");

						scanf("%s", pans);
						getchar();
					} while ( *pans != 'a' && *pans != 'c' && *pans != 'c'
							&& *pans != 'd' && *pans != 'e' && *pans != 'f'
							&& *pans != 'g' && *pans != 'h' && *pans != 'i'
							&& *pans != 'j' && *pans != 'k' && *pans != 'q' );

					d = *pans;

					if (d == 'a') {

						printf("Old value: %d, new value: ", times_ns_supg->type);
						scanf("%d", &times_ns_supg->type);
						getchar();

					}
					else if (d == 'c') {

						printf("Old value: %d, new value: ", times_ns_supg->cur_step);
						scanf("%d", &times_ns_supg->cur_step);
						getchar();

					}
					else if (d == 'd') {

						printf("Old value: %d, new value: ", times_ns_supg->final_step);
						scanf("%d", &times_ns_supg->final_step);
						getchar();

					}
					else if (d == 'e') {

						printf("Old value: %lg, new value: ", times_ns_supg->cur_dtime);
						scanf("%lg", &times_ns_supg->cur_dtime);
						getchar();

					}
					else if (d == 'f') {

						printf("Old value: %lg, new value: ", times_ns_supg->prev_dtime);
						scanf("%lg", &times_ns_supg->prev_dtime);
						getchar();

					}
					else if (d == 'g') {

						printf("Old value: %lg, new value: ", times_ns_supg->cur_time);
						scanf("%lg", &times_ns_supg->cur_time);
						getchar();

					}
					else if (d == 'h') {

						printf("Old value: %lg, new value: ", times_ns_supg->final_time);
						scanf("%lg", &times_ns_supg->final_time);
						getchar();

					}
					else if (d == 'i') {

						printf("Old value: %d, new value: ", times_ns_supg->conv_type);
						scanf("%d", &times_ns_supg->conv_type);
						getchar();

					}
					else if (d == 'j') {

						printf("Old value: %lg, new value: ", times_ns_supg->conv_meas);
						scanf("%lg", &times_ns_supg->conv_meas);
						getchar();

					}
					else if (d == 'k') {

						printf("Old value: %lg, new value: ",
								times_ns_supg->alpha);
						scanf("%lg", &times_ns_supg->alpha);
						getchar();

					}

				} while (d != 'q');

			}
			else if (c == 'a') {

				do {

					do {
						/* define a menu */
						printf("\nChoose variable to change:\n");
						printf("\tt - strategy number\n");
						printf("\ti - time interval between adaptations\n");
						printf("\tm - maximal generation level for elements\n");
						printf("\te - global treshold value for adaptation\n");
						printf("\tr - ratio for indicating derefinements\n");
						printf("\tq - quit changing adaptation data\n");

						scanf("%s", pans);
						getchar();
					} while ( *pans != 't' && *pans != 'i' && *pans != 'm' && *pans != 'd'
							&& *pans != 'e' && *pans != 'r' && *pans != 'q' );

					d = *pans;

					if (d == 't') {

						printf("Old value: %d, new value: ", adpts_ns_supg->type);
						scanf("%d", &adpts_ns_supg->type);
						getchar();

					}
					else if (d == 'i') {

						printf("Old value: %d, new value: ", adpts_ns_supg->interval);
						scanf("%d", &adpts_ns_supg->interval);
						getchar();

					}
					else if (d == 'm') {

						printf("Old value: %d, new value: ", adpts_ns_supg->maxgen);
						scanf("%d", &adpts_ns_supg->maxgen);
						getchar();

					}
					else if (d == 'e') {

						printf("Old value: %lg, new value: ", adpts_ns_supg->eps);
						scanf("%lg", &adpts_ns_supg->eps);
						getchar();

					}
					else if (d == 'r') {

						printf("Old value: %lg, new value: ", adpts_ns_supg->ratio);
						scanf("%lg", &adpts_ns_supg->ratio);
						getchar();

					}

				} while (d != 'q');

			}

			else if (c == 'l') {

				do {

					do {
						/* define a menu */
						printf("\nChoose variable to change:\n");
						printf("\tt - solver type\n");
						printf("\ti - maximal number of iterations\n");
						printf("\tc - convergence criterion number\n");
						printf("\te - convergence treshold value\n");
						printf("\tm - monitoring level\n");
						printf("\tq - quit changing linear solver data\n");

						scanf("%s", pans);
						getchar();
					} while ( *pans != 't' && *pans != 'm' && *pans != 'c'
							&& *pans != 'e' && *pans != 'p' && *pans != 'k'
							&& *pans != 'b' && *pans != 'q' && *pans != 'q' );

					d = *pans;

					if (d == 't') {

						printf("Old value: %d, new value: ", linss_ns_supg->type);
						scanf("%d", &linss_ns_supg->type);
						getchar();

					}
					else if (d == 'i') {

						printf("Old value: %d, new value: ", linss_ns_supg->max_iter);
						scanf("%d", &linss_ns_supg->max_iter);
						getchar();

					}
					else if (d == 'c') {

						printf("Old value: %d, new value: ", linss_ns_supg->conv_type);
						scanf("%d", &linss_ns_supg->conv_type);
						getchar();

					}
					else if (d == 'e') {

						printf("Old value: %lg, new value: ", linss_ns_supg->conv_meas);
						scanf("%lg", &linss_ns_supg->conv_meas);
						getchar();

					}
					else if (d == 'm') {

						printf("Old value: %d, new value: ", linss_ns_supg->monitor);
						scanf("%d", &linss_ns_supg->monitor);
						getchar();

					}

				} while (d != 'q');

			}

		} while (c != 'q');

	}
	/* select the proper problem data structure */
	else if (Problem_id == PDC_HEAT_ID) {

		do {

			do {
				/* define a menu */
				printf("\nChoose a group of data:\n");
				printf("\tc - general control data \n");
				printf("\tt - time integration parameters \n");
				printf("\ta - adaptation parameters \n");
				printf("\tl - linear solver parameters \n");
				printf("\tq - quit changing data for problem %d\n", Problem_id);

				scanf("%s", pans);
				getchar();
			} while ( *pans != 'c' && *pans != 't' && *pans != 'a'
					&& *pans != 'l' && *pans != 'q' && *pans != 'q' );

			c = *pans;

			if (c == 'c') {

				do {

					do {
						/* define a menu */
						printf("\nChoose variable to change:\n");
						printf("\tv - \n");
						printf("\tq - quit changing general control data\n");

						scanf("%s", pans);
						getchar();
					} while ( *pans != 'v' && *pans != 'v' && *pans != 'v'
							&& *pans != 'v' && *pans != 'v' && *pans != 'q' );

					d = *pans;

					if (d == 'v') {

						//printf("Old value: %lf, new value: ",
						//	   ctrls_heat->);
						//scanf("%lg",&ctrls_heat->); getchar();

					}

				} while (d != 'q');

			}
			else if (c == 't') {

				do {

					do {
						/* define a menu */
						printf("\nChoose variable to change:\n");
						printf("\ta - method identifier\n");
						printf("\tc - current time-step number\n");
						printf("\td - final time-step number\n");
						printf("\te - current time-step length\n");
						printf("\tf - previous time-step length\n");
						printf("\tg - current time\n");
						printf("\th - final time\n");
						printf("\ti - convergence in time criterion number\n");
						printf("\tj - convergence in time treshold value\n");
						printf("\tk - implicitness parameter alpha (theta)\n");
						printf("\tq - quit changing time integration data\n");

						scanf("%s", pans);
						getchar();
					} while ( *pans != 'a' && *pans != 'c' && *pans != 'c'
							&& *pans != 'd' && *pans != 'e' && *pans != 'f'
							&& *pans != 'g' && *pans != 'h' && *pans != 'i'
							&& *pans != 'j' && *pans != 'k' && *pans != 'q' );

					d = *pans;

					if (d == 'a') {

						printf("Old value: %d, new value: ", times_heat->type);
						scanf("%d", &times_heat->type);
						getchar();

					}
					else if (d == 'c') {

						printf("Old value: %d, new value: ", times_heat->cur_step);
						scanf("%d", &times_heat->cur_step);
						getchar();

					}
					else if (d == 'd') {

						printf("Old value: %d, new value: ", times_heat->final_step);
						scanf("%d", &times_heat->final_step);
						getchar();

					}
					else if (d == 'e') {

						printf("Old value: %lg, new value: ", times_heat->cur_dtime);
						scanf("%lg", &times_heat->cur_dtime);
						getchar();

					}
					else if (d == 'f') {

						printf("Old value: %lg, new value: ", times_heat->prev_dtime);
						scanf("%lg", &times_heat->prev_dtime);
						getchar();

					}
					else if (d == 'g') {

						printf("Old value: %lg, new value: ", times_heat->cur_time);
						scanf("%lg", &times_heat->cur_time);
						getchar();

					}
					else if (d == 'h') {

						printf("Old value: %lg, new value: ", times_heat->final_time);
						scanf("%lg", &times_heat->final_time);
						getchar();

					}
					else if (d == 'i') {

						printf("Old value: %d, new value: ", times_heat->conv_type);
						scanf("%d", &times_heat->conv_type);
						getchar();

					}
					else if (d == 'j') {

						printf("Old value: %lg, new value: ", times_heat->conv_meas);
						scanf("%lg", &times_heat->conv_meas);
						getchar();

					}
					else if (d == 'k') {

						printf("Old value: %lg, new value: ",
								times_heat->alpha);
						scanf("%lg", &times_heat->alpha);
						getchar();

					}

				} while (d != 'q');

			}
			else if (c == 'a') {

				do {

					do {
						/* define a menu */
						printf("\nChoose variable to change:\n");
						printf("\tt - strategy number\n");
						printf("\ti - time interval between adaptations\n");
						printf("\tm - maximal generation level for elements\n");
						printf("\te - global treshold value for adaptation\n");
						printf("\tr - ratio for indicating derefinements\n");
						printf("\tq - quit changing adaptation data\n");

						scanf("%s", pans);
						getchar();
					} while ( *pans != 't' && *pans != 'i' && *pans != 'm' && *pans != 'd'
							&& *pans != 'e' && *pans != 'r' && *pans != 'q' );

					d = *pans;

					if (d == 't') {

						printf("Old value: %d, new value: ", adpts_heat->type);
						scanf("%d", &adpts_heat->type);
						getchar();

					}
					else if (d == 'i') {

						printf("Old value: %d, new value: ", adpts_heat->interval);
						scanf("%d", &adpts_heat->interval);
						getchar();

					}
					else if (d == 'm') {

						printf("Old value: %d, new value: ", adpts_heat->maxgen);
						scanf("%d", &adpts_heat->maxgen);
						getchar();

					}
					else if (d == 'e') {

						printf("Old value: %lg, new value: ", adpts_heat->eps);
						scanf("%lg", &adpts_heat->eps);
						getchar();

					}
					else if (d == 'r') {

						printf("Old value: %lg, new value: ", adpts_heat->ratio);
						scanf("%lg", &adpts_heat->ratio);
						getchar();

					}

				} while (d != 'q');

			}

			else if (c == 'l') {

				do {

					do {
						/* define a menu */
						printf("\nChoose variable to change:\n");
						printf("\tt - solver type\n");
						printf("\ti - maximal number of iterations\n");
						printf("\tc - convergence criterion number\n");
						printf("\te - convergence treshold value\n");
						printf("\tm - monitoring level\n");
						printf("\tq - quit changing linear solver data\n");

						scanf("%s", pans);
						getchar();
					} while ( *pans != 't' && *pans != 'm' && *pans != 'c'
							&& *pans != 'e' && *pans != 'p' && *pans != 'k'
							&& *pans != 'b' && *pans != 'q' && *pans != 'q' );

					d = *pans;

					if (d == 't') {

						printf("Old value: %d, new value: ", linss_heat->type);
						scanf("%d", &linss_heat->type);
						getchar();

					}
					else if (d == 'i') {

						printf("Old value: %d, new value: ", linss_heat->max_iter);
						scanf("%d", &linss_heat->max_iter);
						getchar();

					}
					else if (d == 'c') {

						printf("Old value: %d, new value: ", linss_heat->conv_type);
						scanf("%d", &linss_heat->conv_type);
						getchar();

					}
					else if (d == 'e') {

						printf("Old value: %lg, new value: ", linss_heat->conv_meas);
						scanf("%lg", &linss_heat->conv_meas);
						getchar();

					}
					else if (d == 'm') {

						printf("Old value: %d, new value: ", linss_heat->monitor);
						scanf("%d", &linss_heat->monitor);
						getchar();

					}

				} while (d != 'q');

			}

		} while (c != 'q');

	}

	return;
}

int  pdr_ns_supg_heat_bc_free()
{

	return (0);
}

int  pdr_ns_supg_heat_material_free()
{

	return (0);
}
