#include "../../../include/modfem/qtheat/Problem.hpp"

#include <QSettings>
#include <QDir>
#include <QDataStream>

#include <cstdlib>
#include <iostream>

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>

/* utilities - including simple time measurement library */
#include <modfem/uth_intf.h>		/* USES */
#include <modfem/uth_system.h>
/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>		/* USES */
/* interface for all approximation modules */
#include <modfem/aph_intf.h>		/* USES */
/* interface for all solver modules */
#include <modfem/sih_intf.h>		/* USES */
/* problem dependent module interface */
#include <modfem/pdh_intf.h>		/* USES */
#include <modfem/uth_bc.h>
/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>	/* IMPLEMENTS */
/* interface for thread management modules */
#include <modfem/tmh_intf.h>

#ifdef PARALLEL
/* interface of parallel mesh manipulation modules */
#include "mmph_intf.h"		/* USES */
/* interface for all parallel approximation modules */
#include "apph_intf.h"		/* USES */
/* interface for parallel communication modules */
#include "pch_intf.h"		/* USES */
#endif

/* visualization module */
//#include "mod_fem_viewer.h"	/* USES */

/* problem module's types and functions */
#include <modfem/pd_heat/pdh_heat.h>		/* USES */
/* types and functions related to problem structures */
#include <modfem/pd_heat/pdh_heat_problem.h>
// bc and material header files are included in problem header files
#include <modfem/pd_heat/pdh_heat_materials.h>

/***************************************/
/* DECLARATIONS OF INTERNAL PROCEDURES */
/***************************************/
/* Rules: */
/* - name always begins with pdr_heat */
/* - argument names start uppercase */

/*------------------------------------------------------------
pdr_heat_post_process - simple post-processing
------------------------------------------------------------*/
double pdr_heat_post_process(
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output
);

/*------------------------------------------------------------
pdr_heat_write_profile - to dump a set of values along a line
------------------------------------------------------------*/
int pdr_heat_write_profile(
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

namespace modfem {
namespace qtheat {

Problem::Problem(QObject * parent):
	QObject(parent),
	m(new Members{".",
					0,
					0,
					0,
					0,
					0, {}, std::make_unique<Qt3DRender::QBuffer>(), std::make_unique<Mesh>()})
{
	QSettings settings;
	setDirectory(settings.value("ModFEM/QtHeat.0/directory", ".").toString());

	resetMeshData(); //temp
}

Problem::~Problem()
{
	QSettings settings;
	settings.setValue("ModFEM/QtHeat.0/directory", directory());
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

Mesh * Problem::mesh() const
{
	return m->mesh.get();
}

Qt3DRender::QBuffer * Problem::buffer() const
{
	return m->buffer.get();
}

QByteArray Problem::meshData() const
{
	return m->meshData;
}

void Problem::setDirectoryFromURL(const QUrl & url)
{
	setDirectory(url.toLocalFile());
}

void Problem::init()
{
	QDir problemDir(directory());
	char * argv[] = {};
	int argc = 0;
	CUTEHMI_DEBUG("Initializing heat problem...");
	FILE * interactiveInput = nullptr;
	FILE * interactiveOutput = nullptr;
	utr_set_interactive(problemDir.path().toLocal8Bit().data(), argc, argv, & interactiveInput, & interactiveOutput);
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
	int iaux = pdr_heat_init(problemDir.path().toLocal8Bit().data(), interactiveInput, interactiveOutput);
	std::cout.flush();
	if (iaux == EXIT_FAILURE)
		CUTEHMI_CRITICAL("Failed to initialize heat problem.");
	else {
#ifdef PARALLEL
		if (pcr_my_proc_id() == pcr_print_master()) {
#endif
			printf("\n<<<---***---!!!--- Performance data begin for HEAT problem ---!!!---***--->>>\n");
			printf("Time for initializing problem, mesh, field, material, BC etc. data: %lf\n",
					time_clock() - time_begin);
			printf("<<<---***---!!!-------- Performance data end --------!!!---***--->>>\n\n");
#ifdef PARALLEL
		}
#endif

		utv_time.problem_initialization += time_clock() - time_begin;

		setProblemId(pdv_heat_current_problem_id);
		setMeshId(pdr_ctrl_i_params(problemId(), 2));
		setFieldId(pdr_ctrl_i_params(problemId(), 3));
		setSolutionCount(pdr_ctrl_i_params(problemId(), 4));
		setEquationCount(pdr_ctrl_i_params(problemId(), 5));

		m->mesh->init(meshId());
	}
}

void Problem::resetMeshData()
{
	//temp
//	QDataStream stream(& m->meshData, QIODevice::WriteOnly);
//	stream.setFloatingPointPrecision(QDataStream::SinglePrecision);

//	stream << (qint32) - 5.0f << (qint32) - 5.0f << (qint32)0.0f;
//	stream << (qint32)0.0f << (qint32) - 5.0f << (qint32)0.0f;
//	stream << -5.0f << -5.0f << 0.0f;
//	m->buffer->setData(tempData);
	double point = -5.0;
	m->meshData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = -5.0;
	m->meshData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = -5.0;
	m->meshData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = 0.0;
	m->meshData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = 0.0;
	m->meshData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = -5.0;
	m->meshData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = 0.0;
	m->meshData.append(reinterpret_cast<char *>(& point), sizeof (double));
	emit meshDataChanged();
	//endtemp

	///@todo implement.
}

void Problem::resetBuffer()
{
	//temp
	QByteArray tempData;
	QDataStream stream(& tempData, QIODevice::WriteOnly);

	stream << -5.0f << -5.0f << 0.0f;
	stream << 0.0f << -5.0f << 0.0f;
	stream << -5.0f << -5.0f << 0.0f;
	m->buffer->setData(tempData);
	emit bufferChanged();
	//endtemp

	///@todo implement.
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
pdr_heat_post_process
------------------------------------------------------------*/
double pdr_heat_post_process(
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
	//int el_nodes[MMC_MAXELVNO + 1];	/* list of nodes of El */
	//double node_coor[3 * MMC_MAXELVNO];	/* coord of nodes of El */
	int el_nodes[APC_MAX_GEO_DOFS + 1];      /* list of nodes of El */
	int el_nodes_type[APC_MAX_GEO_DOFS + 1];      /* list of nodes type of El */
	double node_coor[3 * APC_MAX_GEO_DOFS]; /* coord of nodes of El */
	double dofs_loc[APC_MAXELSD];	/* element solution dofs */
	double dofs_loc2[APC_MAXELSD];	/* element solution dofs */
	int i, j, iel, ki, iaux, mat_num, nel, sol_vec_id, nreq;
	int list_el[20];
	int problem_id, field_id, mesh_id;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	problem_id = pdv_heat_current_problem_id;

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
	//mmr_el_node_coor(mesh_id, nel, el_nodes, node_coor);
	apr_get_el_geo_dofs(field_id, nel, el_nodes, el_nodes_type, node_coor);

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
pdr_heat_write_profile
------------------------------------------------------------*/
int pdr_heat_write_profile(
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

	problem_id = pdv_heat_current_problem_id;

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

/*------------------------------------------------------------
pdr_heat_initial_condition
------------------------------------------------------------*/
double pdr_heat_initial_condition(
		int Field_id, // field_id - each problem should know its field id
		double * Xcoor,  // point coordinates
		int Sol_comp_id // solution component
)
{

	/*kbw
	  printf("specified initial temperature at point %lf, %lf, %lf : %lf\n",
		 Xcoor[0], Xcoor[1], Xcoor[2],
		 pdv_heat_problem.ctrl.ambient_temperature);
	/*kew*/
	/*kbw
	  if(fabs(Xcoor[0]-0.6)<0.25 && fabs(Xcoor[1]-0.6)<0.25){
		return(1000.0);
	  }
	  else{
		return (pdv_heat_problem.ctrl.ambient_temperature);
	  }
	/*kew*/
#ifndef PHASE_TRANSFORMATION
	return (pdv_heat_problem.ctrl.ambient_temperature);
#else
	return (pdv_heat_problem.ctrl.austenitizing_temperature);
#endif
}


/*------------------------------------------------------------
pdr_heat_init - read problem data
------------------------------------------------------------*/
int pdr_heat_init(char * Work_dir, FILE * Interactive_input, FILE * Interactive_output)
{

	FILE * testforfile;
	char filename[300], arg[300];
	int nr_sol; // number of solution vectors - determined by time integration
	int mesh_id, field_id, iaux;


	pdr_heat_problem_clear(& pdv_heat_problem);

	nr_sol = 3;
	strcpy(pdv_heat_problem.ctrl.work_dir, Work_dir);
	pdv_heat_problem.ctrl.interactive_input = Interactive_input;
	pdv_heat_problem.ctrl.interactive_output = Interactive_output;
	sprintf(filename, "%s", "problem_heat.dat");
	pdr_heat_problem_read(Work_dir, filename, Interactive_output, & pdv_heat_problem, nr_sol);

	// check data

	fprintf(Interactive_output,
			"\nHEAT problem %d settings :\n", pdv_heat_problem.ctrl.name);

	fprintf(Interactive_output, "\nCONTROL PARAMETERS:\n");

	fprintf(Interactive_output, "\tmesh_type:\t\t\t\t%s\n",
			pdv_heat_problem.ctrl.mesh_type);
	fprintf(Interactive_output, "\tmesh_file_in:\t\t\t\t%s\n",
			pdv_heat_problem.ctrl.mesh_filename);
	fprintf(Interactive_output, "\tmesh_file_out:\t\t\t\t%s\n\n",
			pdv_heat_problem.ctrl.mesh_dmp_filepattern);

	fprintf(Interactive_output, "\tfield_file_in:\t\t\t\t%s\n",
			pdv_heat_problem.ctrl.field_filename);
	fprintf(Interactive_output, "\tfield_file_out:\t\t\t\t%s\n\n",
			pdv_heat_problem.ctrl.field_dmp_filepattern);

#ifdef PHASE_TRANSFORMATION
	fprintf(Interactive_output, "\tphase_transformation_field_file_in:\t%s\n",
			pdv_heat_problem.ctrl.phase_transformation_field_filename);
	fprintf(Interactive_output, "\tphase_transformation_field_file_out:\t%s\n",
			pdv_heat_problem.ctrl.phase_transformation_field_dmp_filepattern);
	fprintf(Interactive_output, "\tphases_field_file_in:\t\t\t%s\n",
			pdv_heat_problem.ctrl.phases_field_filename);
	fprintf(Interactive_output, "\tphases_field_file_out:\t\t\t%s\n",
			pdv_heat_problem.ctrl.phases_field_dmp_filepattern);
#endif

	fprintf(Interactive_output, "\n\tpenalty for Dirichlet BCs:\t\t%lf\n\n",
			pdv_heat_problem.ctrl.penalty);


	/****************************************/
	/* initialization of material data      */
	/****************************************/
	if (strlen(pdv_heat_problem.ctrl.material_filename) == 0) {

		// e.g. for non-dimensional form of NS equations material files are not used
		pdv_heat_problem.ctrl.ref_temperature = -1.0;

	}
	else {

		iaux = pdr_heat_material_read(pdv_heat_problem.ctrl.work_dir,
						pdv_heat_problem.ctrl.material_filename,
						pdv_heat_problem.ctrl.interactive_output);
		if (iaux == EXIT_FAILURE) exit(-1);

		fprintf(Interactive_output, "\nmaterials configuration file:\t\t\t%s\n",
				pdv_heat_problem.ctrl.material_filename);
//    fprintf(Interactive_output, "\tnumber of materials:\t\t\t%d\n\n",
//	    pdv_heat_problem.materials.materials_num);

		fprintf(Interactive_output, "HEAT material database read\n");

	}

	// reference temperature used to indicate whether material database is used
	if (pdv_heat_problem.ctrl.ref_temperature <= 0.0) {

		fprintf(Interactive_output, "\tmaterial data not temperature dependent\n");
		fprintf(Interactive_output, "\n\tthermal_conductivity:\t\t%lf\n",
				pdv_heat_problem.ctrl.thermal_conductivity);
		fprintf(Interactive_output, "\n\tdensity:\t\t%lf\n",
				pdv_heat_problem.ctrl.density);
		fprintf(Interactive_output, "\n\tspecific_heat:\t\t%lf\n",
				pdv_heat_problem.ctrl.specific_heat);

	}
	else {

		fprintf(Interactive_output, "\tmaterial data from material database\n");
		pdv_heat_problem.ctrl.thermal_conductivity = -1.0;
		pdv_heat_problem.ctrl.density = -1.0;
		pdv_heat_problem.ctrl.specific_heat = -1.0;

		fprintf(Interactive_output, "\treference_temperature:\t\t\t%lf\n",
				pdv_heat_problem.ctrl.ref_temperature);

	}

	fprintf(Interactive_output, "\tambient temperature:\t\t\t%lf\n",
			pdv_heat_problem.ctrl.ambient_temperature);
#ifdef PHASE_TRANSFORMATION
	fprintf(Interactive_output, "\taustenitizing temperature:\t\t%lf\n",
			pdv_heat_problem.ctrl.austenitizing_temperature);
#endif

	fprintf(Interactive_output, "\nTIME INTEGRATION PARAMETERS:\n");

	fprintf(Interactive_output, "\ttype:\t\t\t\t\t%d\n",
			pdv_heat_problem.time.type);
	fprintf(Interactive_output, "\timplicitness parameter:\t\t\t%lf\n\n",
			pdv_heat_problem.time.alpha);

	fprintf(Interactive_output, "\tcurrent timestep:\t\t\t%d\n",
			pdv_heat_problem.time.cur_step);
	fprintf(Interactive_output, "\tcurrent time:\t\t\t\t%lf\n",
			pdv_heat_problem.time.cur_time);
	fprintf(Interactive_output, "\n\tcurrent timestep_length:\t\t%lf\n",
			pdv_heat_problem.time.cur_dtime);
	fprintf(Interactive_output, "\tprevious timestep_length:\t\t%lf\n\n",
			pdv_heat_problem.time.prev_dtime);

	fprintf(Interactive_output, "\tfinal time:\t\t\t\t%lf\n",
			pdv_heat_problem.time.final_time);
	fprintf(Interactive_output, "\tfinal timestep:\t\t\t\t%d\n",
			pdv_heat_problem.time.final_step);

	fprintf(Interactive_output, "\n\tconvergence criterion type:\t\t%d\n",
			pdv_heat_problem.time.conv_type);
	fprintf(Interactive_output, "\terror tolerance (n-epsilon):\t\t%lf\n\n",
			pdv_heat_problem.time.conv_meas);

	fprintf(Interactive_output, "\tmonitoring level:\t\t\t%d\n\n",
			pdv_heat_problem.time.monitor);

	fprintf(Interactive_output, "\tgraph_dump_intv:\t\t\t%d\n",
			pdv_heat_problem.time.intv_graph);
	fprintf(Interactive_output, "\tfull_dump_intv:\t\t\t\t%d\n\n",
			pdv_heat_problem.time.intv_dumpout);

	fprintf(Interactive_output, "\nNONLINEAR SOLVER PARAMETERS:\n");

	fprintf(Interactive_output, "\ttype:\t\t\t\t\t%d\n",
			pdv_heat_problem.nonl.type);

	fprintf(Interactive_output, "\tmax_nonl_iter:\t\t\t\t%d\n",
			pdv_heat_problem.nonl.max_iter);

	fprintf(Interactive_output, "\n\tconvergence criterion type:\t\t%d\n",
			pdv_heat_problem.nonl.conv_type);
	fprintf(Interactive_output, "\terror tolerance (k-epsilon):\t\t%lf\n",
			pdv_heat_problem.nonl.conv_meas);
	fprintf(Interactive_output, "\tmonitoring level:\t\t\t%d\n\n",
			pdv_heat_problem.nonl.monitor);


	fprintf(Interactive_output, "\nLINEAR SOLVER PARAMETERS:\n");

	fprintf(Interactive_output, "\n\tsolver type:\t\t\t\t%d\n",
			pdv_heat_problem.lins.type);

	fprintf(Interactive_output, "\tsolver_file:\t\t\t\t%s\n",
			pdv_heat_problem.ctrl.solver_filename);

	if (pdv_heat_problem.lins.type != 0) {
		fprintf(Interactive_output, "\n\tmax_lins_iter:\t\t\t\t%d\n",
				pdv_heat_problem.lins.max_iter);
		fprintf(Interactive_output, "\tconvergence criterion type:\t\t%d\n",
				pdv_heat_problem.lins.conv_type);
		fprintf(Interactive_output, "\terror tolerance:\t\t\t%.15lf\n",
				pdv_heat_problem.lins.conv_meas);
	}

	fprintf(Interactive_output, "\tmonitoring level:\t\t\t%d\n\n",
			pdv_heat_problem.lins.monitor);


	fprintf(Interactive_output, "\nADAPTATION PARAMETERS:\n");

	fprintf(Interactive_output, "\tadapt_type:\t\t\t\t%d\n",
			pdv_heat_problem.adpt.type);
	fprintf(Interactive_output, "\tadapt_interval:\t\t\t\t%d\n",
			pdv_heat_problem.adpt.interval);
	fprintf(Interactive_output, "\tadapt_eps:\t\t\t\t%lf\n",
			pdv_heat_problem.adpt.eps);
	fprintf(Interactive_output, "\tadapt_ratio:\t\t\t\t%lf\n\n",
			pdv_heat_problem.adpt.ratio);

	fprintf(Interactive_output, "\tmonitoring level:\t\t\t%d\n\n",
			pdv_heat_problem.adpt.monitor);


	/****************************************/
	/* initialization of bc data            */
	/****************************************/

	iaux = pdr_heat_bc_read(pdv_heat_problem.ctrl.work_dir,
					pdv_heat_problem.ctrl.bc_filename,
					pdv_heat_problem.ctrl.interactive_output,
					&pdv_heat_problem.bc);
	if (iaux == EXIT_FAILURE) exit(-1);

	fprintf(Interactive_output, "\nboundary conditions configuration file:\t\t%s\n",
			pdv_heat_problem.ctrl.bc_filename);
	fprintf(Interactive_output, "\tnumber of BCs:\t\t\t\t%d\n\n",
			pdr_heat_get_bc_assign_count(&pdv_heat_problem.bc));


	fprintf(Interactive_output, "HEAT BC OK\n");


	/****************************************/
	/* initialization of mesh data          */
	/****************************************/

	mesh_id = utr_initialize_mesh( Interactive_output, Work_dir,
					pdv_heat_problem.ctrl.mesh_type[0],
					pdv_heat_problem.ctrl.mesh_filename);
	pdv_heat_problem.ctrl.mesh_id = mesh_id;

	fprintf(Interactive_output, "Read mesh file %s\n",
			pdv_heat_problem.ctrl.mesh_filename);

	utr_mesh_insert_new_bcnums(Work_dir, mesh_id);
	utr_bc_to_insert_completed();

#ifdef DEBUG
	{
		int currfa = 0;
		int fa_bnum;
		/* check if every boundary has been assigned boundary condtion */
		while (currfa = mmr_get_next_face_all(pdv_heat_problem.ctrl.mesh_id,
								currfa)) {
			fa_bnum = mmr_fa_bc(pdv_heat_problem.ctrl.mesh_id, currfa);
			//fprintf(Interactive_output, "BC HEAT %d set for boundary %d\n",
			//	    pdr_heat_get_bc_type(&pdv_heat_problem.bc, fa_bnum), fa_bnum);
			if (fa_bnum > 0) {
				if ((pdr_heat_get_bc_type(&pdv_heat_problem.bc, fa_bnum) == BC_HEAT_NONE)) {
					fprintf(Interactive_output, "BC HEAT not set for boundary:\t%d\n", fa_bnum);
					fprintf(Interactive_output, "Check bc config file - Exiting.\n");
					exit(-1);
				}
			}
		}
	}
#endif


	mmr_set_max_gen(mesh_id, pdv_heat_problem.adpt.maxgen);
	mmr_set_max_gen_diff(mesh_id, 1); // one irregularity of meshes enforced

	fprintf(Interactive_output, "\nAfter reading initial mesh data.\n\n");
	fprintf(Interactive_output,
			"Mesh entities (number of active, maximal index):\n");
	fprintf(Interactive_output, "Elements: nrel %d, nmel %d\n",
			mmr_get_nr_elem(mesh_id), mmr_get_max_elem_id(mesh_id));
	fprintf(Interactive_output, "Faces:    nrfa %d, nmfa %d\n",
			mmr_get_nr_face(mesh_id), mmr_get_max_face_id(mesh_id));
	fprintf(Interactive_output, "Edges:    nred %d, nmed %d\n",
			mmr_get_nr_edge(mesh_id), mmr_get_max_edge_id(mesh_id));
	fprintf(Interactive_output, "Nodes:    nrno %d, nmno %d\n",
			mmr_get_nr_node(mesh_id), mmr_get_max_node_id(mesh_id));

	fprintf(Interactive_output,
			"\nMaximal generation level set to %d, maximal generation difference set to %d\n",
			mmr_get_max_gen(mesh_id), mmr_get_max_gen_diff(mesh_id));

	/****************************************/
	/* initialization of approximation field data */
	/****************************************/
	int pdeg;

	// Get PDEG values
	pdeg = pdv_heat_problem.ctrl.pdeg;

	// first heat field
	if (strcmp(pdv_heat_problem.ctrl.field_filename, "z") == 0) {
		fprintf(Interactive_output,
				"\nInitializing heat field with 0\n");
		// 's' - for standard continuous basis functions
		// 'z' - for zeroing field values
		pdv_heat_problem.ctrl.field_id = utr_initialize_field(
						Interactive_output, 's', 'z',
						pdv_heat_problem.ctrl.mesh_id,
						pdv_heat_problem.ctrl.nreq,
						pdv_heat_problem.ctrl.nr_sol, pdeg, NULL, NULL);
	}
	else if (strcmp(pdv_heat_problem.ctrl.field_filename, "i") == 0) {
		fprintf(Interactive_output,
				"\nInitializing heat field with initial_condition function\n");
		// 's' - for standard continuous basis functions
		// 'i' - for initializing using function pdr_heat_heat_initial_condition
		pdv_heat_problem.ctrl.field_id = utr_initialize_field(
						Interactive_output, 's', 'i',
						pdv_heat_problem.ctrl.mesh_id,
						pdv_heat_problem.ctrl.nreq,
						pdv_heat_problem.ctrl.nr_sol, pdeg, NULL,
						pdr_heat_initial_condition);
	}
	else {
		sprintf(arg, "%s/%s", Work_dir, pdv_heat_problem.ctrl.field_filename);
		testforfile = fopen(arg, "r");
		if (testforfile != NULL) {
			fclose(testforfile);
			fprintf(Interactive_output, "\nInput field file %s.",
					pdv_heat_problem.ctrl.field_filename);
			// 's' - for standard continuous basis functions
			// 'i' - for initializing using function pdr_heat_heat_initial_condition
			pdv_heat_problem.ctrl.field_id = utr_initialize_field(
							Interactive_output, 's', 'r',
							pdv_heat_problem.ctrl.mesh_id,
							pdv_heat_problem.ctrl.nreq,
							pdv_heat_problem.ctrl.nr_sol, pdeg, arg, NULL);
		} else {
			fprintf(Interactive_output,
					"\nInput field file %s not found - setting field to 0\n",
					pdv_heat_problem.ctrl.field_filename);
			// 's' - for standard continuous basis functions
			// 'z' - for zeroing field values
			pdv_heat_problem.ctrl.field_id = utr_initialize_field(
							Interactive_output, 's', 'z',
							pdv_heat_problem.ctrl.mesh_id,
							pdv_heat_problem.ctrl.nreq,
							pdv_heat_problem.ctrl.nr_sol, pdeg, NULL, NULL);
		}
	}
	apr_check_field(pdv_heat_problem.ctrl.field_id);

#ifdef PHASE_TRANSFORMATION
	// second phase_transformation field
	// zeroing mean no phase transformation (PT_HEAT_NONE)
	if (strcmp(pdv_heat_problem.ctrl.phase_transformation_field_filename, "z") == 0) {
		// 's' - for standard continuous basis functions
		// 'z' - for zeroing field values
		pdv_heat_problem.ctrl.phase_transformation_field_id = utr_initialize_field(
						Interactive_output, 's', 'z', pdv_heat_problem.ctrl.mesh_id,
						pdv_heat_problem.ctrl.phase_transformation_nreq, pdv_heat_problem.ctrl.phase_transformation_nr_sol,
						pdeg, NULL, NULL);
	}
	else if (strcmp(pdv_heat_problem.ctrl.phase_transformation_field_filename, "i") == 0) {
		// 's' - for standard continuous basis functions
		// 'i' - for initializing using function pdr_plast_flow_initial_condition
		pdv_heat_problem.ctrl.phase_transformation_field_id = utr_initialize_field(
						Interactive_output, 's', 'i', pdv_heat_problem.ctrl.mesh_id,
						pdv_heat_problem.ctrl.phase_transformation_nreq, pdv_heat_problem.ctrl.phase_transformation_nr_sol,
						pdeg, NULL, pdr_heat_initial_condition);
	}
	else {
		sprintf(arg, "%s/%s", Work_dir, pdv_heat_problem.ctrl.phase_transformation_field_filename);
		testforfile = fopen(arg, "r");
		if (testforfile != NULL) {
			fclose(testforfile);
			fprintf(Interactive_output, "Input field file %s.",
					pdv_heat_problem.ctrl.phase_transformation_field_filename);
			// 's' - for standard continuous basis functions
			// 'i' - for initializing using function pdr_strain_initial_condition
			pdv_heat_problem.ctrl.phase_transformation_field_id = utr_initialize_field(
							Interactive_output, 's', 'r', pdv_heat_problem.ctrl.mesh_id,
							pdv_heat_problem.ctrl.phase_transformation_nreq, pdv_heat_problem.ctrl.phase_transformation_nr_sol,
							pdeg, arg, NULL);
		} else {
			fprintf(Interactive_output,
					"Input field file %s not found - setting field to 0\n", pdv_heat_problem.ctrl.phase_transformation_field_filename);
			// 's' - for standard continuous basis functions
			// 'z' - for zeroing field values
			pdv_heat_problem.ctrl.phase_transformation_field_id = utr_initialize_field(
							Interactive_output, 's', 'z', pdv_heat_problem.ctrl.mesh_id,
							pdv_heat_problem.ctrl.phase_transformation_nreq, pdv_heat_problem.ctrl.phase_transformation_nr_sol,
							pdeg, NULL, NULL);
		}
	}
	apr_check_field(pdv_heat_problem.ctrl.phase_transformation_field_id);

	// third phases field
	if (strcmp(pdv_heat_problem.ctrl.phases_field_filename, "z") == 0) {
		// 's' - for standard continuous basis functions
		// 'z' - for zeroing field values
		pdv_heat_problem.ctrl.phases_field_id = utr_initialize_field(
						Interactive_output, 's', 'z', pdv_heat_problem.ctrl.mesh_id,
						pdv_heat_problem.ctrl.phases_nreq, pdv_heat_problem.ctrl.phases_nr_sol,
						pdeg, NULL, NULL);
	}
	else if (strcmp(pdv_heat_problem.ctrl.phases_field_filename, "i") == 0) {
		// 's' - for standard continuous basis functions
		// 'i' - for initializing using function pdr_plast_flow_initial_condition
		pdv_heat_problem.ctrl.phases_field_id = utr_initialize_field(
						Interactive_output, 's', 'i', pdv_heat_problem.ctrl.mesh_id,
						pdv_heat_problem.ctrl.phases_nreq, pdv_heat_problem.ctrl.phases_nr_sol,
						pdeg, NULL, pdr_heat_initial_condition);
	}
	else {
		sprintf(arg, "%s/%s", Work_dir, pdv_heat_problem.ctrl.phases_field_filename);
		testforfile = fopen(arg, "r");
		if (testforfile != NULL) {
			fclose(testforfile);
			fprintf(Interactive_output, "Input field file %s.",
					pdv_heat_problem.ctrl.phases_field_filename);
			// 's' - for standard continuous basis functions
			// 'i' - for initializing using function pdr_strain_initial_condition
			pdv_heat_problem.ctrl.phases_field_id = utr_initialize_field(
							Interactive_output, 's', 'r', pdv_heat_problem.ctrl.mesh_id,
							pdv_heat_problem.ctrl.phases_nreq, pdv_heat_problem.ctrl.phases_nr_sol,
							pdeg, arg, NULL);
		} else {
			fprintf(Interactive_output,
					"Input field file %s not found - setting field to 0\n", pdv_heat_problem.ctrl.phases_field_filename);
			// 's' - for standard continuous basis functions
			// 'z' - for zeroing field values
			pdv_heat_problem.ctrl.phases_field_id = utr_initialize_field(
							Interactive_output, 's', 'z', pdv_heat_problem.ctrl.mesh_id,
							pdv_heat_problem.ctrl.phases_nreq, pdv_heat_problem.ctrl.phases_nr_sol,
							pdeg, NULL, NULL);
		}
	}
	apr_check_field(pdv_heat_problem.ctrl.phases_field_id);

	fprintf(Interactive_output, "\nInitializing phases field\n");
	pdr_phases_field_init(pdv_heat_problem.ctrl.field_id);

	//****** ACHTUNG this line must be done for last initialized field ***
	pdv_heat_problem.ctrl.fields_nb = pdv_heat_problem.ctrl.phases_field_id;
#endif

	// third heat_dtdt field
	/* fprintf(Interactive_output,  */
	/* 	    "\nInitializing heat_dtdt field with 0\n");  */
	/* pdv_heat_dtdt_problem.ctrl.field_id = utr_initialize_field( */
	/*                           Interactive_output, 's', 'z',  */
	/* 			    pdv_heat_dtdt_problem.ctrl.mesh_id,  */
	/* 			    pdv_heat_dtdt_problem.ctrl.nreq,  */
	/* 			    pdv_heat_dtdt_problem.ctrl.nr_sol, pdeg,NULL,NULL); */

	/* fprintf(Interactive_output, "\n\tAS: num_dof = % d", apr_get_ent_nrdofs(pdv_heat_dtdt_problem.ctrl.field_id,APC_VERTEX,1)); */
	/* fprintf(Interactive_output, "\n\tAS: nreq = % d", apr_get_nreq(pdv_heat_dtdt_problem.ctrl.field_id)); */
	/* fprintf(Interactive_output, "\n\tAS: numshap = % d", apr_get_ent_numshap(pdv_heat_dtdt_problem.ctrl.field_id, APC_VERTEX, 1)); */
	/* fprintf(Interactive_output, "\n\tAS: pdeg = % d", apr_get_ent_pdeg(pdv_heat_dtdt_problem.ctrl.field_id, APC_VERTEX, 1)); */
	/* fprintf(Interactive_output, "\n\tAS: field_id = %d, mesh_id = %d, nreq = %d, nr_sol = %d, pdeg = %d\n", */
	/* 	  pdv_heat_dtdt_problem.ctrl.field_id, */
	/* 	  pdv_heat_dtdt_problem.ctrl.mesh_id,  */
	/* 	  pdv_heat_dtdt_problem.ctrl.nreq,  */
	/* 	  pdv_heat_dtdt_problem.ctrl.nr_sol, pdeg); */

	/* apr_check_field(pdv_heat_dtdt_problem.ctrl.field_id); */

	return (0);
}

/*---------------------------------------------------------
pdr_change_data - to change some of control data
---------------------------------------------------------*/
void pdr_change_data(
		int Problem_id	/* in: data structure to be used  */
)
{


	/* local variables */
	pdt_heat_adpts * adpts_heat = &pdv_heat_problem.adpt;
	pdt_heat_times * times_heat = &pdv_heat_problem.time;
	pdt_heat_linss * linss_heat = &pdv_heat_problem.lins;
	pdt_heat_ctrls * ctrls_heat = &pdv_heat_problem.ctrl;
	char c, d, pans[100]; /* string variable to read menu */

	/*++++++++++++++++ executable statements ++++++++++++++++*/


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


	return;
}
