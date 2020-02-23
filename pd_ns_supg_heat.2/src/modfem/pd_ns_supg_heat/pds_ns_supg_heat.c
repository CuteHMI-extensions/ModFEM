#include <modfem/pd_ns_supg_heat/pdh_ns_supg_heat.h>

#include <modfem/aph_intf.h> /* USES */
#include <modfem/uth_log.h>

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/**************************************/
/* GLOBAL CONSTANTS                   */
/**************************************/
/* Rules:
/* - constants always uppercase and start with PDC_ */

/* from pdh_inf.h */
const int PDC_ELEMENT = APC_ELEMENT;
const int PDC_FACE = APC_FACE;
const int PDC_EDGE = APC_EDGE;
const int PDC_VERTEX = APC_VERTEX;
const int PDC_MIXED_ELEMENT = APC_MIXED_ELEMENT;
const int PDC_MIXED_FACE = APC_MIXED_FACE;
const int PDC_MIXED_EDGE = APC_MIXED_EDGE;
const int PDC_MIXED_VERTEX = APC_MIXED_VERTEX;

const int PDC_NO_COMP = APC_NO_COMP;  /* do not compute stiff mat and rhs vect */
const int PDC_COMP_SM = APC_COMP_SM;  /* compute entries to stiff matrix only */
const int PDC_COMP_RHS = APC_COMP_RHS;/* compute entries to rhs vector only */
const int PDC_COMP_BOTH = APC_COMP_BOTH; /* compute entries for sm and rhsv */


/**************************************/
/* GLOBAL VARIABLES                   */
/**************************************/
/* Rules:
/* - name always begins with pdv_ */

// ID of the current problem
// on purpose initialized to 0 which is wrong value !
// later should be replaced by one of the two proper values:
// ns_supg -> problem_id = PDC_NS_SUPG_ID = 1
// heat -> problem_id = PDC_HEAT_ID = 2
int pdv_ns_supg_heat_current_problem_id = 0;	/* ID of the current problem */
// problem structure for ns_supg module
pdt_ns_supg_problem pdv_ns_supg_problem;
// problem structure for heat module
pdt_heat_problem pdv_heat_problem;

#ifdef PARALLEL
int pdv_ns_supg_exchange_table_index = -1;
int pdv_heat_exchange_table_index = -1;
#endif

/*---------------------------------------------------------
pdr_err_indi - to return error indicator for an element
----------------------------------------------------------*/
double pdr_err_indi(		/* returns error indicator for an element */
  int Problem_id,	/* in: data structure to be used  */
  int Mode,	/* in: mode of operation */
  int El	/* in: element number */
	)
{

  if (Mode == PDC_ADAPT_EXPL) {

	return pdr_ns_supg_heat_err_indi_explicit(Problem_id, El);

  } else if (Mode == PDC_ADAPT_ZZ) {

	return pdr_ns_supg_heat_err_indi_ZZ(Problem_id, El);

  } else {

	printf("Unknown error indicator in pdr_err_indi!\n");

  }

  return (0.0);
}

/*------------------------------------------------------------
  pdr_get_problem_structure - to return pointer to problem structure
------------------------------------------------------------*/
void* pdr_get_problem_structure(int Problem_id)
{
  if(Problem_id==PDC_NS_SUPG_ID){
	return (&pdv_ns_supg_problem);
  } else if(Problem_id==PDC_HEAT_ID){
	return (&pdv_heat_problem);
  } else{
	printf("Wrong problem_id in pdr_get_problem_structure!");
	exit(1);
  }

  return (NULL);
}

/*------------------------------------------------------------
  pdr_problem_name - to return the problem's name
------------------------------------------------------------*/
int pdr_problem_name(
				  /* returns: >=0 - success code, <0 - error code */
  int Problem_id,  /* in: problem ID or
				 PDC_USE_CURRENT_PROBLEM_ID for the current problem */
  char* Problem_name /* out: the name of the problem solved */
  )
{
  if(Problem_id==PDC_NS_SUPG_ID){
	strcpy(Problem_name, pdv_ns_supg_problem.ctrl.name);
  } else if(Problem_id==PDC_HEAT_ID){
	strcpy(Problem_name, pdv_heat_problem.ctrl.name);
  } else{
	printf("Wrong problem_id in pdr_problem_name!");
	exit(1);
  }
  return 0;
}


/*------------------------------------------------------------
pdr_ctrl_i_params - to return one of control parameters
------------------------------------------------------------*/
int pdr_ctrl_i_params(int Problem_id, int Num)
{

  pdt_heat_ctrls *ctrl_heat = &pdv_heat_problem.ctrl;
  pdt_ns_supg_ctrls *ctrl_ns_supg = &pdv_ns_supg_problem.ctrl;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Problem_id==PDC_NS_SUPG_ID){

	/* if (Num == 1) { */
	/*   return (ctrl_ns_supg->name); */ // changed to string
	/* } else  */
	if (Num == 2) {
	  return (ctrl_ns_supg->mesh_id);
	} else if (Num == 3) {
	  return (ctrl_ns_supg->field_id);
	} else if (Num == 4) {
	  return (ctrl_ns_supg->nr_sol);
	} else if (Num == 5) {
	  return (ctrl_ns_supg->nreq);
	} else if (Num == 6) {
	  return (ctrl_ns_supg->solver_id);
	} else {

	}
  } else if(Problem_id==PDC_HEAT_ID) {

	/* if (Num == 1) { */
	/*   return (ctrl_heat->name); */ // changed to string
	/* } else  */
	if (Num == 2) {
	  return (ctrl_heat->mesh_id);
	} else if (Num == 3) {
	  return (ctrl_heat->field_id);
	} else if (Num == 4) {
	  return (ctrl_heat->nr_sol);
	} else if (Num == 5) {
	  return (ctrl_heat->nreq);
	} else if (Num == 6) {
	  return (ctrl_heat->solver_id);
	} else {

	}


  }

  return (-1);
}

/*------------------------------------------------------------
pdr_ctrl_d_params - to return one of control parameters
------------------------------------------------------------*/
double pdr_ctrl_d_params(int Problem_id, int Num)
{
  pdt_heat_ctrls *ctrl_heat = &pdv_heat_problem.ctrl;
  pdt_ns_supg_ctrls *ctrl_ns_supg = &pdv_ns_supg_problem.ctrl;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Problem_id==PDC_NS_SUPG_ID){

	if (Num == 10) {
	  return (ctrl_ns_supg->gravity_field[0]);
	} else if (Num == 11) {
	  return (ctrl_ns_supg->gravity_field[1]);
	} else if (Num == 12) {
	  return (ctrl_ns_supg->gravity_field[2]);
	} else if (Num == 20) {
	  return (ctrl_ns_supg->ref_temperature);
	}
  } else if(Problem_id==PDC_HEAT_ID){
	if (Num == 20) {
	  return (ctrl_heat->ref_temperature);
	}
  }

  return(0.0);
}

/*------------------------------------------------------------
pdr_adapt_i_params - to return parameters of adaptation
------------------------------------------------------------*/
int pdr_adapt_i_params(int Problem_id, int Num)
{

  pdt_heat_adpts *adpts_heat = &pdv_heat_problem.adpt;
  pdt_ns_supg_adpts *adpts_ns_supg = &pdv_ns_supg_problem.adpt;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Problem_id==PDC_NS_SUPG_ID){
	if (Num == 1)
	  return (adpts_ns_supg->type);
	else if (Num == 2)
	  return (adpts_ns_supg->interval);
	else if (Num == 3)
	  return (adpts_ns_supg->maxgen);
	else if (Num == 7)
	  return (adpts_ns_supg->monitor);
	else {
	  printf("Wrong parameter number in adapt_i_params!");
	  exit(1);
	}
  } else if(Problem_id==PDC_HEAT_ID){
	if (Num == 1)
	  return (adpts_heat->type);
	else if (Num == 2)
	  return (adpts_heat->interval);
	else if (Num == 3)
	  return (adpts_heat->maxgen);
	else if (Num == 7)
	  return (adpts_heat->monitor);
	else {
	  printf("Wrong parameter number in adapt_i_params!");
	  exit(1);
	}
  } else{
	printf("Wrong problem_id in adapt_i_params!");
	exit(1);
  }


  return (-1);
}


/*------------------------------------------------------------
pdr_adapt_d_params - to return parameters of adaptation
------------------------------------------------------------*/
double pdr_adapt_d_params(int Problem_id, int Num)
{


  pdt_heat_adpts *adpts_heat = &pdv_heat_problem.adpt;
  pdt_ns_supg_adpts *adpts_ns_supg = &pdv_ns_supg_problem.adpt;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Problem_id==PDC_NS_SUPG_ID){
	if (Num == 5)
	  return (adpts_ns_supg->eps);
	else if (Num == 6)
	  return (adpts_ns_supg->ratio);
	else {
	  printf("Wrong parameter number in adapt_d_params!");
	  exit(1);
	}
  } else if(Problem_id==PDC_HEAT_ID){
	if (Num == 5)
	  return (adpts_heat->eps);
	else if (Num == 6)
	  return (adpts_heat->ratio);
	else {
	  printf("Wrong parameter number in adapt_d_params!");
	  exit(1);
	}
  } else{
	printf("Wrong problem_id in adapt_i_params!");
	exit(1);
  }

  return (-1);
}

/*------------------------------------------------------------
pdr_time_i_params - to return parameters of timeation
------------------------------------------------------------*/
int pdr_time_i_params(int Problem_id, int Num)
{

  pdt_heat_times *times_heat = &pdv_heat_problem.time;
  pdt_ns_supg_times *times_ns_supg = &pdv_ns_supg_problem.time;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Problem_id==PDC_NS_SUPG_ID){
	if (Num == 1)
	  return (times_ns_supg->type);
	else if (Num == 3)
	  return (times_ns_supg->cur_step);
	else if (Num == 4)
	  return (times_ns_supg->final_step);
	else if (Num == 9)
	  return (times_ns_supg->conv_type);
	else if (Num == 11)
	  return (times_ns_supg->monitor);
	else {
	  printf("Wrong parameter number in time_i_params!");
	  exit(1);
	}
  } else if(Problem_id==PDC_HEAT_ID){
	if (Num == 1)
	  return (times_heat->type);
	else if (Num == 3)
	  return (times_heat->cur_step);
	else if (Num == 4)
	  return (times_heat->final_step);
	else if (Num == 9)
	  return (times_heat->conv_type);
	else if (Num == 11)
	  return (times_heat->monitor);
	else {
	  printf("Wrong parameter number in time_i_params!");
	  exit(1);
	}
  } else{
	printf("Wrong problem_id in time_i_params!");
	exit(1);
  }


  return (-1);
}


/*------------------------------------------------------------
pdr_time_d_params - to return parameters of timeation
------------------------------------------------------------*/
double pdr_time_d_params(int Problem_id, int Num)
{


  pdt_heat_times *times_heat = &pdv_heat_problem.time;
  pdt_ns_supg_times *times_ns_supg = &pdv_ns_supg_problem.time;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Problem_id==PDC_NS_SUPG_ID){
	if (Num == 2)
	  return (times_ns_supg->alpha);
	else {
	  printf("Wrong parameter number in time_d_params!");
	  exit(1);
	}
  } else if(Problem_id==PDC_HEAT_ID){
	if (Num == 2)
	  return (times_heat->alpha);
	else {
	  printf("Wrong parameter number in time_d_params!");
	  exit(1);
	}
  } else{
	printf("Wrong problem_id in time_i_params!");
	exit(1);
  }

  return (-1);
}

/*---------------------------------------------------------
pdr_set_time_i_params - to change parameters of time discretization
---------------------------------------------------------*/
void pdr_set_time_i_params(
		int Problem_id,	     /* in: data structure to be used  */
	int Num,             /* in: parameter number in control structure */
	int Value            /* in: parameter value */
	)
{
/* auxiliary variables */
  pdt_ns_supg_times *times_ns_supg = &pdv_ns_supg_problem.time;
  pdt_heat_times *times_heat = &pdv_heat_problem.time;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  if(Problem_id==PDC_NS_SUPG_ID){

	if(Num==1) times_ns_supg->type=Value;
	else if(Num==2) times_ns_supg->cur_step=Value;
	else if(Num==3) times_ns_supg->final_step=Value;
	else if(Num==8) times_ns_supg->conv_type=Value;
	else if(Num==10) times_ns_supg->monitor=Value;
	else if(Num==11) times_ns_supg->intv_dumpout=Value;
	else if(Num==12) times_ns_supg->intv_graph=Value;
	else {
	  printf("Wrong parameter number in set_time_i_params!");
	  exit(1);
	}

  }
  else if(Problem_id==PDC_HEAT_ID){

	if(Num==1) times_heat->type=Value;
	else if(Num==2) times_heat->cur_step=Value;
	else if(Num==3) times_heat->final_step=Value;
	else if(Num==8) times_heat->conv_type=Value;
	else if(Num==10) times_heat->monitor=Value;
	else if(Num==11) times_heat->intv_dumpout=Value;
	else if(Num==12) times_heat->intv_graph=Value;
	else {
	  printf("Wrong parameter number in set_time_i_params!");
	  exit(1);
	}

  }

  return;
}

/*---------------------------------------------------------
pdr_set_time_d_params - to change parameters of time discretization
---------------------------------------------------------*/
void pdr_set_time_d_params(
		int Problem_id,	     /* in: data structure to be used  */
	int Num,             /* in: parameter number in control structure */
	double Value         /* in: parameter value */
	)
{
/* auxiliary variables */
  pdt_ns_supg_times *times_ns_supg = &pdv_ns_supg_problem.time;
  pdt_heat_times *times_heat = &pdv_heat_problem.time;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  if(Problem_id==PDC_NS_SUPG_ID){

	if(Num==4) times_ns_supg->cur_time=Value;
	else if(Num==5) times_ns_supg->final_time=Value;
	else if(Num==6) times_ns_supg->cur_dtime=Value;
	else if(Num==7) times_ns_supg->prev_dtime=Value;
	else if(Num==9) times_ns_supg->conv_meas=Value;
	else {
	  printf("Wrong parameter number in set_time_d_params!");
	  exit(1);
	}

  }
  else if(Problem_id==PDC_HEAT_ID){

	if(Num==4) times_heat->cur_time=Value;
	else if(Num==5) times_heat->final_time=Value;
	else if(Num==6) times_heat->cur_dtime=Value;
	else if(Num==7) times_heat->prev_dtime=Value;
	else if(Num==9) times_heat->conv_meas=Value;
	else {
	  printf("Wrong parameter number in set_time_d_params!");
	  exit(1);
	}

  }

  return;
}


/*---------------------------------------------------------
pdr_lins_i_params - to return parameters of linear equations solver
---------------------------------------------------------*/
int pdr_lins_i_params( /* returns: integer linear solver parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_ns_supg_linss *linss_ns_supg = &pdv_ns_supg_problem.lins;
  pdt_heat_linss *linss_heat = &pdv_heat_problem.lins;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  if(Problem_id==PDC_NS_SUPG_ID){

	if(Num==1) return(linss_ns_supg->type);
	else if(Num==2) return(linss_ns_supg->max_iter);
	else if(Num==3) return(linss_ns_supg->conv_type);
	else if(Num==5) return(linss_ns_supg->monitor);
	else {
	  printf("Wrong parameter number in lins_i_params!");
	  exit(-1);
	}

  }
  else if(Problem_id==PDC_HEAT_ID){

	if(Num==1) return(linss_heat->type);
	else if(Num==2) return(linss_heat->max_iter);
	else if(Num==3) return(linss_heat->conv_type);
	else if(Num==5) return(linss_heat->monitor);
	else {
	  printf("Wrong parameter number in lins_i_params!");
	  exit(-1);
	}

  }

/* error condition - that point should not be reached */
  return(-1);
}

/*---------------------------------------------------------
pdr_lins_d_params - to return parameters of linear equations solver
---------------------------------------------------------*/
double pdr_lins_d_params( /* returns: real linear solver parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_ns_supg_linss *linss_ns_supg = &pdv_ns_supg_problem.lins;
  pdt_heat_linss *linss_heat = &pdv_heat_problem.lins;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  if(Problem_id==PDC_NS_SUPG_ID){


	if(Num==4) return(linss_ns_supg->conv_meas);
	else {
	  printf("Wrong parameter number in lins_d_params!");
	  exit(-1);
	}

  }
  else if(Problem_id==PDC_HEAT_ID){


	if(Num==4) return(linss_heat->conv_meas);
	else {
	  printf("Wrong parameter number in lins_d_params!");
	  exit(-1);
	}

  }

/* error condition - that point should not be reached */
  return(-1);
}

/*------------------------------------------------------------
pdr_ns_supg_initial_condition
------------------------------------------------------------*/
double pdr_ns_supg_initial_condition(
  int Field_id, // field_id - each problem should know its field id
  double *Xcoor,   // point coordinates
  int Sol_comp_id // solution component
)
{

  return (0.0);
}

/*------------------------------------------------------------
pdr_heat_initial_condition
------------------------------------------------------------*/
double pdr_heat_initial_condition(
  int Field_id, // field_id - each problem should know its field id
  double *Xcoor,   // point coordinates
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
	return (pdv_heat_problem.ctrl.ambient_temperature);
}

/*------------------------------------------------------------
pdr_ns_supg_heat_init - read problem data
------------------------------------------------------------*/
int pdr_ns_supg_heat_init(
  char *Work_dir,
  FILE *Interactive_input,
  FILE *Interactive_output
)
{

  FILE *testforfile;
  char filename[300], arg[300];
  int nr_sol; // number of solution vectors - determined by time integration
  int mesh_id, field_id, iaux;

  // ns_supg uses two problem structures - one from ns_supg
  pdv_ns_supg_heat_current_problem_id = PDC_NS_SUPG_ID;
  pdr_ns_supg_problem_clear(&pdv_ns_supg_problem);

  // and the second from heat
  pdv_ns_supg_heat_current_problem_id = PDC_HEAT_ID;
  pdr_heat_problem_clear(&pdv_heat_problem);


  // first we read ns_supg problem structure
  pdv_ns_supg_heat_current_problem_id = PDC_NS_SUPG_ID;
  nr_sol = 3;
  strcpy(pdv_ns_supg_problem.ctrl.work_dir, Work_dir);
  pdv_ns_supg_problem.ctrl.interactive_input = Interactive_input;
  pdv_ns_supg_problem.ctrl.interactive_output = Interactive_output;
  sprintf(filename, "%s", "problem_ns_supg.dat");

  pdr_ns_supg_problem_read(Work_dir, filename, Interactive_output,
			   &pdv_ns_supg_problem, nr_sol);


  // second we read heat problem structure
  pdv_ns_supg_heat_current_problem_id = PDC_HEAT_ID;
  nr_sol = 3;
  strcpy(pdv_heat_problem.ctrl.work_dir, Work_dir);
  pdv_heat_problem.ctrl.interactive_input = Interactive_input;
  pdv_heat_problem.ctrl.interactive_output = Interactive_output;
  sprintf(filename, "%s", "problem_heat.dat");
  pdr_heat_problem_read(Work_dir, filename, Interactive_output,
			   &pdv_heat_problem, nr_sol);

  // check data

  fprintf(Interactive_output,
	  "\nNS_SUPG problem %s and HEAT problem %s settings :\n",
	  pdv_ns_supg_problem.ctrl.name, pdv_heat_problem.ctrl.name);

  fprintf(Interactive_output, "\nCONTROL PARAMETERS:\n");

  fprintf(Interactive_output, "\n\tmesh_type:\t\t\t\t%s\n",
	  pdv_ns_supg_problem.ctrl.mesh_type);
  if(0 !=  strcmp(pdv_heat_problem.ctrl.mesh_type,
		 pdv_ns_supg_problem.ctrl.mesh_type)){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	exit(-1);
	//printf("HEAT problem parameter overridden.\n");
	//strcpy(pdv_heat_problem.ctrl.mesh_type, pdv_ns_supg_problem.ctrl.mesh_type);
  }

  fprintf(Interactive_output, "\tmesh_file_in:\t\t\t\t%s\n",
	  pdv_ns_supg_problem.ctrl.mesh_filename);
  if(0 !=  strcmp(pdv_heat_problem.ctrl.mesh_filename,
		 pdv_ns_supg_problem.ctrl.mesh_filename)){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	exit(-1);
	//printf("HEAT problem parameter overridden.\n");
	//strcpy(pdv_heat_problem.ctrl.mesh_filename,
	//	   pdv_ns_supg_problem.ctrl.mesh_filename);
  }

  fprintf(Interactive_output, "\n\tns_supg field_file_in:\t\t\t%s\n",
	  pdv_ns_supg_problem.ctrl.field_filename);
  fprintf(Interactive_output, "\theat field_file_in:\t\t\t%s\n",
	  pdv_heat_problem.ctrl.field_filename);


  fprintf(Interactive_output, "\tns_supg mesh_file_out:\t\t\t%s\n",
	  pdv_ns_supg_problem.ctrl.mesh_dmp_filepattern);
  if(0 !=  strcmp(pdv_heat_problem.ctrl.mesh_dmp_filepattern,
		 pdv_ns_supg_problem.ctrl.mesh_dmp_filepattern)){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	exit(-1);
	//printf("HEAT problem parameter overridden.\n");
	//strcpy(pdv_heat_problem.ctrl.mesh_dmp_filepattern,
	//	   pdv_ns_supg_problem.ctrl.mesh_dmp_filepattern);
  }
  fprintf(Interactive_output, "\tns_supg field_file_out:\t\t\t%s\n",
	  pdv_ns_supg_problem.ctrl.field_dmp_filepattern);
  fprintf(Interactive_output, "\theat field_file_out:\t\t\t%s\n",
	  pdv_heat_problem.ctrl.field_dmp_filepattern);

  fprintf(Interactive_output,"\n\tns_supg penalty for Dirichlet BCs:\t%lf\n",
	  pdv_ns_supg_problem.ctrl.penalty);

  fprintf(Interactive_output,"\theat penalty for Dirichlet BCs:\t\t%lf\n",
	  pdv_heat_problem.ctrl.penalty);

  fprintf(Interactive_output, "\tns_supg gravity:\t\t\t%lf %lf %lf\n",
	  pdv_ns_supg_problem.ctrl.gravity_field[0],
	  pdv_ns_supg_problem.ctrl.gravity_field[1],
	  pdv_ns_supg_problem.ctrl.gravity_field[2]);

  /****************************************/
  /* initialization of material data      */
  /****************************************/
  fprintf(Interactive_output, "\n\tmaterials_file:\t\t\t\t%s\n",
	  pdv_ns_supg_problem.ctrl.material_filename);
  if(0 !=  strcmp(pdv_heat_problem.ctrl.material_filename,
		  pdv_ns_supg_problem.ctrl.material_filename)){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	exit(-1);
	//printf("HEAT problem parameter overridden.\n");
	//strcpy(pdv_heat_problem.ctrl.material_filename,
	//	   pdv_ns_supg_problem.ctrl.material_filename);
  }

  if(strlen(pdv_ns_supg_problem.ctrl.material_filename)==0){

	// e.g. for non-dimensional form of NS equations material files are not used
	pdv_ns_supg_problem.ctrl.ref_temperature = -1.0;

  }
  else{

	iaux = pdr_ns_supg_material_read(pdv_ns_supg_problem.ctrl.work_dir,
					 pdv_ns_supg_problem.ctrl.material_filename,
					 pdv_ns_supg_problem.ctrl.interactive_output);
	if (iaux == EXIT_FAILURE) exit(-1);

//    fprintf(Interactive_output, "\tnumber of materials:\t\t\t%d\n",
//	    pdv_ns_supg_problem.materials.materials_num);

//    fprintf(Interactive_output, "\nNS_SUPG material database read\n");

  }

  if(pdv_ns_supg_problem.ctrl.ref_temperature <= 0.0){

	// reference temperature used to indicate whether material database is used
	fprintf(Interactive_output, "\tns_supg material data not temperature dependent\n");
	fprintf(Interactive_output, "\n\tns_supg density:\t\t%lf\n",
		pdv_ns_supg_problem.ctrl.density);
	fprintf(Interactive_output, "\tns_supg dynamic viscosity:\t\t%lf\n",
		pdv_ns_supg_problem.ctrl.dynamic_viscosity);

  }
  else{

	fprintf(Interactive_output, "\tns_supg density and dynamic viscosity from material database\n");
	pdv_ns_supg_problem.ctrl.density = -1.0;
	pdv_ns_supg_problem.ctrl.dynamic_viscosity = -1.0;

	fprintf(Interactive_output, "\tns_supg reference_temperature:\t\t%lf\n",
		pdv_ns_supg_problem.ctrl.ref_temperature);

  }

  fprintf(Interactive_output, "\tns_supg reference_length:\t\t%lf\n",
	  pdv_ns_supg_problem.ctrl.ref_length);
  fprintf(Interactive_output, "\tns_supg reference_velocity:\t\t%lf\n",
	  pdv_ns_supg_problem.ctrl.ref_velocity);

  if(pdv_ns_supg_problem.ctrl.ref_temperature <= 0.0){
	fprintf(Interactive_output, "\tns_supg Reynolds number:\t\t%lf\n",
	  pdv_ns_supg_problem.ctrl.ref_length*pdv_ns_supg_problem.ctrl.ref_velocity*pdv_ns_supg_problem.ctrl.density/pdv_ns_supg_problem.ctrl.dynamic_viscosity);
  }

  /****************************************/
  /* initialization of heat material data */
  /****************************************/
  if(strlen(pdv_heat_problem.ctrl.material_filename)==0){

	// e.g. for non-dimensional form of heat equation material files are not used
	pdv_heat_problem.ctrl.ref_temperature = -1.0;

  }
  else{

	iaux = pdr_heat_material_read(pdv_heat_problem.ctrl.work_dir,
				  pdv_heat_problem.ctrl.material_filename,
				  pdv_heat_problem.ctrl.interactive_output);
	if (iaux == EXIT_FAILURE) exit(-1);

	fprintf(Interactive_output, "\nHEAT materials database read\n");

  }

  // ns_supg and heat can have different reference temperatures!!!
  // reference temperatures are used to indicate whether material data are constant
  if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){

	fprintf(Interactive_output, "\theat material data not temperature dependent\n");
	fprintf(Interactive_output, "\n\theat thermal_conductivity:\t\t%lf\n",
		pdv_heat_problem.ctrl.thermal_conductivity);
	fprintf(Interactive_output, "\n\theat density:\t\t%lf\n",
		pdv_heat_problem.ctrl.density);

	if(pdv_heat_problem.ctrl.density != pdv_ns_supg_problem.ctrl.density){
	  printf("Inconsistency between ns_supg and heat for density! Exiting.\n");
	  exit(-1);
	}

	fprintf(Interactive_output, "\n\theat specific_heat:\t\t%lf\n",
		pdv_heat_problem.ctrl.specific_heat);

  }
  else{

	fprintf(Interactive_output, "\theat material data from material database\n");
	pdv_heat_problem.ctrl.thermal_conductivity = -1.0;
	pdv_heat_problem.ctrl.density = -1.0;
	pdv_heat_problem.ctrl.specific_heat = -1.0;

	fprintf(Interactive_output, "\theat reference_temperature:\t\t%lf\n",
		pdv_heat_problem.ctrl.ref_temperature);

  }

  fprintf(Interactive_output, "\theat ambient temperature:\t\t%lf\n",
	  pdv_heat_problem.ctrl.ambient_temperature);


  fprintf(Interactive_output, "\nTIME INTEGRATION PARAMETERS:\n");

  fprintf(Interactive_output, "\ttype:\t\t\t\t\t%d\n",
	  pdv_ns_supg_problem.time.type);
  if(pdv_heat_problem.time.type != pdv_ns_supg_problem.time.type){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.type = pdv_ns_supg_problem.time.type;
  }
  fprintf(Interactive_output, "\timplicitness parameter:\t\t\t%lf\n\n",
	  pdv_ns_supg_problem.time.alpha);
  if(pdv_heat_problem.time.alpha != pdv_ns_supg_problem.time.alpha){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.alpha = pdv_ns_supg_problem.time.alpha;
  }

  fprintf(Interactive_output, "\tcurrent timestep:\t\t\t%d\n",
	  pdv_ns_supg_problem.time.cur_step);
  if(pdv_heat_problem.time.cur_step != pdv_ns_supg_problem.time.cur_step){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.cur_step = pdv_ns_supg_problem.time.cur_step;
  }
  fprintf(Interactive_output, "\tcurrent time:\t\t\t\t%lf\n",
	  pdv_ns_supg_problem.time.cur_time);
  if(pdv_heat_problem.time.cur_time != pdv_ns_supg_problem.time.cur_time){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.cur_time = pdv_ns_supg_problem.time.cur_time;
  }
  fprintf(Interactive_output, "\n\tcurrent timestep_length:\t\t%lf\n",
	  pdv_ns_supg_problem.time.cur_dtime);
  if(pdv_heat_problem.time.cur_dtime != pdv_ns_supg_problem.time.cur_dtime){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.cur_dtime = pdv_ns_supg_problem.time.cur_dtime;
  }
  fprintf(Interactive_output, "\tprevious timestep_length:\t\t%lf\n\n",
	  pdv_ns_supg_problem.time.prev_dtime);
  if(pdv_heat_problem.time.prev_dtime != pdv_ns_supg_problem.time.prev_dtime){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.prev_dtime = pdv_ns_supg_problem.time.prev_dtime;
  }

  fprintf(Interactive_output, "\tfinal time:\t\t\t\t%lf\n",
	  pdv_ns_supg_problem.time.final_time);
  if(pdv_heat_problem.time.final_time != pdv_ns_supg_problem.time.final_time){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.final_time = pdv_ns_supg_problem.time.final_time;
  }
  fprintf(Interactive_output, "\tfinal timestep:\t\t\t\t%d\n",
	  pdv_ns_supg_problem.time.final_step);
  if(pdv_heat_problem.time.final_step != pdv_ns_supg_problem.time.final_step){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.final_step = pdv_ns_supg_problem.time.final_step;
  }

  fprintf(Interactive_output, "\n\tns_supg convergence criterion type:\t%d\n",
	  pdv_ns_supg_problem.time.conv_type);
  fprintf(Interactive_output, "\tns_supg error tolerance (n-epsilon):\t%lf\n",
	  pdv_ns_supg_problem.time.conv_meas);
  fprintf(Interactive_output, "\n\theat convergence criterion type:\t%d\n",
	  pdv_heat_problem.time.conv_type);
  fprintf(Interactive_output, "\theat error tolerance (n-epsilon):\t%lf\n\n",
	  pdv_heat_problem.time.conv_meas);

  fprintf(Interactive_output, "\tmonitoring level:\t\t\t%d\n\n",
	  pdv_ns_supg_problem.time.monitor);
  if(pdv_heat_problem.time.monitor != pdv_ns_supg_problem.time.monitor){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.monitor = pdv_ns_supg_problem.time.monitor;
  }

  fprintf(Interactive_output, "\tgraph_dump_intv:\t\t\t%d\n",
	  pdv_ns_supg_problem.time.intv_graph);
  if(pdv_heat_problem.time.intv_graph != pdv_ns_supg_problem.time.intv_graph){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.intv_graph = pdv_ns_supg_problem.time.intv_graph;
  }
  fprintf(Interactive_output, "\tfull_dump_intv:\t\t\t\t%d\n\n",
	  pdv_ns_supg_problem.time.intv_dumpout);
  if(pdv_heat_problem.time.intv_dumpout!=pdv_ns_supg_problem.time.intv_dumpout){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.time.intv_dumpout = pdv_ns_supg_problem.time.intv_dumpout;
  }

  fprintf(Interactive_output, "\nPARAMETERS FOR TIME STEP LENGTH CONTROL BASED ON NON_LINEAR CONVERGENCE:\n");
  fprintf(Interactive_output, "\nALL PARAMETERS TAKEN FROM NS_SUPG (HEAT VALUES OVERRIDEN):\n");
  fprintf(Interactive_output, "\tcontrol (on/off switch):\t%d\n",
	  pdv_ns_supg_problem.time.time_step_length_nonl_control);
  pdv_heat_problem.time.time_step_length_nonl_control = pdv_ns_supg_problem.time.time_step_length_nonl_control;
  if(pdv_ns_supg_problem.time.time_step_length_nonl_control == 1){
	fprintf(Interactive_output, "\tnonl_iter_max:\t\t\t%d\n",
		pdv_ns_supg_problem.time.time_step_length_nonl_iter_max);
	pdv_heat_problem.time.time_step_length_nonl_iter_max = pdv_ns_supg_problem.time.time_step_length_nonl_iter_max;
	fprintf(Interactive_output, "\tnonl_iter_min:\t\t\t%d\n",
		pdv_ns_supg_problem.time.time_step_length_nonl_iter_min);
	pdv_heat_problem.time.time_step_length_nonl_iter_min = pdv_ns_supg_problem.time.time_step_length_nonl_iter_min;
	fprintf(Interactive_output, "\tincrease_mult:\t\t\t%lf\n",
		pdv_ns_supg_problem.time.time_step_length_nonl_iter_increase_mult);
	pdv_heat_problem.time.time_step_length_nonl_iter_increase_mult = pdv_ns_supg_problem.time.time_step_length_nonl_iter_increase_mult;
	fprintf(Interactive_output, "\tdecrease_mult:\t\t\t%lf\n",
		pdv_ns_supg_problem.time.time_step_length_nonl_iter_decrease_mult);
	pdv_heat_problem.time.time_step_length_nonl_iter_decrease_mult = pdv_ns_supg_problem.time.time_step_length_nonl_iter_decrease_mult;
	fprintf(Interactive_output, "\ttime_step_min:\t\t\t%lf\n",
		pdv_ns_supg_problem.time.time_step_length_min);
	pdv_heat_problem.time.time_step_length_min = pdv_ns_supg_problem.time.time_step_length_min;
	fprintf(Interactive_output, "\ttime_step_max:\t\t\t%lf\n",
		pdv_ns_supg_problem.time.time_step_length_max);
	pdv_heat_problem.time.time_step_length_max = pdv_ns_supg_problem.time.time_step_length_max;
  }

  fprintf(Interactive_output, "\nNONLINEAR SOLVER PARAMETERS:\n");

  fprintf(Interactive_output, "\ttype:\t\t\t\t\t%d\n",
	  pdv_ns_supg_problem.nonl.type);
  if(pdv_heat_problem.nonl.type != pdv_ns_supg_problem.nonl.type){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.nonl.type = pdv_ns_supg_problem.nonl.type;
  }

  fprintf(Interactive_output, "\tmax_nonl_iter:\t\t\t\t%d\n",
	  pdv_ns_supg_problem.nonl.max_iter);
  if(pdv_heat_problem.nonl.max_iter != pdv_ns_supg_problem.nonl.max_iter){
	printf("Inconsistency between ns_supg and heat for the above parameter!\n");
	printf("HEAT problem parameter overridden.\n");
	pdv_heat_problem.nonl.max_iter = pdv_ns_supg_problem.nonl.max_iter;
  }

  fprintf(Interactive_output, "\n\tns_supg convergence criterion type:\t%d\n",
	  pdv_ns_supg_problem.nonl.conv_type);
  fprintf(Interactive_output, "\tns_supg error tolerance (k-epsilon):\t%lf\n",
	  pdv_ns_supg_problem.nonl.conv_meas);
  fprintf(Interactive_output, "\tns_supg monitoring level:\t\t%d\n",
	  pdv_ns_supg_problem.nonl.monitor);

  fprintf(Interactive_output, "\n\theat convergence criterion type:\t%d\n",
	  pdv_heat_problem.nonl.conv_type);
  fprintf(Interactive_output, "\theat error tolerance (k-epsilon):\t%lf\n",
	  pdv_heat_problem.nonl.conv_meas);
  fprintf(Interactive_output, "\theat monitoring level:\t\t\t%d\n\n",
	  pdv_heat_problem.nonl.monitor);


  fprintf(Interactive_output, "\nLINEAR SOLVER PARAMETERS:\n");

  fprintf(Interactive_output, "\tns_supg solver type:\t\t\t%d\n",
	  pdv_ns_supg_problem.lins.type);

  fprintf(Interactive_output, "\tns_supg solver_file:\t\t\t%s\n",
	  pdv_ns_supg_problem.ctrl.solver_filename);

  if(pdv_ns_supg_problem.lins.type!=0
	 && (0 == strlen(pdv_heat_problem.ctrl.solver_filename))){
	fprintf(Interactive_output, "\n\tns_supg max_lins_iter:\t\t\t%d\n",
		pdv_ns_supg_problem.lins.max_iter);
	fprintf(Interactive_output, "\tns_supg convergence criterion type:\t%d\n",
		pdv_ns_supg_problem.lins.conv_type);
	fprintf(Interactive_output, "\tns_supg error tolerance:\t\t%.15lf\n",
		pdv_ns_supg_problem.lins.conv_meas);
		// pdr_lins_d_params(PDC_NS_SUPG_ID, 4));
	fprintf(Interactive_output, "\tns_supg monitoring level:\t\t%d\n",
		pdv_ns_supg_problem.lins.monitor);
  }


  fprintf(Interactive_output, "\n\theat solver type:\t\t\t%d\n",
	  pdv_heat_problem.lins.type);

  fprintf(Interactive_output, "\theat solver_file:\t\t\t%s\n",
	  pdv_heat_problem.ctrl.solver_filename);

  if(pdv_heat_problem.lins.type!=0
	 && (0 == strlen(pdv_heat_problem.ctrl.solver_filename))){
	fprintf(Interactive_output, "\n\theat max_lins_iter:\t\t\t%d\n",
		pdv_heat_problem.lins.max_iter);
	fprintf(Interactive_output, "\theat convergence criterion type:\t%d\n",
		pdv_heat_problem.lins.conv_type);
	fprintf(Interactive_output, "\theat error tolerance:\t\t\t%.15lf\n",
		pdv_heat_problem.lins.conv_meas);
			// pdr_lins_d_params(PDC_HEAT_ID, 4));
	fprintf(Interactive_output, "\theat monitoring level:\t\t\t%d\n\n",
		pdv_heat_problem.lins.monitor);
  }



  fprintf(Interactive_output, "\nADAPTATION PARAMETERS:\n");

  fprintf(Interactive_output, "\tns_supg adapt_type:\t\t\t%d\n",
	  pdv_ns_supg_problem.adpt.type);
  fprintf(Interactive_output, "\tns_supg adapt_interval:\t\t\t%d\n",
	  pdv_ns_supg_problem.adpt.interval);
  fprintf(Interactive_output, "\tns_supg adapt_eps:\t\t\t%lf\n",
	  pdv_ns_supg_problem.adpt.eps);
  fprintf(Interactive_output, "\tns_supg adapt_ratio:\t\t\t%lf\n\n",
	  pdv_ns_supg_problem.adpt.ratio);

  fprintf(Interactive_output, "\tns_supg monitoring level:\t\t%d\n\n",
	  pdv_ns_supg_problem.adpt.monitor);

  fprintf(Interactive_output, "\theat adapt_type:\t\t\t%d\n",
	  pdv_heat_problem.adpt.type);
  fprintf(Interactive_output, "\theat adapt_interval:\t\t\t%d\n",
	  pdv_heat_problem.adpt.interval);
  fprintf(Interactive_output, "\theat adapt_eps:\t\t\t\t%lf\n",
	  pdv_heat_problem.adpt.eps);
  fprintf(Interactive_output, "\theat adapt_ratio:\t\t\t%lf\n\n",
	  pdv_heat_problem.adpt.ratio);

  fprintf(Interactive_output, "\theat monitoring level:\t\t\t%d\n\n",
	  pdv_heat_problem.adpt.monitor);


  /****************************************/
  /* initialization of bc data            */
  /****************************************/

  iaux = pdr_ns_supg_bc_read(pdv_ns_supg_problem.ctrl.work_dir,
				 pdv_ns_supg_problem.ctrl.bc_filename,
				 pdv_ns_supg_problem.ctrl.interactive_output,
				 &pdv_ns_supg_problem.bc);
  if (iaux == EXIT_FAILURE) exit(-1);

  fprintf(Interactive_output, "\nns_supg boundary conditions configuration file:\t%s\n",
	  pdv_ns_supg_problem.ctrl.bc_filename);
  fprintf(Interactive_output, "\tnumber of pressure pins:\t\t%d\n",
	  pdr_ns_supg_get_pressure_pins_count(&pdv_ns_supg_problem.bc));
  fprintf(Interactive_output, "\tnumber of velocity pins:\t\t%d\n",
	  pdr_ns_supg_get_velocity_pins_count(&pdv_ns_supg_problem.bc));
  fprintf(Interactive_output, "\tnumber of BCs:\t\t\t\t%d\n\n",
	  pdr_ns_supg_get_bc_assign_count(&pdv_ns_supg_problem.bc));

  fprintf(Interactive_output, "NS_SUPG BC OK\n");

  iaux = pdr_heat_bc_read(pdv_heat_problem.ctrl.work_dir,
				 pdv_heat_problem.ctrl.bc_filename,
				 pdv_heat_problem.ctrl.interactive_output,
				 &pdv_heat_problem.bc);
  mf_check(iaux != EXIT_FAILURE, "Heat BC read failed!") ;

  fprintf(Interactive_output, "\nheat boundary conditions configuration file:\t%s\n",
	  pdv_heat_problem.ctrl.bc_filename);
  fprintf(Interactive_output, "\tnumber of BCs:\t\t\t\t%d\n\n",
	  pdr_heat_get_bc_assign_count(&pdv_heat_problem.bc));


  fprintf(Interactive_output, "HEAT BC OK\n");

  /****************************************/
  /* initialization of mesh data          */
  /****************************************/

  // mesh is read using ns_supg problem parameters (there is one mesh for
  // both problems)
  pdv_ns_supg_heat_current_problem_id = PDC_NS_SUPG_ID;
  mesh_id = utr_initialize_mesh( Interactive_output, Work_dir,
				 pdv_ns_supg_problem.ctrl.mesh_type[0],
				 pdv_ns_supg_problem.ctrl.mesh_filename);
  pdv_ns_supg_problem.ctrl.mesh_id = mesh_id;
  pdv_heat_problem.ctrl.mesh_id = mesh_id;

  fprintf(Interactive_output, "Read mesh file %s\n",
	  pdv_ns_supg_problem.ctrl.mesh_filename);


  utr_mesh_insert_new_bcnums(Work_dir,mesh_id);
  fprintf(stdout, "DADAa4\n");

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

#ifdef DEBUG
  {
  int currfa = 0;
  int fa_bnum;
  /* check if every boundary has been assigned boundary condtion */
  while (currfa = mmr_get_next_face_all(pdv_ns_supg_problem.ctrl.mesh_id,
					currfa)) {
	fa_bnum = mmr_fa_bc(pdv_ns_supg_problem.ctrl.mesh_id, currfa);
	//fprintf(Interactive_output, "BC NS_SUPG %d set for boundary %d\n",
	//	    pdr_ns_supg_get_bc_type(&pdv_ns_supg_problem.bc, fa_bnum), fa_bnum);
	if (fa_bnum > 0) {
	  if ((pdr_ns_supg_get_bc_type(&pdv_ns_supg_problem.bc, fa_bnum) == BC_NS_SUPG_NONE)) {
	fprintf(Interactive_output, "BC NS_SUPG not set for boundary:\t%d\n", fa_bnum);
	fprintf(Interactive_output, "Check bc config file - Exiting.\n");
	exit(-1);
	  }
	}
  }
  }
#endif

  mmr_set_max_gen(mesh_id, pdv_ns_supg_problem.adpt.maxgen);
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
  pdeg = pdv_ns_supg_problem.ctrl.pdeg;

  // first ns_supg field
  if (strcmp(pdv_ns_supg_problem.ctrl.field_filename,"z")==0){

	fprintf(Interactive_output,
		"\nInitializing ns_supg field with 0\n");
	//  - for standard continuous basis functions
	// 'z' - for zeroing field values

	pdv_ns_supg_problem.ctrl.field_id = utr_initialize_field(
							   Interactive_output, 's', 'z',
					   pdv_ns_supg_problem.ctrl.mesh_id,
				   pdv_ns_supg_problem.ctrl.nreq,
				   pdv_ns_supg_problem.ctrl.nr_sol, pdeg, NULL,NULL);


  }
  else if (strcmp(pdv_ns_supg_problem.ctrl.field_filename,"i")==0) {
	fprintf(Interactive_output,
		"\nInitializing ns_supg field with initial_condition function\n");
  // 's' - for standard continuous basis functions
  // 'i' - for initializing using function pdr_ns_supg_heat_initial_condition
	pdv_ns_supg_problem.ctrl.field_id = utr_initialize_field(
				   Interactive_output, 's', 'i',
					   pdv_ns_supg_problem.ctrl.mesh_id,
				   pdv_ns_supg_problem.ctrl.nreq,
				   pdv_ns_supg_problem.ctrl.nr_sol, pdeg,  NULL,
				   pdr_ns_supg_initial_condition);

  }
  else{

	sprintf(arg, "%s/%s", Work_dir, pdv_ns_supg_problem.ctrl.field_filename);
	char decompressed_file[255]={0};
	utr_io_decompress_file(NULL, arg, decompressed_file);
	strcpy(arg, decompressed_file);
	testforfile = fopen(arg, "r");
	if (testforfile != NULL) {
	  fclose(testforfile);
	  fprintf(Interactive_output, "Input field file %s.",
		  pdv_ns_supg_problem.ctrl.field_filename);
  // 's' - for standard continuous basis functions
  // 'i' - for initializing using function pdr_ns_supg_heat_initial_condition
	  pdv_ns_supg_problem.ctrl.field_id = utr_initialize_field(
							   Interactive_output, 's', 'r',
							   pdv_ns_supg_problem.ctrl.mesh_id,
							   pdv_ns_supg_problem.ctrl.nreq,
							   pdv_ns_supg_problem.ctrl.nr_sol, pdeg,  arg,NULL);

	} else {
	  fprintf(Interactive_output,
		  "Input field file %s not found - setting field to 0\n",
		  pdv_ns_supg_problem.ctrl.field_filename);
	  // 's' - for standard continuous basis functions
	  // 'z' - for zeroing field values
	  pdv_ns_supg_problem.ctrl.field_id = utr_initialize_field(
							   Interactive_output, 's', 'z',
					   pdv_ns_supg_problem.ctrl.mesh_id,
				   pdv_ns_supg_problem.ctrl.nreq,
				   pdv_ns_supg_problem.ctrl.nr_sol, pdeg,  NULL,NULL);
	}
  }
  apr_check_field(pdv_ns_supg_problem.ctrl.field_id);


  // second heat field


  if (strcmp(pdv_heat_problem.ctrl.field_filename,"z")==0){
	fprintf(Interactive_output,
		"\nInitializing heat field with 0\n");
	// 's' - for standard continuous basis functions
	// 'z' - for zeroing field values
	pdv_heat_problem.ctrl.field_id = utr_initialize_field(
							   Interactive_output, 's', 'z',
					   pdv_heat_problem.ctrl.mesh_id,
				   pdv_heat_problem.ctrl.nreq,
				   pdv_heat_problem.ctrl.nr_sol, pdeg, NULL,NULL);
  }
  else if (strcmp(pdv_heat_problem.ctrl.field_filename,"i")==0) {
	fprintf(Interactive_output,
		"\nInitializing heat field with initial_condition function\n");
  // 's' - for standard continuous basis functions
  // 'i' - for initializing using function pdr_heat_heat_initial_condition
	pdv_heat_problem.ctrl.field_id = utr_initialize_field(
							   Interactive_output, 's', 'i',
					   pdv_heat_problem.ctrl.mesh_id,
				   pdv_heat_problem.ctrl.nreq,
				   pdv_heat_problem.ctrl.nr_sol, pdeg,  NULL,
				   pdr_heat_initial_condition);
  }
  else{
	sprintf(arg, "%s/%s", Work_dir, pdv_heat_problem.ctrl.field_filename);
	char decompressed_file[255]={0};
	utr_io_decompress_file(NULL, arg, decompressed_file);
	strcpy(arg, decompressed_file);
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
							   pdv_heat_problem.ctrl.nr_sol, pdeg,  arg,NULL);
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
				   pdv_heat_problem.ctrl.nr_sol, pdeg, NULL,NULL);
	}
  }
  apr_check_field(pdv_heat_problem.ctrl.field_id);

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

return(0);
}

#ifdef __cplusplus
}
#endif
