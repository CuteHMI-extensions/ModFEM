#include <modfem/pd_heat/pdh_heat.h>

#include <modfem/aph_intf.h> /* USES */

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

// Normally global variables for a problem were defined in "pds_xxx_main.c" files, but they're moved experimentally here to avoid linking errors.

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

int pdv_heat_current_problem_id = PDC_HEAT_ID;
pdt_heat_problem pdv_heat_problem;

#ifdef PARALLEL
int pdv_exchange_table_index = -1;
#endif

// Similarly functions defined in main are moved here.

/*---------------------------------------------------------
pdr_err_indi - to return error indicator for an element
----------------------------------------------------------*/
double pdr_err_indi(		/* returns error indicator for an element */
  int Problem_id,	/* in: data structure to be used  */
  int Mode,	/* in: mode of operation */
  int El	/* in: element number */
	)
{

  if (Mode == PDC_ADAPT_ZZ) {

	return pdr_heat_err_indi_ZZ(Problem_id, El);

  } else {

	printf("Unknown error indicator in pdr_err_indi!\n");

  }

  return (0.0);
}

int pdr_ctrl_i_params(int Problem_id, int Num)
{
  pdt_heat_ctrls * ctrl_heat = & pdv_heat_problem.ctrl;

  if (Num == 2)
	  return ctrl_heat->mesh_id;
  else if (Num == 3)
	  return ctrl_heat->field_id;
  else if (Num == 4)
	  return ctrl_heat->nr_sol;
  else if (Num == 5)
	  return ctrl_heat->nreq;
  else if (Num == 6)
	  return ctrl_heat->solver_id;
  else
	  printf("Wrong parameter ID %d in pdr_ctrl_i_params for heat problem %d! Exiting.\n", Problem_id, Num);

  return -1;
}

void* pdr_get_problem_structure(int Problem_id)
{
  return (&pdv_heat_problem);
}

int pdr_adapt_i_params(int Problem_id, int Num)
{

  pdt_heat_adpts *adpts_heat = &pdv_heat_problem.adpt;

/*++++++++++++++++ executable statements ++++++++++++++++*/

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


  return (-1);
}


/*------------------------------------------------------------
pdr_adapt_d_params - to return parameters of adaptation
------------------------------------------------------------*/
double pdr_adapt_d_params(int Problem_id, int Num)
{


  pdt_heat_adpts *adpts_heat = &pdv_heat_problem.adpt;

/*++++++++++++++++ executable statements ++++++++++++++++*/

	if (Num == 5)
	  return (adpts_heat->eps);
	else if (Num == 6)
	  return (adpts_heat->ratio);
	else {
	  printf("Wrong parameter number in adapt_d_params!");
	  exit(1);
	}

  return (-1);
}

/*------------------------------------------------------------
  pdr_module_introduce - to return the problem module's name
------------------------------------------------------------*/
int pdr_module_introduce(
				  /* returns: >=0 - success code, <0 - error code */
  char* Problem_module_name /* out: the name of the problem module */
  )
{

  char* string = "HEAT";

  strcpy(Problem_module_name,string);

  return(1);
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

/* auxiliary variables */
  pdt_heat_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  // if only one problem allowed
  strcpy(Problem_name, pdv_heat_problem.ctrl.name);

  return(1);
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
  pdt_heat_linss *linss_heat = &pdv_heat_problem.lins;

/*++++++++++++++++ executable statements ++++++++++++++++*/


	if(Num==1) return(linss_heat->type);
	else if(Num==2) return(linss_heat->max_iter);
	else if(Num==3) return(linss_heat->conv_type);
	else if(Num==5) return(linss_heat->monitor);
	else {
	  printf("Wrong parameter number in lins_i_params!");
	  exit(-1);
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
  pdt_heat_linss *linss_heat = &pdv_heat_problem.lins;

/*++++++++++++++++ executable statements ++++++++++++++++*/


	if(Num==4) return(linss_heat->conv_meas);
	else {
	  printf("Wrong parameter number in lins_d_params!");
	  exit(-1);
	}


/* error condition - that point should not be reached */
  return(-1);
}

#ifdef __cplusplus
}
#endif
