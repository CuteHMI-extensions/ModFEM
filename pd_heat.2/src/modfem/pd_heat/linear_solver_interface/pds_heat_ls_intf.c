/************************************************************************
File pds_heat_ls_intf.c - interface between the problem dependent
		 module and linear solver modules (direct and iterative)

Contains definitions of routines:

Implementation of pdh_intf.h:

  pdr_get_list_ent - to return to the solver module :
						  1. the list of integration entities - entities
							 for which stiffness matrices and load vectors are
							 provided by the FEM code
						  2. the list of DOF entities - entities with which
							 there are dofs associated by the given approximation
  pdr_get_list_ent_coarse - the same as above but for COARSE level and
							given the corresponding lists from the fine level
  pdr_get_max_num_grid_levels - for limiting nr_levels in multigrid
							  based on mesh and field data
  pdr_create_assemble_stiff_mat - to create element stiffness matrices
								 and assemble them to the global SM
  pdr_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
								   to the global SM using assembly table
  pdr_assemble_local_stiff_mat - to assemble an element stiffness matrix
								   to the global SM
  pdr_comp_stiff_mat - to construct a stiffness matrix and a load vector for
					  some given mesh entity
  pdr_comp_stiff_mat_uncon - to construct UNCONSTRAINED stiffness matrix and
					  a load vector for some given mesh entity
  pdr_read_sol_dofs - to read a vector of dofs associated with a given
				   mesh entity from approximation field data structure
  pdr_write_sol_dofs - to write a vector of dofs associated with a given
				   mesh entity to approximation field data structure
  pdr_L2_proj_sol - to project solution between elements of different generations
  pdr_renum_coeff - to return a coefficient being a basis for renumbering
  pdr_get_ent_pdeg - to return the degree of approximation index
					  associated with a given mesh entity
  pdr_dof_ent_sons - to return info whether the entity is owned
					 and a list of dof entity sons for owned entity
  pdr_proj_sol_lev - to project solution between mesh levels
  pdr_vec_norm - to compute a norm of global vector in parallel
  pdr_sc_prod - to compute a scalar product of two global vectors
  pdr_create_exchange_tables - to create tables to exchange dofs
  pdr_exchange_dofs - to exchange dofs between processors

  pdr_select_el_coeff - to select coefficients returned to approximation
						routines for element integrals in weak formulation
		   (the procedure indicates which terms are non-zero in weak form)

  pdr_el_coeff - to return coefficients for internal integrals

Special procedure:
  pdr_heat_give_me_velocity_at_point - to provide the velocity and its
	gradient at a particular point given its local coordinates within an element
HEAT MODULE ASKS FOR IMPLEMENTATION - it has to be provided by procedures
defined in ls_intf directory of the problem module that uses heat as submodule

pdr_prepare_pde_coeff_streaming - to prepare an array with PDE coeffs for all elements
								  on L_elem_id list

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
	2011    - Przemyslaw Plaszewski
	2011    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)

*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>

/* problem dependent module interface */
#include <modfem/pdh_intf.h>		/* IMPLEMENTS */
#include <modfem/pdh_control_intf.h>		/* IMPLEMENTS */
/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>		/* USES */
/* interface for all approximation modules */
#include <modfem/aph_intf.h>		/* USES */
/* utilities - including simple time measurement library */
#include <modfem/uth_intf.h>		/* USES */
#include <modfem/uth_system.h> /* USES */

#ifdef PARALLEL
/* interface of parallel mesh manipulation modules */
#include "mmph_intf.h"		/* USES */
/* interface for all parallel approximation modules */
#include "apph_intf.h"		/* USES */
/* interface for parallel communication modules */
#include "pch_intf.h"		/* USES */
#endif

/* problem module's types and functions */
#include <modfem/pd_heat/pdh_heat.h> /* USES */
/* types and functions related to problem structures */
#include <modfem/pd_heat/pdh_heat_problem.h>
// bc and material header files are included in problem header files
/* weakform functions */
#include <modfem/pd_heat/pdh_heat_weakform.h>  	/* USES */

/* Constants */
// float versus double (MUST BE COMPATIBLE WITH KERNEL SWITCH!!!!)
// data type for integration and assembly
//#define SCALAR float
#define SCALAR double



/**************************************/
/* PDH_INTF.H PROCEDURES              */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */

/*------------------------------------------------------------
  pdr_get_list_ent - to return to the solver module :
						  1. the list of integration entities - entities
							 for which stiffness matrices and load vectors are
							 provided by the FEM code
						  2. the list of DOF entities - entities with which
							 there are dofs associated by the given approximation
------------------------------------------------------------*/
int pdr_get_list_ent(  /* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in:  problem (and solver) identification */
  int *Nr_int_ent,	/* out: number of integration entitites */
  int **List_int_ent_type,	/* out: list of types of integration entitites */
  int **List_int_ent_id,	/* out: list of IDs of integration entitites */
  int *Nr_dof_ent,	/* out: number of dof entities (entities with which there
			   are dofs associated by the given approximation) */
  int **List_dof_ent_type,	/* out: list of types of integration entitites */
  int **List_dof_ent_id,	/* out: list of IDs of integration entitites */
  int **List_dof_ent_nrdofs,	/* out: list of no of dofs for 'dof' entity */
  int *Nrdofs_glob,	/* out: global number of degrees of freedom (unknowns) */
  int *Max_dofs_per_dof_ent	/* out: maximal number of dofs per dof entity */
  )
{

  // default procedure for returning lists of integration and dof entities
  // for standard and discontinuous Galerkin approximations
  utr_get_list_ent(Problem_id, Nr_int_ent, List_int_ent_type, List_int_ent_id,
	 Nr_dof_ent, List_dof_ent_type, List_dof_ent_id, List_dof_ent_nrdofs,
		 Nrdofs_glob, Max_dofs_per_dof_ent);



  return (1);
}

/*------------------------------------------------------------
  pdr_get_list_ent_coarse - the same as above but for COARSE level and
							given the corresponding lists from the fine level
------------------------------------------------------------*/
int pdr_get_list_ent_coarse( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in:  problem (and solver) identification */
  int Nr_int_ent_fine,	/* in: number of integration entitites */
  int *List_int_ent_type_fine,	/* in: list of types of integration entitites */
  int *List_int_ent_id_fine,	/* in: list of IDs of integration entitites */
  int Nr_dof_ent_fine,	/* in: number of dof entities (entities with which there
			   are dofs associated by the given approximation) */
  int *List_dof_ent_type_fine,	/* in: list of types of integration entitites */
  int *List_dof_ent_id_fine,	/* in: list of IDs of integration entitites */
  int *List_dof_ent_nrdof_fine,	/* in: list of no of dofs for 'dof' entity */
  int Nrdof_glob_fine,	/* in: global number of degrees of freedom (unknowns) */
  int Max_dof_per_ent_fine,	/* in: maximal number of dofs per dof entity */
  int *Pdeg_coarse_p,	/* in: degree of approximation for coarse space */
  int *Nr_int_ent_p,	/* out: number of integration entitites */
  int **List_int_ent_type,	/* out: list of types of integration entitites */
  int **List_int_ent_id,	/* out: list of IDs of integration entitites */
  int *Nr_dof_ent_p,	/* out: number of dof entities (entities with which there
			   are dofs associated by the given approximation) */
  int **List_dof_ent_type,	/* out: list of types of integration entitites */
  int **List_dof_ent_id,	/* out: list of IDs of integration entitites */
  int **List_dof_ent_nrdofs,	/* out: list of no of dofs for 'dof' entity */
  int *Nrdof_glob,	/* out: global number of degrees of freedom (unknowns) */
  int *Max_dof_per_ent	/* out: maximal number of dofs per dof entity */
  )
{
  printf("pdr_get_list_ent_coarse NOT IMPLEMENTED for pdd_heat!");
  exit (-1);
}

/*------------------------------------------------------------
pdr_get_max_num_grid_levels - for limiting nr_levels in multigrid
							  based on mesh and field data
------------------------------------------------------------*/
int pdr_get_max_num_grid_levels(
  int Problem_id
)
{
  return(utr_get_max_num_grid_levels(Problem_id));
}


/*------------------------------------------------------------
 pdr_create_assemble_stiff_mat - to create element stiffness matrices
								 and assemble them to the global SM
------------------------------------------------------------*/
int pdr_create_assemble_stiff_mat(
  int Problem_id,
  int Level_id,
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int* Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
  int Nr_int_ent,
  int* L_int_ent_type,
  int* L_int_ent_id,
  int Nr_colors_elems,
  int* L_color_index_elems,
  int Nr_colors_faces,
  int* L_color_index_faces,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int* Pos_first_dof_int_ent,
  int* Local_to_global,
  int Max_dofs_int_ent
)
{

  utr_create_assemble_stiff_mat(Problem_id, Level_id, Comp_type, Pdeg_coarse_p,
				Nr_int_ent, L_int_ent_type, L_int_ent_id,
				Nr_colors_elems, L_color_index_elems,
				Nr_colors_faces, L_color_index_faces,
				Asse_pos_first_dof_int_ent, Assembly_table,
				Pos_first_dof_int_ent, Local_to_global,
				Max_dofs_int_ent);

  return(1);
}


/*------------------------------------------------------------
  pdr_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
								   to the global SM using assembly table
------------------------------------------------------------*/
int pdr_assemble_local_stiff_mat_with_table(
						 /* returns: >=0 - success code, <0 - error code */
  int Problem_id,        /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent,        /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* Assembly_table_int_ent, /* part of the global assembly table */
  int* Local_to_global_int_ent, /* part of the global local_to_global table */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs        /* in: flag to rewrite or sum up entries */
						 /*   'T' - true, rewrite entries when assembling */
						 /*   'F' - false, sum up entries when assembling */
)
{
  int i=6;
  int solver_id = pdr_ctrl_i_params(Problem_id,i);

  sir_assemble_local_stiff_mat_with_table(solver_id, Level_id, Comp_type, Nr_dof_ent,
				   Assembly_table_int_ent, Local_to_global_int_ent,
				   Stiff_mat, Rhs_vect, Rewr_dofs);

  return(1);

}



/*------------------------------------------------------------
  pdr_assemble_local_stiff_mat - to assemble an element stiffness matrix
								   to the global SM
------------------------------------------------------------*/
int pdr_assemble_local_stiff_mat(
						 /* returns: >=0 - success code, <0 - error code */
  int Problem_id,        /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent,        /* in: number of global dof blocks */
						 /*     associated with the local stiffness matrix */
  int* L_dof_ent_type,   /* in: list of dof blocks' IDs */
  int* L_dof_ent_id,     /* in: list of dof blocks' IDs */
  int* L_dof_ent_nrdofs, /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs        /* in: flag to rewrite or sum up entries */
						 /*   'T' - true, rewrite entries when assembling */
						 /*   'F' - false, sum up entries when assembling */
)
{
  int i=6;
  int solver_id = pdr_ctrl_i_params(Problem_id,i);

  sir_assemble_local_stiff_mat(solver_id, Level_id, Comp_type,
				   Nr_dof_ent, L_dof_ent_type,
				   L_dof_ent_id,L_dof_ent_nrdofs,
				   Stiff_mat, Rhs_vect, Rewr_dofs);

  return(1);

}

/*------------------------------------------------------------
  pdr_read_sol_dofs - to read a vector of dofs associated with a given
				   mesh entity from approximation field data structure
------------------------------------------------------------*/
int pdr_read_sol_dofs( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in: solver ID (used to identify the subproblem) */
  int Dof_ent_type, int Dof_ent_id, int Nrdof, double *Vect_dofs
						/* in: dofs to be written */
	)
{

  int i, field_id;
  int vect_id = 0;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/*kbw
  printf("in pdr_read_sol_dofs before apr_read_ent_dofs\n");
  printf("problem_id %d, Dof_ent_type %d, Dof_ent_id %d, nrdof %d\n",
  Problem_id, Dof_ent_type, Dof_ent_id,  Nrdof);
/*kew */

  i=3; field_id=pdr_ctrl_i_params(Problem_id,i);
  apr_read_ent_dofs(field_id, Dof_ent_type, Dof_ent_id,
			Nrdof, vect_id, Vect_dofs);

/*kbw
  printf("in pdr_read_sol_dofs after apr_read_ent_dofs\n");
  for(i=0;i<Nrdof;i++){
  printf("%10.6lf",Vect_dofs[i]);
  }
  printf("\n");
/*kew */


  return (1);
}

/*------------------------------------------------------------
  pdr_write_sol_dofs - to write a vector of dofs associated with a given
				   mesh entity to approximation field data structure
------------------------------------------------------------*/
int pdr_write_sol_dofs(	/* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in: solver ID (used to identify the subproblem) */
  int Dof_ent_type, int Dof_ent_id, int Nrdof,
  double *Vect_dofs	/* in: dofs to be written */
			)
{
  int vect_id = 0;

  int field_id=pdr_ctrl_i_params(Problem_id,3);
  apr_write_ent_dofs(field_id, Dof_ent_type, Dof_ent_id,
			 Nrdof, vect_id, Vect_dofs);

  return (1);
}

/*---------------------------------------------------------
  pdr_L2_proj_sol - to project solution between elements of different generations
----------------------------------------------------------*/
int pdr_L2_proj_sol(
  int Problem_id,	/* in: problem ID */
  int El,	/* in: element number */
  int *Pdeg,	/* in: element degree of approximation */
  double *Dofs,	/* out: workspace for degress of freedom of El */
  /*    NULL - write to  data structure */
  int *El_from,	/* in: list of elements to provide function */
  int *Pdeg_from,	/* in: degree of polynomial for each El_from */
  double *Dofs_from	/* in: Dofs of El_from or... */
			)
{

  int field_id, i;

  i=3; field_id=pdr_ctrl_i_params(Problem_id,i);
  i = -1; 	/* mode: -1 - projection from father to son */
  apr_L2_proj(field_id, i, El, Pdeg, Dofs, El_from, Pdeg_from, Dofs_from, NULL);

  return (0);
}

/*---------------------------------------------------------
pdr_renum_coeff - to return a coefficient being a basis for renumbering
----------------------------------------------------------*/
int pdr_renum_coeff(
  int Problem_id,	/* in: problem ID */
  int Ent_type,	/* in: type of mesh entity */
  int Ent_id,	/* in: mesh entity ID */
  double *Ren_coeff	/* out: renumbering coefficient */
			)
{
// no problem dependent renumbering
  *Ren_coeff = 1.0;
  return (1);
}


/*------------------------------------------------------------
  pdr_get_ent_pdeg - to return the degree of approximation index
					  associated with a given mesh entity
------------------------------------------------------------*/
int pdr_get_ent_pdeg(  /* returns: >0 - approximation index,
				   <0 - error code */
  int Problem_id,	/* in: approximation field ID  */
  int Ent_type,	/* in: type of mesh entity */
  int Ent_id	/* in: mesh entity ID */
	)
{
  int i, field_id;
  i=3; field_id=pdr_ctrl_i_params(Problem_id,i);
  return (apr_get_ent_pdeg(field_id, Ent_type, Ent_id));

}

/*---------------------------------------------------------
  pdr_dof_ent_sons - to return info whether the entity is owned
					 and a list of dof entity sons for owned entity
---------------------------------------------------------*/
int pdr_dof_ent_sons( /* returns: info whether the entity is owned */
  int Problem_id,     /* in: problem ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id,        /* in: mesh entity ID */
  int *Ent_sons     /* out: list of dof entity sons (for owned entities) */
					 /* 	Ent_sons[0] - number of sons */
  )
{
  return( utr_dof_ent_sons( Problem_id, Ent_type, Ent_id, Ent_sons) );
}


/*---------------------------------------------------------
  pdr_proj_sol_lev - to L2 project solution dofs between mesh levels
---------------------------------------------------------*/
int pdr_proj_sol_lev(		/* returns: >=0 - success; <0 - error code */
  int Problem_id,	/* in: problem ID */
  int Solver_id,	/* in: solver data structure to be used */
  int Ilev_from,	/* in: level number to project from */
  double *Vec_from,	/* in: vector of values to project */
  int Ilev_to,	/* in: level number to project to */
  double *Vec_to	/* out: vector of projected values */
  )
{

  printf("pdr_proj_sol_lev NOT IMPLEMENTED for pdd_heat!");
  exit (-1);
}

/*---------------------------------------------------------
  pdr_vec_norm - to compute a norm of global vector (in parallel)
---------------------------------------------------------*/
double pdr_vec_norm(		/* returns: L2 norm of global Vector */
  int Problem_id,	/* in: problem ID */
  int Solver_id,	/* in: solver data structure to be used */
  int Level_id,	/* in: level number */
  int Nrdof,	/* in: number of vector components */
  double *Vector	/* in: local part of global Vector */
  )
{

  double vec_norm = 0.0;
  int i, field_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef PARALLEL
  // simplified setting, only one problem and one field
  vec_norm = appr_sol_vec_norm(pdv_exchange_table_index,
				   Level_id, Nrdof, Vector);
#endif

  return (vec_norm);
}


/*---------------------------------------------------------
  pdr_sc_prod - to compute a scalar product of two global vectors
---------------------------------------------------------*/
double pdr_sc_prod(	/* retruns: scalar product of Vector1 and Vector2 */
  int Problem_id,	/* in: problem ID */
  int Solver_id,	/* in: solver data structure to be used */
  int Level_id,	/* in: level number */
  int Nrdof,	/* in: number of vector components */
  double *Vector1,	/* in: local part of global Vector */
  double *Vector2	/* in: local part of global Vector */
  )
{

  double sc_prod = 0.0;
  int i, field_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef PARALLEL
  // simplified setting, only one problem and one field
  sc_prod = appr_sol_sc_prod(pdv_exchange_table_index,
				 Level_id, Nrdof, Vector1, Vector2);
#endif

  return (sc_prod);
}

/*---------------------------------------------------------
  pdr_create_exchange_tables - to create tables to exchange dofs
---------------------------------------------------------*/
int pdr_create_exchange_tables(
				/* returns: >=0 -success code, <0 -error code */
  int Problem_id,	/* in: problem ID */
  int Solver_id,	/* in: solver data structure to be used */
  int Level_id,	/* in: level ID */
  int Nr_dof_ent,	/* in: number of DOF entities in the level */
  /* all four subsequent arrays are indexed by block IDs with 1(!!!) offset */
  int *L_dof_ent_type,	/* in: list of DOF entities associated with DOF blocks */
  int *L_dof_ent_id,	/* in: list of DOF entities associated with DOF blocks */
  int *L_bl_nrdof,	/* in: list of nrdofs for each dof block */
  int *L_bl_posg,	/* in: list of positions within the global */
  /*     vector of dofs for each dof block */
  int *L_elem_to_bl,	/* in: list of DOF blocks associated with DOF entities */
  int *L_face_to_bl,	/* in: list of DOF blocks associated with DOF entities */
  int *L_edge_to_bl,	/* in: list of DOF blocks associated with DOF entities */
  int *L_vert_to_bl	/* in: list of DOF blocks associated with DOF entities */
  )
{

  int i, field_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef PARALLEL
  // simplified setting, only one problem and one field
  appr_create_exchange_tables(pdv_exchange_table_index,
				  Level_id, Nr_dof_ent, L_dof_ent_type,
				  L_dof_ent_id, L_bl_nrdof, L_bl_posg,
				  L_elem_to_bl, L_face_to_bl,
				  L_edge_to_bl, L_vert_to_bl);
#endif

  return (0);

}

/*---------------------------------------------------------
  pdr_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
int pdr_exchange_dofs(
  int Problem_id,	/* in: problem ID */
  int Solver_id,	/* in: solver data structure to be used */
  int Level_id,	/* in: level number */
  double *Vec_dofs	/* in: vector of dofs to be exchanged */
  )
{

  int i, field_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef PARALLEL
  // simplified setting, only one problem and one field
  appr_exchange_dofs(pdv_exchange_table_index, Level_id, Vec_dofs);
#endif

  return (0);
}


/*------------------------------------------------------------
  pdr_select_el_coeff_vect - to select coefficients returned to approximation
						routines for element integrals in weak formulation
		   (the procedure indicates which terms are non-zero in weak form)
------------------------------------------------------------*/
int pdr_select_el_coeff_vect( // returns success indicator
  int Problem_id,
  int *Coeff_vect_ind	/* out: coefficient indicator */
				  )
{

  pdr_heat_select_el_coeff_vect(Problem_id, Coeff_vect_ind);

  return(1);

}

/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
/*------------------------------------------------------------
  pdr_select_el_coeff - to select coefficients returned to approximation
						routines for element integrals in weak formulation
		   (the procedure indicates which terms are non-zero in weak form)
------------------------------------------------------------*/
double *pdr_select_el_coeff( /* returns: pointer !=NULL to indicate selection */
  int Problem_id,
  double **Mval,	/* out: mass matrix coefficient */
  double **Axx,double **Axy,double **Axz, /* out:diffusion coefficients, e.g.*/
  double **Ayx,double **Ayy,double **Ayz, /* Axy denotes scalar or matrix */
  double **Azx,double **Azy,double **Azz, /* related to terms with dv/dx*du/dy */
  /* second order derivatives in weak formulation (scalar for scalar problems */
  /* matrix for vector problems) */
  double **Bx,double **By,double **Bz,	/* out: convection coefficients */
  /* Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
  double **Tx,double **Ty,double **Tz,	/* out: convection coefficients */
  /* Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
  double **Cval,/* out: reaction coefficients - for terms without derivatives */
  /*  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
  double **Lval,/* out: rhs coefficient for time term, Lval denotes scalar */
  /* or matrix corresponding to time derivative - similar as mass matrix but  */
  /* with known solution at the previous time step (usually denoted by u_n) */
  double **Qx,/* out: rhs coefficients for terms with derivatives */
  double **Qy,/* Qy denotes scalar or matrix corresponding to terms with dv/dy */
  double **Qz,/* derivatives in weak formulation */
  double **Sval	/* out: rhs coefficients without derivatives (source terms) */
  )
{


  double *select_coeff=NULL;

  select_coeff=pdr_heat_select_el_coeff(Problem_id, Mval,
					Axx,Axy,Axz,Ayx,Ayy,Ayz,Azx,Azy,Azz,
					Bx,By,Bz,Tx,Ty,Tz,Cval,Lval,Qx,Qy,Qz,Sval);

  return(select_coeff);

}


/*------------------------------------------------------------
pdr_el_coeff - to return coefficients for internal integrals
------------------------------------------------------------*/
int pdr_el_coeff(
  int Problem_id,
  int Elem,	/* in: element number */
  int Mat_num,	/* in: material number */
  double Hsize,	/* in: size of an element */
  int Pdeg,	/* in: local degree of polynomial */
  double *X_loc,      /* in: local coordinates of point within element */
  double *Base_phi,   /* in: basis functions */
  double *Base_dphix, /* in: x-derivatives of basis functions */
  double *Base_dphiy, /* in: y-derivatives of basis functions */
  double *Base_dphiz, /* in: z-derivatives of basis functions */
  double *Xcoor,	/* in: global coordinates of a point */
  double* Uk_val, 	/* in: computed solution from previous iteration */
  double* Uk_x, 	/* in: gradient of computed solution Uk_val */
  double* Uk_y,   	/* in: gradient of computed solution Uk_val */
  double* Uk_z,   	/* in: gradient of computed solution Uk_val */
  double* Un_val, 	/* in: computed solution from previous time step */
  double* Un_x, 	/* in: gradient of computed solution Un_val */
  double* Un_y,   	/* in: gradient of computed solution Un_val */
  double* Un_z,   	/* in: gradient of computed solution Un_val */
  double* Mval,	/* out: mass matrix coefficient */
  double *Axx, double *Axy, double *Axz,  /* out:diffusion coefficients */
  double *Ayx, double *Ayy, double *Ayz,  /* e.g. Axy denotes scalar or matrix */
  double *Azx, double *Azy, double *Azz,  /* related to terms with dv/dx*du/dy */
  /* second order derivatives in weak formulation (scalar for scalar problems */
  /* matrix for vector problems) */
  double *Bx, double *By, double *Bz,	/* out: convection coefficients */
  /* Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
  double *Tx, double *Ty, double *Tz,	/* out: convection coefficients */
  /* Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
  double *Cval,	/* out: reaction coefficients - for terms without derivatives */
  /*  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
  double *Lval,	/* out: rhs coefficient for time term, Lval denotes scalar */
  /* or matrix corresponding to time derivative - similar as mass matrix but  */
  /* with known solution at the previous time step (usually denoted by u_n) */
  double *Qx, /* out: rhs coefficients for terms with derivatives */
  double *Qy, /* Qy denotes scalar or matrix corresponding to terms with dv/dy */
  double *Qz, /* derivatives in weak formulation */
  double *Sval	/* out: rhs coefficients without derivatives (source terms) */
  )
{

	/* get velocity at gauss point */
	double vel[3];

	char problem_name[100];
	pdr_problem_name(Problem_id, problem_name);
	if(strcmp(problem_name, "cube_whirl") == 0){

	// we can get velocity at global point (specifying Elem == -1)
	  int el_temp = -1;
	  pdr_heat_get_velocity_at_point(Problem_id, el_temp, Xcoor, Base_phi,
					 NULL, NULL, NULL,
					 //Base_dphix, Base_dphiy, Base_dphiz,
					 vel, NULL, NULL, NULL);

	}
	else {

	  // or we can get velocity at local point within element (Elem >= 0)
	  pdr_heat_get_velocity_at_point(Problem_id, Elem, X_loc, Base_phi,
					 NULL, NULL, NULL,
					 //Base_dphix, Base_dphiy, Base_dphiz,
					 vel, NULL, NULL, NULL);
	}

/*kbw
	printf("pdr_heat_get_velocity at point %lf, %lf, %lf : %lf, %lf, %lf \n",
	   Xcoor[0], Xcoor[1], Xcoor[2], vel[0], vel[1], vel[2]);
/*kew*/

// There are two possible ways to navigate through problem parameters
// 1. get problem structure and then get directly parameters
// (the method is FASTER but requires the knowledge of problem structure type)

	pdt_heat_problem *problem =
	  (pdt_heat_problem *)pdr_get_problem_structure(Problem_id);
	int field_id = problem->ctrl.field_id;
	int mesh_id = problem->ctrl.mesh_id;
	// heating-cooling problem
	//int field_dtdt_id = pdr_ctrl_i_params(PDC_HEAT_DTDT_ID, 3);
	//double sol_dofs_dtdt[APC_MAXELSD];	/* solution dofs */
	int num_eq;

	double delta_t = problem->time.cur_dtime;
	double implicit = problem->time.alpha;

/*kbw
  printf("delta_t %lf, implicitness alpha = %lf\n",delta_t, implicit);
/*kew*/
  //int nreq = problem->ctrl.nreq;

// 2. get particular parameters using interface functions
// (the method is slower but can be used in modules that do not know the type
//  of problem structure for heat problem)

  //int nreq = pdr_ctrl_i_params(Problem_id, 5);

  /* select the proper field */
  //field_id = pdr_ctrl_i_params(Problem_id, 3);
  //mesh_id = apr_get_mesh_id(field_id);
  //nreq =apr_get_nreq(field_id);

  // nreq substituted as constant to allow compilers for constants propagation
	int nreq = PDC_HEAT_NREQ;
#ifdef DEBUG
	if(nreq != apr_get_nreq(field_id)){
	  printf("wrong parameter HEAT_NREQ in pdr_el_coeff 1\n");
	  printf("%d != %d. Exiting !!!",nreq, apr_get_nreq(field_id));
	  exit(-1);
	}
	if(nreq != pdr_ctrl_i_params(Problem_id,5)){
	  printf("wrong parameter HEAT_NREQ in pdr_el_coeff 2\n");
	  printf("%d != %d. Exiting !!!",nreq, pdr_ctrl_i_params(Problem_id,5));
	  exit(-1);
	}
	if(nreq != problem->ctrl.nreq){
	  printf("wrong parameter HEAT_NREQ in pdr_el_coeff 3\n");
	  printf("%d != %d. Exiting !!!",nreq, problem->ctrl.nreq);
	  exit(-1);
	}
#endif

	/*! ----------------------------------------------------------------------! */
	/*! -------------------- MATERIAL DATA AT GAUSS POINT ----------------- --! */
	/*! ----------------------------------------------------------------------! */
	double ref_temperature = problem->ctrl.ref_temperature;
	double thermal_conductivity;
	double specific_heat;
	double density;

	// for heat problems with constant material parameters (ref_temperature <= 0)
	if(ref_temperature<=0){
	  thermal_conductivity = problem->ctrl.thermal_conductivity;
	  density = problem->ctrl.density;
	  specific_heat = problem->ctrl.specific_heat;
	}
	// for heat problems with material parameters temperature dependent (ref_temperature > 0)
	else{
	  utt_material_query_params qparams;
	  utt_material_query_result qresult;
	  int i;

	  //1.set query parameters (which material and temperature)
	  qparams.group_idx = Mat_num;	//query by material index ...
	  //qparams.material_idx = 1;	//query by material index ...

	  /* printf("\nAS: (pdr_el_coeff) mat = %d", Mat_num); */
	  qparams.name = "";	//... not by material name
	  double tk = Uk_val[0]; // temperature is the only unknown
	  qparams.temperature = tk;	//temperature from last iteration
	  qparams.cell_id = Elem;
	  for( i=0; i<3; i++ ){
	qparams.xg[i] = Xcoor[i];
	  }
	  //2.get query results
	  pdr_heat_material_query( &qparams, &qresult);
	  //3.set values to those obtained with query
	  thermal_conductivity = qresult.thermal_conductivity;
	  specific_heat = qresult.specific_heat;
	  density = qresult.density;
	  //double thermal_diffusivity = tconductivity / ( density * specific_heat );
	  //double texpansion = qresult.thermal_expansion_coefficient;

	}
	  /*kbw
  printf("in supg_heat: conductivity %lf, specific_heat %lf, density %lf\n",
	 thermal_conductivity,specific_heat,density );
/*kew*/

	/* get coefficients for heat weak formulation */
	/* pdr_heat_el_coeff(Problem_id, Elem, Mat_num, Hsize, Pdeg, X_loc, */
	/* 		      Base_phi, Base_dphix, Base_dphiy, Base_dphiz, */
	/* 		      Xcoor, Uk_val, Uk_x, Uk_y, Uk_z, Un_val, Un_x, Un_y, Un_z, */
	/* 		      Mval, Axx, Axy, Axz, Ayx, Ayy, Ayz, Azx, Azy, Azz, */
	/* 		      Bx, By, Bz, Tx, Ty, Tz, Cval, Lval, Qx, Qy, Qz, Sval,  */
	/* 		      vel, thermal_diffusivity, delta_t, implicit); */
	double daux = density*specific_heat;
	pdr_heat_el_coeff(Problem_id, Elem, Mat_num, Hsize, Pdeg, X_loc,
			  Base_phi, Base_dphix, Base_dphiy, Base_dphiz,
			  Xcoor, Uk_val, Uk_x, Uk_y, Uk_z, Un_val, Un_x, Un_y, Un_z,
			  Mval, Axx, Axy, Axz, Ayx, Ayy, Ayz, Azx, Azy, Azz,
			  Bx, By, Bz, Tx, Ty, Tz, Cval, Lval, Qx, Qy, Qz, Sval,
			  vel, thermal_conductivity, daux, delta_t, implicit);

  return (0);
}


/*------------------------------------------------------------
  pdr_comp_stiff_mat - to provide a solver with a stiffness matrix
					  and a load vector corresponding to the specified
					  mesh entity, together with information on how to
					  assemble entries into the global stiffness matrix
					  and the global load vector
------------------------------------------------------------*/
int pdr_comp_stiff_mat(	 /* returns: >=0 - success code, <0 - error code */
  int Problem_id,	/* in: approximation field ID  */
  int Int_ent_type,	/* in: unique identifier of the integration entity */
  int Int_ent_id,	/* in: unique identifier of the integration entity */
  int Comp_sm,	/* in: indicator for the scope of computations: */
  /*   PDC_NO_COMP  - do not compute anything */
  /*   PDC_COMP_SM - compute entries to stiff matrix only */
  /*   PDC_COMP_RHS - compute entries to rhs vector only */
  /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg_vec,	/* in: enforced degree of polynomial (if > 0 ) */
  int *Nr_dof_ent,	/* in: size of arrays, */
			/* out: number of mesh entities with which dofs and */
			/*      stiffness matrix blocks are associated */
  int *List_dof_ent_type,	/* out: list of types for 'dof' entities */
  int *List_dof_ent_id,	/* out: list of ids for 'dof' entities */
  int *List_dof_ent_nrdofs,	/* out: list of no of dofs for 'dof' entity */
  int *Nrdofs_loc,	/* in(optional): size of Stiff_mat and Rhs_vect */
		   /* out(optional): actual number of dofs per integration entity */
  double *Stiff_mat,	/* out(optional): stiffness matrix stored columnwise */
  double *Rhs_vect,	/* outpds_elast_ls_std_intf.c(optional): rhs vector */
  char *Rewr_dofs	/* out(optional): flag to rewrite or sum up entries */
			/*   'T' - true, rewrite entries when assembling */
			/*   'F' - false, sum up entries when assembling */
 )
{

  // use generic utility that calls pdr_comp_stiff_mat_uncon for
  // creating problem dependent unconstrained stiff_mat and load_vec
  utr_comp_stiff_mat(Problem_id, Int_ent_type, Int_ent_id, Comp_sm,
			 Pdeg_vec, Nr_dof_ent, List_dof_ent_type,
			 List_dof_ent_id, List_dof_ent_nrdofs,
			 Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);


  return (1);
}

/**-----------------------------------------------------------
  pdr_comp_stiff_mat_uncon - to construct UNCONSTRAINED stiffness matrix and
					  a load vector for some given mesh entity
------------------------------------------------------------*/
int pdr_comp_stiff_mat_uncon( /** returns: >=0 - success code, <0 - error code */
  int Problem_id,     /** in: approximation field ID  */
  int Int_ent_type,    /** in: type of the integration entity */
  int Int_ent_id,    /** in: unique identifier of the integration entity */
  int Comp_sm,       /** in: indicator for the scope of computations: */
					 /**   PDC_NO_COMP  - do not compute anything */
					 /**   PDC_COMP_SM - compute entries to stiff matrix only */
					 /**   PDC_COMP_RHS - compute entries to rhs vector only */
					 /**   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg_vec,        /** in: enforced degree of polynomial (if != NULL ) */
  int* Nrdofs_loc,        /** in(optional): size of Stiff_mat and Rhs_vect */
				/** out(optional): actual number of dofs per integration entity*/
  double* Stiff_mat,      /** out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect,       /** out(optional): rhs vector */
  char* Rewr_dofs        /** out(optional): flag to rewrite or sum up entries */
						 /**   'T' - true, rewrite entries when assembling */
						 /**   'F' - false, sum up entries when assembling */
  )
{

  int pdeg;
  /* element degree of approximation for linear prisms is a single number */
  if (Pdeg_vec == NULL)
	pdeg = 0;
  else
	pdeg = Pdeg_vec[0];

/*kbw
	 printf("in pdr_comp_stiff_mat: problem_id %d, int_ent: type %d, id %d, enforced pdeg %d\n",
	 Problem_id, Int_ent_type, Int_ent_id, pdeg);
/*kew */

  /* definitions of both functions  in pds_heat_weakform */
  if (Int_ent_type == PDC_ELEMENT) {

	pdr_heat_comp_el_stiff_mat(Problem_id, Int_ent_id, Comp_sm, pdeg,
			  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  } else if (Int_ent_type == PDC_FACE) {

	pdr_heat_comp_fa_stiff_mat(Problem_id, Int_ent_id, Comp_sm, pdeg,
				   Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  } else {
	printf("ERROR: Wrong integration entity type in pdr_comp_stiff_mat!\n");
	exit(-1);
  }

  return(0);
}



/*------------------------------------------------------------
  pdr_heat_give_me_velocity_at_point - to provide the velocity and its
	gradient at a particular point given its local coordinates within an element
HEAT MODULE ASKS FOR IMPLEMENTATION - it has to be provided by procedures
defined in ls_intf directory of the problem module that uses heat as submodule
------------------------------------------------------------*/
int pdr_heat_give_me_velocity_at_point(
  int Problem_id,
  int El_id, // element
  double *X_loc, // local coordinates of point
  double *Base_phi, // shape functions at point (if available - to speed up)
  double *Base_dphix, // derivatives of shape functions at point
  double *Base_dphiy, // derivatives of shape functions at point
  double *Base_dphiz, // derivatives of shape functions at point
  double *Velocity, // temperature
  double *DVel_dx, // x-derivative of temperature
  double *DVel_dy, // y-derivative of temperature
  double *DVel_dz // z-derivative of temperature
  )
{

  // we call procedure to specify velocity
  pdr_heat_get_velocity_at_point(Problem_id, El_id, X_loc, Base_phi,
					Base_dphix, Base_dphiy, Base_dphiz,
					Velocity, DVel_dx, DVel_dy, DVel_dz);

  return(0);
}

/*------------------------------------------------------------
pdr_prepare_pde_coeff_streaming - to prepare an array with PDE coeffs for all elements
								  on L_elem_id list
------------------------------------------------------------*/
int pdr_prepare_pde_coeff_streaming(
  int Problem_id,
  int Nr_elems,
  int* L_elem_id,
  int* Size_global_pde_data_p,
  int* Size_per_element_pde_data_p,
  int* Size_per_int_point_pde_data_p,
  SCALAR** El_pde_dat_host_p
					)
{
  int i,j,k;

/*kbw
	printf("\nEntering pdr_prepare_pde_coeff_streaming (Problem %d, Level %d, Nr_elems %d:\n",
	   Problem_id, Level_id,  Nr_elems);
	printf("num_dofs %d, all_el_pde_coeff_size %d\n",
	   * Size_SM_LV_p, * Size_global_pde_data_p);
	printf("one_el_pde_coeff_size %d, one_int_p_pde_coeff_size %d\n",
	   * Size_per_element_pde_data_p, * Size_per_int_point_pde_data_p);
	printf("el_pde_dat_host %lu\n",
	   * El_pde_dat_host_p);

/*kew*/

  // find out the module name
  char module_name[100];
  pdr_module_introduce(module_name);

  // get problem name
  char problem_name[100];
  pdr_problem_name(Problem_id, problem_name);

  //
  // PREPARE PROBLEM DEPENDENT DATA FOR INTEGRATION
  //

  int field_id = pdr_ctrl_i_params(Problem_id, 3);
  int mesh_id = apr_get_mesh_id(field_id);
  int nreq = apr_get_nreq(field_id);

  if(nreq>PDC_MAXEQ){
	printf("nreq (%d) > PDC_MAXEQ (%d) in pdr_create_assemble_stiff_mat_opencl_elem\n",
	   nreq, PDC_MAXEQ);
	printf("Exiting!\n"); exit(-1);
  }

  // get the active PDE coefficient matrices
  /* pde coefficients */
  static int coeff_ind = 0;
  int coeff_ind_vect[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz,
  // 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
  // 14 - tx, 15 - ty, 16 - tz, 17 - cval
  // 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval

  double axx[PDC_MAXEQ*PDC_MAXEQ];
  double axy[PDC_MAXEQ*PDC_MAXEQ];
  double axz[PDC_MAXEQ*PDC_MAXEQ];
  double ayx[PDC_MAXEQ*PDC_MAXEQ];
  double ayy[PDC_MAXEQ*PDC_MAXEQ];
  double ayz[PDC_MAXEQ*PDC_MAXEQ];
  double azx[PDC_MAXEQ*PDC_MAXEQ];
  double azy[PDC_MAXEQ*PDC_MAXEQ];
  double azz[PDC_MAXEQ*PDC_MAXEQ];
  double bx[PDC_MAXEQ*PDC_MAXEQ];
  double by[PDC_MAXEQ*PDC_MAXEQ];
  double bz[PDC_MAXEQ*PDC_MAXEQ];
  double tx[PDC_MAXEQ*PDC_MAXEQ];
  double ty[PDC_MAXEQ*PDC_MAXEQ];
  double tz[PDC_MAXEQ*PDC_MAXEQ];
  double cval[PDC_MAXEQ*PDC_MAXEQ];
  double mval[PDC_MAXEQ*PDC_MAXEQ];
  double qx[PDC_MAXEQ];
  double qy[PDC_MAXEQ];
  double qz[PDC_MAXEQ];
  double sval[PDC_MAXEQ];
  double lval[PDC_MAXEQ];

  pdr_select_el_coeff_vect(Problem_id, coeff_ind_vect);

  // there are two choices:
  // 1. consider all terms (16 arrays and 4 vectors)
  // 2. use coeff_vector_indicator to select which terms are non-zero
  int coeff_array_indicator[16];
  int nr_coeff_arrays = 0;
  for(i=0; i<16; i++) {
	coeff_array_indicator[i]=coeff_ind_vect[i+2];
	if(coeff_ind_vect[i+2]==1) nr_coeff_arrays++;
  }
  // in create_assemble we combine mval (1) with cval (17)
  // we use this for vectorization of kernels
  if(coeff_ind_vect[1]==1) {
	if(coeff_array_indicator[15]!=1) nr_coeff_arrays++;
	coeff_array_indicator[15]=1;
  }
  int coeff_vector_indicator[4];
  int nr_coeff_vectors = 0;
  for(i=0; i<4; i++) {
	coeff_vector_indicator[i]=coeff_ind_vect[i+19];
	if(coeff_ind_vect[i+19]==1) nr_coeff_vectors++;
  }
  // in create_assemble we combine lval (18) with sval (22)
  // we use this for vectorization of kernels
  if(coeff_ind_vect[18]==1) {
	if(coeff_vector_indicator[3]!=1) nr_coeff_vectors++;
	coeff_vector_indicator[3]=1;
  }



/*kbw*/
  printf("Problem ID: %d, problem_name %s\n", Problem_id, problem_name);
  printf("\nthe number of coefficients arrays %d, indicator = \n",nr_coeff_arrays);
  for(i=0;i<16;i++){
	printf("%2d",coeff_array_indicator[i]);
  }
  printf("\n");
  printf("the number of coefficients vectors %d, indicator = \n",nr_coeff_vectors);
  for(i=0;i<4;i++){
	printf("%2d",coeff_vector_indicator[i]);
  }
  printf("\n");
/*kew*/

  //
  // BASED ON PDEG COMPUTE EXECUTION CHARACTERISTICS
  //

  // choose an example element for a given pdeg and color
  int el_id = L_elem_id[0];

  int pdeg = apr_get_el_pdeg(field_id, el_id, NULL);
  int num_shap = apr_get_el_pdeg_numshap(field_id, el_id, &pdeg);
  int num_dofs = num_shap * nreq;

  // get the size of quadrature data (we assume all elements are of the same type)
  int base = apr_get_base_type(field_id, L_elem_id[0]);
  int ngauss;            /* number of gauss points */
  double xg[300];   	 /* coordinates of gauss points in 3D */
  double wg[100];       /* gauss weights */
  apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);

  //int el_nodes[MMC_MAXELVNO+1];        /* list of nodes of El */
  int el_nodes[APC_MAX_GEO_DOFS+1];        /* list of nodes of El */
  int el_nodes_type[APC_MAX_GEO_DOFS+1];        /* list of nodes type of El */
  //mmr_el_node_coor(mesh_id, el_id, el_nodes, NULL);
  apr_get_el_geo_dofs(field_id,el_id,el_nodes,el_nodes_type,NULL);

  /* for geometrically (multi)linear elements number of geometrical  */
  /* degrees of freedom is equal to the number of vertices - classical FEM nodes */
  int num_geo_dofs = el_nodes[0];

  int pdeg_single = pdeg;
  if(pdeg>100) {
	pdeg_single = pdeg/100;
	if(pdeg != pdeg_single*100 + pdeg_single){
	  printf("wrong pdeg %d ( > 100 but not x0x )\n", pdeg);
	  exit(-1);
	}
  }

  /*kbw*/
  printf("problem and element characteristics: nreq %d, pdeg %d, num_shap %d, num_dofs %d!\n",
	 nreq, pdeg, num_shap, num_dofs);
  //kew*/


  // 4. PDE COEFFICIENTS

  // there are two choices:
  // 1. consider all terms (16 arrays and 4 vectors)
  // 2. use coeff_vector_indicator to select which terms are non-zero
  // this options have to be taken into account when rewriting coefficients returned
  // by problem dependent module to coeff array
  // HERE coeff_vector_indicator is used to calculate pde_coeff_size
  int pde_coeff_size = nr_coeff_arrays*nreq*nreq + nr_coeff_vectors*nreq;


  // for different types of problems solved, there is a multitude of options:
  // 1. no coefficients at all (e.g. Laplace)
  // 2. coefficients constant for all elements (elasticity and uniform material)
  // 3. coefficients constant for all integration points (elasticity)
  // 4. coefficients different for all elements and integration points
  //    a. for nonlinear problems with coefficients sent to GPU
  //    b. for multi-scale problems??? (data different for each integration point)
  // also possible
  // coefficients are sent to GPU but SM calculations involve not only coefficients
  // but some other data as well (previous solution, etc.)
  // then additional factors have to be taken into account (time for sending
  // solution or degrees of freedom for computing solution, time for computing
  // solution from degrees of freedom, time for computing SM entries given
  // coefficients and the other data, etc.)

  // FOR SIMPLE PROBLEMS WE SUGGEST TO USE:
  //      REG_JAC, REG_NOJAC, SHM_JAC or SHM_NOJAC - THESE DO NOT REQUIRE
  //      DIFFERENT COEFFICIENTS FOR EACH ELEMENT AND INTEGRATION POINT
  // FOR COMPLEX PROBLEMS WE SUGGEST TO USE:
  //      REG_GL_DER or SHM_GL_DER WITH READY MADE COEFFICIENTS SENT
  //      TO GPU (A LOT OF DATA!!!!!!)

  // different parameters to differentiate between different cases
  int all_el_pde_coeff_size = 0; // size for all elements
  int one_el_pde_coeff_size = 0; // size for one element and all integration point
  int one_int_p_pde_coeff_size = 0; // size for one element and one integration point

  // default - not practical: all coefficients at all integration points sent
  one_int_p_pde_coeff_size  = pde_coeff_size;

  // coefficients common to all elements
  all_el_pde_coeff_size = 2; // Delta_t, Implicitness_coeff
  one_el_pde_coeff_size = 1; // Element size, DOFs from previous iteration andDOFs from previous time step is already send
  one_int_p_pde_coeff_size = 5; // 3 velocity components,  thermal conductivity, Density_times_specific_heat,

  int el_pde_data_size = all_el_pde_coeff_size +
			  Nr_elems*(one_el_pde_coeff_size + ngauss*one_int_p_pde_coeff_size);

  // 5. COMPUTED STIFFNESS MATRIX AND LOAD VECTOR
  int one_el_stiff_mat_size = num_dofs*num_dofs;
  int one_el_load_vec_size = num_dofs;

/*kbw*/
  printf("\nData structure sizes:\n");
  printf("\tPDE coefficients for each element %d or each integration point %d\n",
	 one_el_pde_coeff_size, one_int_p_pde_coeff_size);
  printf("\tPDE coefficients for all elements %d and total %d\n",
	 all_el_pde_coeff_size, el_pde_data_size);
  printf("\tSM for each element %d\n", one_el_stiff_mat_size);
  printf("\tLoad vector for each element %d\n", one_el_load_vec_size);
/*kew*/



  // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
  // finally fill element input data
  int ielem;
  int packed_bytes = 0;
  int final_position = 0; // for testing packing procedure

  int el_pde_data_bytes = el_pde_data_size * sizeof(SCALAR);
  SCALAR* el_pde_data = (SCALAR *)malloc(el_pde_data_bytes);

  //#ifndef OPENCL_HSA
  memset(el_pde_data, 0, el_pde_data_bytes);
  //#endif

  int position_coeff = 0;

  pdt_heat_problem *problem =
	(pdt_heat_problem *)pdr_get_problem_structure(Problem_id);

  double delta_t = problem->time.cur_dtime; // current time-step length
  double implicitness_parameter = problem->time.alpha; // implicitnes parameter alpha

  /* jbw
  printf("[pdr_prepare_pde_coeff_streaming] ALL position_coeff = %d ; implicitness_parameter = %lf delta_t = %lf\n",
	 position_coeff,implicitness_parameter,delta_t);
  /* jbw */

  el_pde_data[position_coeff+0] = implicitness_parameter;
  el_pde_data[position_coeff+1] = delta_t;

  position_coeff += 2;
  packed_bytes += 2*sizeof(SCALAR);

  // LOOP OVER ELEMENTS
  for(ielem=0; ielem<Nr_elems; ielem++)
	{

	  // element ID
	  el_id = L_elem_id[ielem];

	  // data common to all integration points:
	  // element size + DOFs from previous time step;

	  int el_mate = mmr_el_groupID(mesh_id, el_id);

	  double geo_dofs[3*MMC_MAXELVNO];  /* coord of nodes of El */
	  double el_hsize = mmr_el_hsize(mesh_id, el_id, NULL,NULL,NULL);

	  /* geo_dofs - geometry degrees of freedom - used by apr_elem_calc_3D  */
	  mmr_el_node_coor(mesh_id, el_id, el_nodes, geo_dofs);


	  // checking whether this element has the same data as assumed for this color
	  assert( pdeg == apr_get_el_pdeg(field_id, el_id, NULL) );
	  assert( num_shap == apr_get_el_pdeg_numshap(field_id, el_id, &pdeg) );
	  assert( base == apr_get_base_type(field_id, el_id) );

	  // PDE coefficients

	  // problem data_structure
	  pdt_heat_problem *problem =
	(pdt_heat_problem *)pdr_get_problem_structure(Problem_id);

	  double base_phi[APC_MAXELVD];	/* basis functions */
	  double base_dphix[APC_MAXELVD];	/* x-derivatives of basis function */
	  double base_dphiy[APC_MAXELVD];	/* y-derivatives of basis function */
	  double base_dphiz[APC_MAXELVD];	/* y-derivatives of basis function */

	  double sol_dofs_k[APC_MAXELSD];	/* solution dofs, from previous iteration */
	  int sol_vec_id, tmp_ngauss;

	  sol_vec_id = 1; // solution from the previous iteration
	  apr_get_el_dofs(field_id, el_id, sol_vec_id, sol_dofs_k);

	  // data sent:

	  /* jbw
	  printf("[pdr_prepare_pde_coeff_streaming] ELEM ielem/Nr_elems = %d/%d ; position_coeff = %d ; el_hsize = %lf\n",ielem,Nr_elems,position_coeff,el_hsize);
	  /* jbw */

	  el_pde_data[position_coeff] = el_hsize;
	  position_coeff++;


	  /* ----------------------- CALCULATE ON GPU ----------------------------
	  pdeg = 0;
	  apr_set_quadr_3D(base,&pdeg,&tmp_ngauss,xg,wg);
	  pdeg = apr_get_el_pdeg(field_id, el_id, NULL);
	  // XCOOR V1
	  //xcoor_middle[0] = (geo_dofs[0+0]+geo_dofs[3+0]+geo_dofs[6+0]
	  //	       +geo_dofs[9+0]+geo_dofs[12+0]+geo_dofs[15+0])/6.0;
	  //xcoor_middle[1] = (geo_dofs[0+1]+geo_dofs[3+1]+geo_dofs[6+1]
	  //	       +geo_dofs[9+1]+geo_dofs[12+1]+geo_dofs[15+1])/6.0;
	  //xcoor_middle[2] = (geo_dofs[0+2]+geo_dofs[3+2]+geo_dofs[6+2]
	  //	       +geo_dofs[9+2]+geo_dofs[12+2]+geo_dofs[15+2])/6.0;
	  // XCOOR V2
	  xcoor_middle[0] = xg[3*i+0];
	  xcoor_middle[1] = xg[3*i+1];
	  xcoor_middle[2] = xg[3*i+2];
	  apr_elem_calc_3D(2, nreq, &pdeg, base,
			   xcoor_middle, geo_dofs, NULL,
			   base_phi,base_dphix,base_dphiy,base_dphiz,
			   NULL,NULL,NULL,NULL,NULL,NULL);

	  double velocity[3];
	  pdr_heat_get_velocity_at_point(Problem_id, el_id, xcoor_middle, base_phi,
					 NULL, NULL, NULL,
					 //Base_dphix, Base_dphiy, Base_dphiz,
					 velocity, NULL, NULL, NULL);
	  double m_k = 1.0/3.0; // should be changed for higher order elements
	  double norm_u = sqrt( velocity[0] * velocity[0] +
				velocity[1] * velocity[1] +
				velocity[2] * velocity[2] );
	  // h_k computations: FRANCA, TEZDUYAR
	  double h_k = 0.0;

	  if(norm_u < 1.0e-6){ // when there is no velocity field inside the element
	h_k = el_hsize;
	  }
	  else{ // take element size in the direction of velocity
	// this definition may lead to oscillations - change then to standard size
	int idofs;
	for (idofs = 0; idofs < num_shap; idofs++) {
	  h_k += fabs(velocity[0] * base_dphix[idofs] +
			  velocity[1] * base_dphiy[idofs] +
			  velocity[2] * base_dphiz[idofs]);
	}
	h_k = 2.0 * norm_u / h_k;
	  }

	  //el_pde_data[position_coeff] = h_k;
	  //position_coeff++;
	  /* ----------------------- CALCULATE ON GPU ---------------------------- */

	  // data for each integration point:
	  // 3 velocity components,  thermal conductivity, Density_times_specific_heat

	  int ki; double velocity[3];
	  for(ki=0;ki<ngauss;ki++)
	{

	  double xcoor[3];
	  double u_sol;
	  apr_elem_calc_3D(1, nreq, &pdeg, base,
			   &xg[3*ki], geo_dofs, sol_dofs_k,
			   base_phi,NULL,NULL,NULL,
			   xcoor,&u_sol,NULL,NULL,NULL,NULL);

	  /*! ----------------------------------------------------------------------! */
	  /*! -------------------- MATERIAL DATA AT GAUSS POINT ----------------- --! */
	  /*! ----------------------------------------------------------------------! */
	  double ref_temperature = problem->ctrl.ref_temperature;
	  double thermal_conductivity;
	  double specific_heat;
	  double density;

	  // for heat problems with constant material parameters (ref_temperature <= 0)
	  if(ref_temperature<=0)
		{
		  thermal_conductivity = problem->ctrl.thermal_conductivity;
		  density = problem->ctrl.density;
		  specific_heat = problem->ctrl.specific_heat;
		}
	  // for heat problems with material parameters temperature dependent (ref_temperature > 0)
	  else
		{
		  utt_material_query_params qparams;
		  utt_material_query_result qresult;
		  int i;

		  //1.set query parameters (which material and temperature)
		  qparams.group_idx = el_mate;	//query by material index ...
		  //qparams.material_idx = 1;	//query by material index ...

		  qparams.name = "";	//... not by material name
		  double tk = u_sol; // temperature is the only unknown
		  qparams.temperature = tk;	//temperature from last iteration
		  qparams.cell_id = el_id;
		  for( i=0; i<3; i++ )  qparams.xg[i] = xcoor[i];

		  //2.get query results
		  pdr_heat_material_query(&qparams, &qresult);

		  //3.set values to those obtained with query
		  thermal_conductivity = qresult.thermal_conductivity;
		  specific_heat = qresult.specific_heat;
		  density = qresult.density;
		  //double thermal_diffusivity = tconductivity / ( density * specific_heat );
		  //double texpansion = qresult.thermal_expansion_coefficient;
		}

	  double dens_spec_heat = density*specific_heat;

	  char problem_name[100];
	  pdr_problem_name(Problem_id, problem_name);
	  if(strcmp(problem_name, "cube_whirl") == 0){

		// we can get velocity at global point (specifying Elem == -1)
		int el_temp = -1;
		pdr_heat_get_velocity_at_point(Problem_id, el_temp, xcoor, NULL,
					   NULL, NULL, NULL,
					   //Base_dphix, Base_dphiy, Base_dphiz,
					   velocity, NULL, NULL, NULL);

	  }
	  else {

		// or we can get velocity at local point within element (Elem >= 0)
		pdr_heat_get_velocity_at_point(Problem_id, el_id, &xg[3*ki], base_phi,
					   NULL, NULL, NULL,
					   //Base_dphix, Base_dphiy, Base_dphiz,
					   velocity, NULL, NULL, NULL);
	  }

/*kbw
  printf("pdr_heat_get_velocity at point %lf, %lf, %lf : %lf, %lf, %lf \n",
		 xcoor[0], xcoor[1], xcoor[2], velocity[0], velocity[1], velocity[2]);
/*kew*/

	  /* jbw
	  printf("[pdr_prepare_pde_coeff_streaming] GAUSS ielem = %d ; ki/ngauss = %d/%d ;  position_coeff = %d velocity = <%lf, %lf, %lf> ; tc = %lf ; dsh = %lf\n",
		 ielem,ki,ngauss,position_coeff,velocity[0],velocity[1],velocity[2],thermal_conductivity,dens_spec_heat);
	  /* jbw */

	  el_pde_data[position_coeff] = velocity[0];
	  el_pde_data[position_coeff+1] = velocity[1];
	  el_pde_data[position_coeff+2] = velocity[2];
	  position_coeff+=3;

	  el_pde_data[position_coeff] = thermal_conductivity;
	  position_coeff++;

	  el_pde_data[position_coeff] = dens_spec_heat;
	  position_coeff++;
	}

	  // written data:
	  // h_k + dofs previous iteration + dofs previous timestep +
	  // ngauss * ( 3 componets of velocity + thermal_conductivity + dens_spec_heat )
	  packed_bytes += (1 + ngauss*(3+2))*sizeof(SCALAR);

	} // end loop over elements

  *Size_global_pde_data_p = all_el_pde_coeff_size;
  *Size_per_element_pde_data_p = one_el_pde_coeff_size;
  *Size_per_int_point_pde_data_p = one_int_p_pde_coeff_size;
  // total size = Size_global_pde_data + Nr_elems *
  //              (Size_per_element_pde_data + ngauss*Size_per_int_point_pde_data)
  *El_pde_dat_host_p = el_pde_data;

  return(0);
}
