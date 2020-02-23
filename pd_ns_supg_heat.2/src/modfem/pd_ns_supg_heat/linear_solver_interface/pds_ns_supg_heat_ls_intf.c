/************************************************************************
File pds_ns_supg_heat_ls_intf.c - interface between the problem dependent
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
  pdr_ns_supg_give_me_temperature_at_point - to provide the temperature and its
	gradient at a particular point given its local coordinates within an element
NS_SUPG MODULE ASKS FOR IMPLEMENTATION - it has to be provided by procedures
defined in ls_intf directory of the problem module that uses ns_supg as submodule

Two local procedures :

  pdr_ns_supg_heat_comp_el_stiff_mat - to construct a stiffness matrix and a load
							vector for an element (used by pdr_comp_stiff_mat)

  pdr_ns_supg_heat_comp_fa_stiff_mat - to construct a stiffness matrix and a load
						   vector for a face (used by pdr_comp_stiff_mat)

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
#include <modfem/sih_intf.h>

#include <modfem/uth_system.h>

#ifdef PARALLEL
/* interface of parallel mesh manipulation modules */
#include <modfem/mmph_intf.h>		/* USES */
/* interface for all parallel approximation modules */
#include <modfem/apph_intf.h>		/* USES */
/* interface for parallel communication modules */
#include <modfem/pch_intf.h>		/* USES */
#endif

/* problem module's types and functions */
#include <modfem/pd_ns_supg_heat/pdh_ns_supg_heat.h>	/* USES */
/* types and functions related to problem structures */
#include <modfem/pd_ns_supg/pdh_ns_supg_problem.h>
#include <modfem/pd_heat/pdh_heat_problem.h>
// bc and material header files are included in problem header files
/* weakform functions */
#include <modfem/pd_ns_supg/pdh_ns_supg_weakform.h>	/* USES */
#include <modfem/pd_heat/pdh_heat_weakform.h>  	/* USES */


#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------
  pdr_ns_supg_heat_comp_fa_stiff_mat - to construct a stiffness matrix and
						  a load vector for a face
------------------------------------------------------------*/
int pdr_ns_supg_heat_comp_fa_stiff_mat(
			/*returns: >=0 -success code, <0 -error code */
  int Problem_id,	/* in: approximation field ID  */
  int Fa_id,	/* in: unique identifier of the face */
  int Comp_sm,	/* in: indicator for the scope of computations: */
  /*   PDC_NO_COMP  - do not compute anything */
  /*   PDC_COMP_SM - compute entries to stiff matrix only */
  /*   PDC_COMP_RHS - compute entries to rhs vector only */
  /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int Pdeg_in,	/* in: enforced degree of polynomial (if > 0 ) */
  int *Nrdofs_loc,	/* in(optional): size of Stiff_mat and Rhs_vect */
  /* out(optional): actual number of dofs per integration entity */
  double *Stiff_mat,	/* out(optional): stiffness matrix stored columnwise */
  double *Rhs_vect,	/* out(optional): rhs vector */
  char *Rewr_dofs	/* out(optional): flag to rewrite or sum up entries */
  /*   'T' - true, rewrite entries when assembling */
  /*   'F' - false, sum up entries when assembling */
  );

/*------------------------------------------------------------
  pdr_ns_supg_heat_comp_el_stiff_mat - to construct a stiffness matrix and
						  a load vector for an element
------------------------------------------------------------*/
int pdr_ns_supg_heat_comp_el_stiff_mat(
			/*returns: >=0 - success code, <0 -error code */
  int Problem_id,	/* in: approximation field ID  */
  int El_id,	/* in: unique identifier of the element */
  int Comp_sm,	/* in: indicator for the scope of computations: */
  /*   PDC_NO_COMP  - do not compute anything */
  /*   PDC_COMP_SM - compute entries to stiff matrix only */
  /*   PDC_COMP_RHS - compute entries to rhs vector only */
  /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int Pdeg_in,	/* in: enforced degree of polynomial (if > 0 ) */
  int *Nrdofs_loc,	/* in(optional): size of Stiff_mat and Rhs_vect */
  /* out(optional): actual number of dofs per integration entity */
  double *Stiff_mat,	/* out(optional): stiffness matrix stored columnwise */
  double *Rhs_vect,	/* out(optional): rhs vector */
  char *Rewr_dofs	/* out(optional): flag to rewrite or sum up entries */
  /*   'T' - true, rewrite entries when assembling */
  /*   'F' - false, sum up entries when assembling */
  );

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

  if(Problem_id==PDC_NS_SUPG_ID){


	char field_module_name[100];
	int solver_id, field_id, mesh_id, nreq;
	int nel, nfa, nno, int_ent, dof_ent, nr_elem, nr_face, nr_node, i, iaux, iedg;

	/* associate the suitable field with the problem and the solver */
	i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
	i=3; field_id = pdr_ctrl_i_params(Problem_id,i);
	i=6; solver_id = pdr_ctrl_i_params(Problem_id,i);
	//mesh_id = apr_get_mesh_id(field_id);
	nreq=apr_get_nreq(field_id);

	// check the name of the field module
	apr_module_introduce(field_module_name);

	if( strncmp(field_module_name, "STANDARD_LINEAR", 15) != 0 &&
	strncmp(field_module_name,"STANDARD_QUADRATIC",18) != 0 ){
	  printf("utr_get_list_int_ent works only for STANDARD (linear and quadratic) approximations!\n");
	  printf("Exiting.\n");
	  exit(-1);
	}

	/* prepare arrays */
	nr_elem = mmr_get_nr_elem(mesh_id);
	nr_face = mmr_get_nr_face(mesh_id);

	*List_int_ent_type = (int *) malloc( (nr_elem+nr_face+1)*sizeof(int) );
	*List_int_ent_id = (int *) malloc( (nr_elem+nr_face+1)*sizeof(int) );

	int_ent=0;

	/* loop over elements - integration entities*/
	nel=0;
	int is_fluid,face_neig[2];
	while((nel=mmr_get_next_elem_all(mesh_id, nel))!=0){

	  is_fluid = utr_mat_get_material(mmr_el_groupID(mesh_id,nel))->is_fluid;

	  if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE && is_fluid==FLUID) {

	int_ent++;
	(*List_int_ent_type)[int_ent-1]=PDC_ELEMENT;
	(*List_int_ent_id)[int_ent-1]=nel;

	  } // end if element included for integration

	} // end for all elements

	/*loop over faces*/
	nfa=0;
	while((nfa=mmr_get_next_face_all(mesh_id, nfa))!=0) {

	  if(mmr_fa_status(mesh_id,nfa)==MMC_ACTIVE) {
	if(mmr_fa_bc(mesh_id, nfa)>0) {






	  /* check whether the first neighbor is fluid (the second must be BC !!!!!) */
		mmr_fa_eq_neig(mesh_id, nfa, face_neig, NULL, NULL);
		is_fluid = utr_mat_get_material(mmr_el_groupID(mesh_id,face_neig[0]))->is_fluid;
		//is_fluid = utr_is_fluid(face_neig[0]);

		//if(mmr_is_fluid(face_neig[0]) == 1) {is_fluid = 1;}
		//else{ is_fluid = 0;}

		if(is_fluid==FLUID){

		int_ent++;
		(*List_int_ent_type)[int_ent-1]=PDC_FACE;
		(*List_int_ent_id)[int_ent-1]=nfa;

	  }
	}
	  }
	}

	*Nr_int_ent = int_ent;

	// default procedure for returning lists of dof entities
	// for standard and quadratic approximations
	utr_get_list_dof_ent(Problem_id, Nr_int_ent,
			 List_int_ent_type, List_int_ent_id,
			 Nr_dof_ent, List_dof_ent_type, List_dof_ent_id,
			 List_dof_ent_nrdofs, Nrdofs_glob, Max_dofs_per_dof_ent);


  }
  else if(Problem_id==PDC_HEAT_ID){

	// default procedure for returning lists of integration and dof entities
	// for standard and discontinuous Galerkin approximations
	utr_get_list_ent(Problem_id, Nr_int_ent, List_int_ent_type, List_int_ent_id,
			 Nr_dof_ent, List_dof_ent_type, List_dof_ent_id, List_dof_ent_nrdofs,
			 Nrdofs_glob, Max_dofs_per_dof_ent);

  }


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
  printf("pdr_get_list_ent_coarse NOT IMPLEMENTED in pdd_ns_supg_heat!");
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
				   L_dof_ent_id, L_dof_ent_nrdofs,
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
  pdr_dof_ent_sons - to return a list of dof entity sons
---------------------------------------------------------*/
int pdr_dof_ent_sons(		/* returns: success >=0 or <0 - error code */
  int Problem_id,	/* in: problem ID  */
  int Ent_type,	/* in: type of mesh entity */
  int Ent_id,	/* in: mesh entity ID */
  int *Ent_sons	/* out: list of dof entity sons */
				/*         Ent_sons[0] - number of sons */
	)
{

  int field_id, mesh_id;

  field_id = pdr_ctrl_i_params(Problem_id,3);
  mesh_id = apr_get_mesh_id(field_id);
  if (Ent_type == PDC_ELEMENT) {
/*kbw
	printf("pdr_dof_ent_sons: element %d, mesh_id %d\n", Ent_id, mesh_id);
/*kew*/

	mmr_el_fam(mesh_id, Ent_id, Ent_sons, NULL);
/*kbw
	  printf("sons:");
	  for(i=0;i<Ent_sons[0];i++){
		printf("%d  ", Ent_sons[i+1]);
	  }
	  printf("\n");
/*kew*/
  } else if (Ent_type == PDC_FACE) {
	mmr_fa_fam(mesh_id, Ent_id, Ent_sons, NULL);
  } else if (Ent_type == PDC_EDGE) {
	mmr_edge_sons(mesh_id, Ent_id, Ent_sons, NULL);
  } else {
	Ent_sons[0] = 0;
  }

  return (0);
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

  printf("pdr_proj_sol_lev NOT IMPLEMENTED in pdd_ns_supg_heat!");
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

/*kbw
 printf("Problem_id %d, Solver_id %d, Field_id %d, Level_id %d, Nrdof %d\n",
	Problem_id, Solver_id, field_id, Level_id, Nrdof);
/*kew*/

#ifdef PARALLEL
  if(Problem_id==PDC_NS_SUPG_ID){

	vec_norm = appr_sol_vec_norm(pdv_ns_supg_exchange_table_index,
				 Level_id, Nrdof, Vector);

  } else {

	vec_norm = appr_sol_vec_norm(pdv_heat_exchange_table_index,
				 Level_id, Nrdof, Vector);

  }
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
  if(Problem_id==PDC_NS_SUPG_ID){

	sc_prod = appr_sol_sc_prod(pdv_ns_supg_exchange_table_index,
				   Level_id, Nrdof, Vector1, Vector2);

  } else {

	sc_prod = appr_sol_sc_prod(pdv_heat_exchange_table_index,
				   Level_id, Nrdof, Vector1, Vector2);

  }
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
  if(Problem_id==PDC_NS_SUPG_ID){

	appr_create_exchange_tables(pdv_ns_supg_exchange_table_index,
				Level_id, Nr_dof_ent, L_dof_ent_type,
				L_dof_ent_id, L_bl_nrdof, L_bl_posg,
				L_elem_to_bl, L_face_to_bl,
				L_edge_to_bl, L_vert_to_bl);

  } else {

	appr_create_exchange_tables(pdv_heat_exchange_table_index,
				Level_id, Nr_dof_ent, L_dof_ent_type,
				L_dof_ent_id, L_bl_nrdof, L_bl_posg,
				L_elem_to_bl, L_face_to_bl,
				L_edge_to_bl, L_vert_to_bl);

  }
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
  if(Problem_id==PDC_NS_SUPG_ID){
	appr_exchange_dofs(pdv_ns_supg_exchange_table_index, Level_id, Vec_dofs);
  } else {
	appr_exchange_dofs(pdv_heat_exchange_table_index, Level_id, Vec_dofs);
  }
#endif

  return (1);
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

  // ns_supg_heat is a complex module formed using two basic modules - ns_supg
  // and heat. It is assumed that both modules use apr_num_int_el to perform
  // numerical integration. apr_num_int_el calls pdr_select_el_coeff to know
  // which terms are non-zero in the weak formulation. Here we direct only
  // the request to proper problem modules


  /* allocate storage for coefficients and select the needed ones */
  if(Problem_id==PDC_NS_SUPG_ID){
	pdr_ns_supg_select_el_coeff_vect(Problem_id, Coeff_vect_ind);
  }
  else{
	pdr_heat_select_el_coeff_vect(Problem_id, Coeff_vect_ind);
  }

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

  // ns_supg_heat is a complex module formed using two basic modules - ns_supg
  // and heat. It is assumed that both modules use apr_num_int_el to perform
  // numerical integration. apr_num_int_el calls pdr_select_el_coeff to know
  // which terms are non-zero in the weak formulation. Here we direct only
  // the request to proper problem modules

  double *select_coeff=NULL;

  /* allocate storage for coefficients and select the needed ones */
  if(Problem_id==PDC_NS_SUPG_ID){
	select_coeff=pdr_ns_supg_select_el_coeff(Problem_id, Mval,
	 Axx,Axy,Axz,Ayx,Ayy,Ayz,Azx,Azy,Azz,
	 Bx,By,Bz,Tx,Ty,Tz,Cval,Lval,Qx,Qy,Qz,Sval);
  }
  else{
	select_coeff=pdr_heat_select_el_coeff(Problem_id, Mval,
	 Axx,Axy,Axz,Ayx,Ayy,Ayz,Azx,Azy,Azz,
	 Bx,By,Bz,Tx,Ty,Tz,Cval,Lval,Qx,Qy,Qz,Sval);
  }

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

  // ns_supg_heat is a complex module formed using two basic modules - ns_supg
  // and heat. It is assumed that both modules use apr_num_int_el to perform
  // numerical integration. apr_num_int_el calls pdr_el_coeff to get coefficients
  // for different terms in weak formulation. Here we direct only
  // the request to proper problem modules and provide basic parameters

  if(Problem_id==PDC_NS_SUPG_ID){

// There are two possible ways to navigate through problem parameters
// 1. get problem structure and then get directly parameters
// (the method is FASTER but requires the knowledge of problem structure type)

	pdt_ns_supg_problem *problem =
	  (pdt_ns_supg_problem *)pdr_get_problem_structure(Problem_id);
	int field_id = problem->ctrl.field_id;
	int mesh_id = problem->ctrl.mesh_id;
	//int nreq = problem->ctrl.nreq;

// 2. get particular parameters using interface functions
// (the method is slower but can be used in modules that do not know the type
//  of problem structure for ns_supg problem)

	/* select the proper field */
	//field_id = pdr_ctrl_i_params(Problem_id, 3);
	//mesh_id = apr_get_mesh_id(field_id);
	//nreq =apr_get_nreq(field_id);
	//int nreq = pdr_ctrl_i_params(Problem_id, 5);

	/* get temperature at gauss point */
	double tk=0.0;
	pdr_ns_supg_give_me_temperature_at_point(Problem_id, Elem, X_loc, Base_phi,
							 Base_dphix, Base_dphiy, Base_dphiz,
							 &tk, NULL, NULL, NULL);

/*kbw
	printf("temperature at point %lf, %lf, %lf : %lf\n",
	   Xcoor[0], Xcoor[1], Xcoor[2], tk);
/*kew*/

	/*! ----------------------------------------------------------------------! */
	/*! -------------------- MATERIAL DATA AT GAUSS POINT ----------------- --! */
	/*! ----------------------------------------------------------------------! */
	double ref_temperature = problem->ctrl.ref_temperature;
	double dyn_visc, density, ref_density, VOF;

	// for Navier-Stokes with constant density and viscosity (ref_temperature <= 0)
	if(ref_temperature<=0){
	  dyn_visc = problem->ctrl.dynamic_viscosity;
	  density = problem->ctrl.density;
	  ref_density = problem->ctrl.density;
	}
	else{

	  // for Navier-Stokes with density or viscosity temperature dependent
	  utt_material_query_params qparams;
	  utt_material_query_result qresult;
	  int i;

	  //1.set query parameters (which material and temperature)
	  qparams.group_idx = Mat_num;	// query by material index ...
	  /* printf("\nAS: (pdr_ns_supg_el_coeff) mat = %d", Mat_num); */
	  qparams.name = "";	// ... not by material name
	  qparams.temperature = tk;	// current temperature
	  qparams.cell_id = Elem;
	  for( i=0; i<3; i++ ){
		qparams.xg[i] = Xcoor[i];
	  }
	  //2.get query results
	  pdr_ns_supg_material_query(&qparams, &qresult);
	  //3.set values to those obtained with query
	  dyn_visc = qresult.dynamic_viscosity;
	  //double tconductivity = qresult.thermal_conductivity;
	  //double specific_heat = qresult.specific_heat;
	  density = qresult.density;
	  //VOF = qresult.VOF;
	  //double thermal_diffusivity = tconductivity / ( density * specific_heat );
	  //double texpansion = qresult.thermal_expansion_coefficient;
	  //4.get reference density (using reference temperature)
	  // method 1:
	  // double ref_temperature = pdr_ctrl_d_params(Problem_id, 20);
	  // method 2:
	  qparams.temperature = ref_temperature;
	  //ref_density = density at reference temperature
	  //AS: the same as above qparams.xg's
	  pdr_ns_supg_material_query(&qparams, &qresult);
	  ref_density = qresult.density;

	}
#ifdef TURBULENTFLOW
	double tv;
	TurbulentViscosity(Elem,X_loc,&tv);
	dyn_visc+=tv;
#endif //TURBULENTFLOW

	double delta_t = problem->time.cur_dtime;
	double u_ref = problem->ctrl.ref_velocity;

/*kbw
  if(El_id==5){
	printf("\ndtime %lf, ref_density %lf, dynamic_viscosity %lf\n",
	   delta_t, ref_density, dyn_visc);
	printf("Uk: %lf, %lf, %lf\n", Uk_val[0], Uk_val[1], Uk_val[2]);
	printf("Un: %lf, %lf, %lf\n", Un_val[0], Un_val[1], Un_val[2]);
  }
/*kew */

	if(delta_t <= 0.0 || dyn_visc <= 0.0 || ref_density <= 0.0 || density <= 0.0){
	  printf("Error in setting parameters for ns_supg_el_coeff !!!\n");
	  printf("\ndtime %lf, density %lf, ref_density %lf, dynamic_viscosity %lf\n",
		 delta_t, density, ref_density, dyn_visc);
	  printf("Exiting!!!\n");
	}

	/*! ----------------------------------------------------------------------! */
	/*! ------------ CALCULATE BOUSINESSQUE FORCE AT GAUSS POINT -------------! */
	/*! ----------------------------------------------------------------------! */
	/*
	  /* !!!!!!! body force should be in [kg/m^2/sec^2] or equivalent */

	//OPTION 1: get density from texpansion, tk and ref_temperature
	/* f_x = Reference_density*(1.0-texp*(tk-ref_temperature))*g_pdv_cfg.gravity_field[0]; */
	/* f_y = Reference_density*(1.0-texp*(tk-ref_temperature))*g_pdv_cfg.gravity_field[1]; */
	/* f_z = Reference_density*(1.0-texp*(tk-ref_temperature))*g_pdv_cfg.gravity_field[2]; */

	  //OPTION 2: get density from (piecewise linear interpolated) values from file
	  // gravity_field components returned as double control parameters
	  // gravity_filed[i] = pdr_ctrl_d_params(Problem_id, 10+i);
	  //double f_x = density * pdr_ctrl_d_params(Problem_id, 10);
	  //double f_y = density * pdr_ctrl_d_params(Problem_id, 11);
	  //double f_z = density * pdr_ctrl_d_params(Problem_id, 12);

	double f_x, f_y, f_z;
	char problem_name[100];
	pdr_problem_name(Problem_id, problem_name);
	if(strcmp(problem_name, "HDC") == 0){

	  f_x = -(tk-1.0) * problem->ctrl.gravity_field[0];
	  f_y = -(tk-1.0) * problem->ctrl.gravity_field[1];
	  f_z = -(tk-1.0) * problem->ctrl.gravity_field[2];
	  //printf("HDC");
	}
	else{

	  f_x = density * problem->ctrl.gravity_field[0];
	  f_y = density * problem->ctrl.gravity_field[1];
	  f_z = density * problem->ctrl.gravity_field[2];

	}

  // HDC - gravity as Rayleigh number - set constant density = 1 in materials.dat
  //double f_x = -tk * problem->ctrl.gravity_field[0];
  //double f_y = -tk * problem->ctrl.gravity_field[1];
  //double f_z = -tk * problem->ctrl.gravity_field[2];

	 /*! ----------------------------------------------------------------------! */
	 /*! ---- CALCULATE FRICTIONAL FORCE IN THE MUSHY ZONE AT GAUSS POINT -----! */
	 /*! ----------------------------------------------------------------------! */

	 /* // Abderrazak_in_OptLasTechn_2009 */
	 /* // Voller_Swaminathan_NumHeatTransB_1991 */
	 /* double c_mz = 1.6E03;	// constant depend on the morphology of the porous media */
	 /* //double c_mz = 1.0E10;	// constant depend on the morphology of the porous media [1 / m^2] */
	 /* double q_mz = 1E-03;	// constant introduced to avoid division by zero */
	 /* //double S_mz = - c_mz * (1 - VOF) * (1 - VOF) / (VOF * VOF *VOF + q_mz); */
	 /* double S_mz; */
	 /* double fM_x, fM_y, fM_z; */

	 /* if( VOF > 0.0 ) { */
	 /* 	  S_mz = - c_mz * (1 - VOF) * (1 - VOF) / (VOF * VOF *VOF + q_mz); */
	 /* 	  fM_x = S_mz * dyn_visc * Uk_val[0]; */
	 /* 	  fM_y = S_mz * dyn_visc * Uk_val[1]; */
	 /* 	  fM_z = S_mz * dyn_visc * Uk_val[2]; */
	 /* 	} else { */
	 /* 	  fM_x = 0.0; */
	 /* 	  fM_y = 0.0; */
	 /* 	  fM_z = 0.0; */
	 /* } */

/*
	if ( VOF > 0.0 && VOF < 1.0)
	{
	  printf("\n\tAS: Eld_id = %d global coor %lf %lf %lf",El_id, xcoor[0], xcoor[1], xcoor[2]);
	  printf("\n\t\tAS: VOF = %lf, T = %lf [K]", VOF, tk);
	  printf("\n\t\tAS: fB_x = %lf, fB_y = %lf, fB_z = %lf", fB_x, fB_y, fB_z);
	  printf("\n\t\tAS: fM_x = %lf, fM_y = %lf, fM_z = %lf", fM_x, fM_y, fM_z);
	}
*/

/*! ----------------------------------------------------------------------! */
/*! ---------------------- CALCULATE TOTAL SOURCE TERMS ------------------! */
/*! ----------------------------------------------------------------------! */

	/* f_x += fM_x; */
	/* f_y += fM_y; */
	/* f_z += fM_z; */

		//double implicit = pdr_time_d_params(Problem_id, 2);
	  double implicit = problem->time.alpha;

	  /* get coefficients for ns_supg weak formulation */
	  pdr_ns_supg_el_coeff(Problem_id, Elem, Mat_num, Hsize, Pdeg, X_loc,
			   Base_phi, Base_dphix, Base_dphiy, Base_dphiz,
			   Xcoor, Uk_val, Uk_x, Uk_y, Uk_z, Un_val, Un_x, Un_y, Un_z,
			   Mval, Axx, Axy, Axz, Ayx, Ayy, Ayz, Azx, Azy, Azz,
			   Bx, By, Bz, Tx, Ty, Tz, Cval, Lval, Qx, Qy, Qz, Sval,
			   tk, dyn_visc, density, ref_density, u_ref,
			   delta_t, implicit, f_x, f_y, f_z);
  } // end if ns_supg problem
  else{ // heat problem
	/* get velocity at gauss point from ns_supg module */
	double vel[3];
	pdr_ns_supg_get_velocity_at_point(PDC_NS_SUPG_ID, Elem, X_loc, Base_phi,
					  NULL, NULL, NULL,
					  //Base_dphix, Base_dphiy, Base_dphiz,
					  vel, NULL, NULL, NULL);
/*kbw
	printf("velocity at point %lf, %lf, %lf : %lf, %lf, %lf \n",
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
	  printf("wrong parameter HEAT_NREQ in pdr_ns_supg_heat_el_coeff 1\n");
	  printf("%d != %d. Exiting !!!",nreq, apr_get_nreq(field_id));
	  exit(-1);
	}
	if(nreq != pdr_ctrl_i_params(Problem_id,5)){
	  printf("wrong parameter HEAT_NREQ in pdr_ns_supg_heat_el_coeff 2\n");
	  printf("%d != %d. Exiting !!!",nreq, pdr_ctrl_i_params(Problem_id,5));
	  exit(-1);
	}
	if(nreq != problem->ctrl.nreq){
	  printf("wrong parameter HEAT_NREQ in pdr_ns_supg_heat_el_coeff 3\n");
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

	  /* printf("\nAS: (pdr_ns_supg_heat_el_coeff) mat = %d", Mat_num); */
	  qparams.name = "";	//... not by material name
	  double tk = Uk_val[0]; // temperature is the only unknown
	  qparams.temperature = tk;	//temperature from last iteration
	  qparams.cell_id = Elem;
	  for( i=0; i<3; i++ ){
	qparams.xg[i] = Xcoor[i];
	  }
	  //2.get query results
	  pdr_heat_material_query(&qparams, &qresult);
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
  }

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
int pdr_comp_stiff_mat_uncon(
					   /** returns: >=0 - success code, <0 - error code */
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

	pdr_ns_supg_heat_comp_el_stiff_mat(Problem_id, Int_ent_id, Comp_sm, pdeg,
			  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  } else if (Int_ent_type == PDC_FACE) {

	pdr_ns_supg_heat_comp_fa_stiff_mat(Problem_id, Int_ent_id, Comp_sm, pdeg,
				   Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  } else {
	printf("ERROR: Wrong integration entity type in pdr_comp_stiff_mat!\n");
	exit(-1);
  }

  return(0);
}


/*------------------------------------------------------------
  pdr_ns_supg_heat_comp_el_stiff_mat - to construct a stiffness matrix and
						  a load vector for an element
------------------------------------------------------------*/
int pdr_ns_supg_heat_comp_el_stiff_mat(
			/*returns: >=0 - success code, <0 -error code */
  int Problem_id,	/* in: approximation field ID  */
  int El_id,	/* in: unique identifier of the element */
  int Comp_sm,	/* in: indicator for the scope of computations: */
  /*   PDC_NO_COMP  - do not compute anything */
  /*   PDC_COMP_SM - compute entries to stiff matrix only */
  /*   PDC_COMP_RHS - compute entries to rhs vector only */
  /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int Pdeg_in,	/* in: enforced degree of polynomial (if > 0 ) */
  int *Nrdofs_loc,	/* in(optional): size of Stiff_mat and Rhs_vect */
  /* out(optional): actual number of dofs per integration entity */
  double *Stiff_mat,	/* out(optional): stiffness matrix stored columnwise */
  double *Rhs_vect,	/* out(optional): rhs vector */
  char *Rewr_dofs	/* out(optional): flag to rewrite or sum up entries */
  /*   'T' - true, rewrite entries when assembling */
  /*   'F' - false, sum up entries when assembling */
  )
{

  // ns_supg_heat is a complex module formed using two basic modules - ns_supg
  // and heat. There is only one pdr_comp_el_stiff_mat procedure in the code
  // so we can either compute one large stiffness matrix for one solver - in
  // the current procedure (the approach followed in hmt module) or by
  // combination of matrices produced in ns_supg and heat modules
  //  - or use two problems with two solvers. Here, we follow the last approach,
  // so the current pdr_comp_el_stiff_mat calls only the suitable procedures
  // from the ns_supg and heat modules and passes the matrices to the linear
  // solver associated with the proper problem

  // call the procedure that computes UNCONSTRAINED element stiffness matrix
  if (Problem_id == PDC_NS_SUPG_ID) {

	pdr_ns_supg_comp_el_stiff_mat(Problem_id, El_id, Comp_sm, Pdeg_in,
			  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  }
  else{

	pdr_heat_comp_el_stiff_mat(Problem_id, El_id, Comp_sm, Pdeg_in,
			  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  }

  return (1);
}



/*------------------------------------------------------------
  pdr_ns_supg_heat_comp_fa_stiff_mat - to construct a stiffness matrix and
						  a load vector for a face
------------------------------------------------------------*/
int pdr_ns_supg_heat_comp_fa_stiff_mat(
			/*returns: >=0 -success code, <0 -error code */
  int Problem_id,	/* in: approximation field ID  */
  int Fa_id,	/* in: unique identifier of the face */
  int Comp_sm,	/* in: indicator for the scope of computations: */
  /*   PDC_NO_COMP  - do not compute anything */
  /*   PDC_COMP_SM - compute entries to stiff matrix only */
  /*   PDC_COMP_RHS - compute entries to rhs vector only */
  /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int Pdeg_in,	/* in: enforced degree of polynomial (if > 0 ) */
  int *Nrdofs_loc,	/* in(optional): size of Stiff_mat and Rhs_vect */
  /* out(optional): actual number of dofs per integration entity */
  double *Stiff_mat,	/* out(optional): stiffness matrix stored columnwise */
  double *Rhs_vect,	/* out(optional): rhs vector */
  char *Rewr_dofs	/* out(optional): flag to rewrite or sum up entries */
  /*   'T' - true, rewrite entries when assembling */
  /*   'F' - false, sum up entries when assembling */
  )
{

  // ns_supg_heat is a complex module formed using two basic modules - ns_supg
  // and heat. There is only one pdr_comp_fa_stiff_mat procedure in the code
  // so we can either compute one large stiffness matrix for one solver - in
  // the current procedure (the approach followed in hmt module) or by
  // combination of matrices produced in ns_supg and heat modules
  //  - or use two problems with two solvers. Here, we follow the last approach,
  // so the current pdr_comp_fa_stiff_mat calls only the suitable procedures
  // from the ns_supg and heat modules and passes the matrices to the linear
  // solver associated with the proper problem

  // call the procedure that computes UNCONSTRAINED element stiffness matrix
  // that contains terms resulting from boundary conditions on face Fa_id
  if(Problem_id == PDC_NS_SUPG_ID){

	pdr_ns_supg_comp_fa_stiff_mat(Problem_id, Fa_id, Comp_sm, Pdeg_in,
				  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  }
  else{

	pdr_heat_comp_fa_stiff_mat(Problem_id, Fa_id, Comp_sm, Pdeg_in,
				   Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  }


  return (1);
}


/*------------------------------------------------------------
  pdr_ns_supg_give_me_temperature_at_point - to provide the temperature and its
	gradient at a particular point given its local coordinates within an element
NS_SUPG MODULE ASKS FOR IMPLEMENTATION - it has to be provided by procedures
defined in ls_intf directory of the problem module that uses ns_supg as submodule
------------------------------------------------------------*/
int pdr_ns_supg_give_me_temperature_at_point(
  int Problem_id,
  int El_id, // element
  double *X_loc, // local coordinates of point
  double *Base_phi, // shape functions at point (if available - to speed up)
  double *Base_dphix, // derivatives of shape functions at point
  double *Base_dphiy, // derivatives of shape functions at point
  double *Base_dphiz, // derivatives of shape functions at point
  double *Temp, // temperature
  double *DTemp_dx, // x-derivative of temperature
  double *DTemp_dy, // y-derivative of temperature
  double *DTemp_dz // z-derivative of temperature
  )
{

  // we call heat problem to specify temperature
  pdr_heat_get_temperature_at_point(PDC_HEAT_ID, El_id, X_loc, Base_phi,
					Base_dphix, Base_dphiy, Base_dphiz,
					Temp, DTemp_dx, DTemp_dy, DTemp_dz);

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

  // we call ns_supg problem to specify velocity
  pdr_ns_supg_get_velocity_at_point(PDC_NS_SUPG_ID, El_id, X_loc, Base_phi,
					Base_dphix, Base_dphiy, Base_dphiz,
					Velocity, DVel_dx, DVel_dy, DVel_dz);

  return(0);
}


#ifdef __cplusplus
}
#endif
