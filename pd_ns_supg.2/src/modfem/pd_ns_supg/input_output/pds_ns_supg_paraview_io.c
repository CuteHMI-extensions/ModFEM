/************************************************************************
File pds_ns_supg_paraview_io.c - paraview output generation

Contains definition of routines:
  pdr_ns_supg_write_paraview

------------------------------
History:
	2011    - Przemyslaw Plaszewski (pplaszew@agh.edu.pl)
	2011    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h> /* USES */
/* interface for all approximation modules */
#include <modfem/aph_intf.h> /* USES */
/* utility procedures for all problem modules */
#include <modfem/uth_intf.h> // USES
/* problem dependent module interface */
#include <modfem/pdh_intf.h>		/* USES */
// Paralle communication header
#include <modfem/pch_intf.h> // USES
/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>	/* IMPLEMENTS */

#include <modfem/uth_io_results.h>

/* problem module's types and functions */
#include <modfem/pd_ns_supg/pdh_ns_supg.h> /* USES & IMPLEMENTS */
/* types and functions related to problem structures */
#include <modfem/pd_ns_supg/pdh_ns_supg_problem.h>
// bc and material header files are included in problem header files


/**************************************/
/* LOCAL TYPES                              */
/**************************************/
/* Rules:
/* - type name starts always witn pdt_ */

typedef struct{
  double dofs[PDC_MAXEQ];
}pdt_sol_info;

/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */

/*---------------------------------------------------------
pdr_ns_supg_write_paraview - write mesh in ParaView format
---------------------------------------------------------*/
int pdr_ns_supg_write_paraview(
  char* Work_dir,
  FILE *Interactive_input,
  FILE *Interactive_output
)
{
  int i;
  int num_nodes;
  int num_elems;
  double node_coor[3];
  int num_conn = 0;
  int mesh_id, field_id;
  char filename[300];
  int dofs_write[PDC_MAXEQ+1];

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* open the output file */

#ifdef PARALLEL
  sprintf(filename, "ns_supg_%d_proc%d",
		  pdv_ns_supg_problem.time.cur_step, pcr_my_proc_id());
#else
  sprintf(filename, "ns_supg_%d",
	  pdv_ns_supg_problem.time.cur_step);

#endif

  const ute_paraview_flags VTK_FORMAT = UTE_VTK_XML;

  utr_io_result_add_value(RESULT_OUT_FILENAME,filename );

  mesh_id = pdr_ctrl_i_params(PDC_NS_SUPG_ID, 2);

  utr_write_paraview_mesh(mesh_id, Work_dir,filename, VTK_FORMAT);

  utt_paraview_field_descriptor field_desc[2];
  field_desc[0].field_id = pdv_ns_supg_problem.ctrl.field_id;
  field_desc[0].dofs_write[0]=1;
  field_desc[0].dofs_write[1]=3; // pressure is the fourth component (3rd counting from 0)
  field_desc[0].entity_type = UTE_POINT_DATA;
  field_desc[0].quantity_type = UTE_SCALARS;
  field_desc[0].value_type = UTE_DOUBLE;
  field_desc[0].field_name = "pressure";

  field_desc[1].field_id = pdv_ns_supg_problem.ctrl.field_id;
  field_desc[1].dofs_write[0]=3;
  field_desc[1].dofs_write[1]=0;
  field_desc[1].dofs_write[2]=1;
  field_desc[1].dofs_write[3]=2; // velocity components
  field_desc[1].entity_type = UTE_POINT_DATA;
  field_desc[1].quantity_type = UTE_VECTORS;
  field_desc[1].value_type = UTE_DOUBLE;
  field_desc[1].field_name = "velocity";

  utr_write_paraview_fields(Work_dir,filename,pdv_ns_supg_problem.time.cur_time,
							2, field_desc, VTK_FORMAT);


   // write bc info
  utr_write_paraview_bc(mesh_id,Work_dir,filename,VTK_FORMAT);

  // other options directly using field data structures as in
  // utr_write_paraview_field


  return 0;
}

/* THE CODE BELOW IS LEFT AS AN EXAMPLE FOR EXPLICIT WRITING OF PARAVIEW FILES */
/* ITS PARTS CAN BE ADAPTED TO THE PROBLEM SOLVED AND MOVED ABOVE */

/*---------------------------------------------------------
pdr_paraview_write_field - write field in ParaView format
---------------------------------------------------------*/
/* int pdr_paraview_write_field(int Problem_id, int Field_id, char *Filename) */
/* { */
/*     FILE *fp; */
/*     int i; */
/*     int mesh_id; */
/*     pdt_sol_info *sol_infos; */
/*     int num_nodes, *check; */
/*     /\* material data (queries) *\/ */
/*     pdt_material_query_params qparams; */
/*     pdt_material_query_result qresult; */
/*     /\* material data *\/ */
/*     double viscosity, tconductivity, specific_heat, texpansion, density, enthalpy, eresistivity, VOF; */


/*   pdt_problem *problem = &pdv_problems[Problem_id]; */
/*   pdt_ctrls *ctrls = &pdv_problems[Problem_id].ctrl; */

/* /\* open the output file *\/ */
/*     fp = fopen(Filename, "w"); */
/*     if (fp == NULL) { */
/* 	fprintf(ctrls->interactive_output, "Cannot open file '%s' for Paraview field data\n", Filename); */
/* 	return (-1); */
/*     } */

/*     mesh_id = apr_get_mesh_id(Field_id); */
/*     //num_nodes = mmr_mesh_i_params(mesh_id, 1); */
/*     num_nodes = mmr_get_nr_node(mesh_id); */
/*     sol_infos = (pdt_sol_info *) malloc((num_nodes + 1) * sizeof(pdt_sol_info)); */
/*     check = (int *) malloc((num_nodes + 1) * sizeof(int)); */

/*     for (i = 1; i <= num_nodes; i++) { */
/* 	check[i] = apr_read_ent_dofs(Field_id, APC_VERTEX, i, 5, 1, sol_infos[i].dofs); */
/*     } */

/*     fprintf(fp, "POINT_DATA %d\n", num_nodes); */
/*     fprintf(fp, "SCALARS pressure double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/* 	if (check[i] != 2) */
/* 	    fprintf(fp, "%.12lg\n", sol_infos[i].dofs[3]); */
/*     } */

/*     fprintf(fp, "SCALARS temperature double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2) */
/*             fprintf(fp, "%.12lg\n", sol_infos[i].dofs[4]); */
/*     }  */

/*     //1.set query parameters (which material and temperature) */
/*     qparams.material_idx = 0;   //query by material index ... */
/*     qparams.name = "";          //... not by material name //TODO: get material by mesh coordinates */

/*     fprintf(fp, "SCALARS density double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2){ */
/*             qparams.temperature = sol_infos[i].dofs[4]; //temperature from last iteration */
/*             //2.get query results */
/*             pdr_material_query(&problem->materials, &qparams, &qresult); */
/*             //3.set values to those obtained with query */
/*             density = qresult.density; */
/*             fprintf(fp, "%.12lg\n", density); */
/* 	} */
/*     } */

/*     fprintf(fp, "SCALARS viscosity double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2){ */
/*             qparams.temperature = sol_infos[i].dofs[4]; //temperature from last iteration */
/*             //2.get query results */
/*             pdr_material_query(&problem->materials, &qparams, &qresult); */
/*             //3.set values to those obtained with query */
/*             viscosity = qresult.viscosity; */
/*             fprintf(fp, "%.12lg\n", viscosity); */
/* 	} */
/*     } */

/*     fprintf(fp, "SCALARS thermal_conductivity double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2){ */
/*             qparams.temperature = sol_infos[i].dofs[4]; //temperature from last iteration */
/*             //2.get query results */
/*             pdr_material_query(&problem->materials, &qparams, &qresult); */
/*             //3.set values to those obtained with query */
/*             tconductivity = qresult.thermal_conductivity; */
/*             fprintf(fp, "%.12lg\n", tconductivity); */
/* 	} */
/*     } */

/*     fprintf(fp, "SCALARS specific_heat double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2){ */
/*             qparams.temperature = sol_infos[i].dofs[4]; //temperature from last iteration */
/*             //2.get query results */
/*             pdr_material_query(&problem->materials, &qparams, &qresult); */
/*             //3.set values to those obtained with query */
/*             specific_heat = qresult.specific_heat; */
/*             fprintf(fp, "%.12lg\n", specific_heat); */
/* 	} */
/*     } */

/*     fprintf(fp, "SCALARS thermal_expansion_coefficient double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2){ */
/*             qparams.temperature = sol_infos[i].dofs[4]; //temperature from last iteration */
/*             //2.get query results */
/*             pdr_material_query(&problem->materials, &qparams, &qresult); */
/*             //3.set values to those obtained with query */
/*             texpansion = qresult.thermal_expansion_coefficient; */
/*             fprintf(fp, "%.12lg\n", texpansion); */
/* 	} */
/*     } */

/*     fprintf(fp, "SCALARS enthalpy double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2){ */
/*             qparams.temperature = sol_infos[i].dofs[4]; //temperature from last iteration */
/*             //2.get query results */
/*             pdr_material_query(&problem->materials, &qparams, &qresult); */
/*             //3.set values to those obtained with query */
/*             enthalpy = qresult.enthalpy; */
/*             fprintf(fp, "%.12lg\n", enthalpy); */
/* 	} */
/*     } */

/*     fprintf(fp, "SCALARS electrical_resistivity double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2){ */
/*             qparams.temperature = sol_infos[i].dofs[4]; //temperature from last iteration */
/*             //2.get query results */
/*             pdr_material_query(&problem->materials, &qparams, &qresult); */
/*             //3.set values to those obtained with query */
/*             eresistivity = qresult.electrical_resistivity; */
/*             fprintf(fp, "%.12lg\n", eresistivity); */
/* 	} */
/*     } */

/*     fprintf(fp, "SCALARS VOF double 1\n"); */
/*     fprintf(fp, "LOOKUP_TABLE default\n"); */

/*     for (i = 1; i <= num_nodes; i++) { */
/*         if (check[i] != 2){ */
/*             qparams.temperature = sol_infos[i].dofs[4]; //temperature from last iteration */
/*             //2.get query results */
/*             pdr_material_query(&problem->materials, &qparams, &qresult); */
/*             //3.set values to those obtained with query */
/*             VOF = qresult.VOF; */
/*             fprintf(fp, "%.12lg\n", VOF); */
/* 	} */
/*     } */

/*     fprintf(fp, "VECTORS velocity double\n"); */
/*     for (i = 1; i <= num_nodes; i++) { */
/* 	if (check[i] != 2) */
/* 	    fprintf(fp, "%.12lg %.12lg %.12lg\n", sol_infos[i].dofs[0], sol_infos[i].dofs[1], sol_infos[i].dofs[2]); */
/*     } */

/*     fclose(fp); */

/*     free(sol_infos); */
/*     free(check); */

/*     return 0; */
/* } */

