/************************************************************************
File pds_heat_paraview_io.c - paraview output generation

Contains definition of routines:
  pdr_heat_write_paraview

------------------------------
History:
	2011    - Przemyslaw Plaszewski (pplaszew@agh.edu.pl)
	2011    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
	2016    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h> /* USES */
/* interface for all approximation modules */
#include <modfem/aph_intf.h> /* USES */
/* problem dependent module interface */
#include <modfem/pdh_intf.h>		/* USES */
/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>	/* USES */
// Parallel communication
#include <modfem/pch_intf.h> // USES

#include <modfem/uth_io_results.h>

/* problem module's types and functions */
#include <modfem/pd_heat/pdh_heat.h> /* USES & IMPLEMENTS */
/* types and functions related to problem structures */
#include <modfem/pd_heat/pdh_heat_problem.h>
// bc and material header files are included in problem header files
#include <modfem/uth_mat.h>

/**************************************/
/* LOCAL TYPES                              */
/**************************************/
/* Rules:
/* - type name starts always witn pdt_ */

typedef struct {
	double dofs[PDC_MAXEQ];
} pdt_sol_info;

typedef struct {
	int elNodes[7];
	int numNodes;
	int elType;
} PVElInfo;

typedef struct {
	double dofs[4];
} PVSolInfo;

/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */

/// pdr_heat_write_paraview_field - to dump temperature dependent field in Paraview format
/**     \param Mesh_id - id of the mesh associated with given field Field_id
	\param Field_id - id of the field to dump
	\param Filename - c-string with name of the file to write on disk
	\param Desc - c-string c-array with name of the field values
	(at least Desc[0]!= NULL should be passed)
 */

int pdr_heat_write_temperature_dpendent_fields(
		int Field_id,
		const char * Filename,
		int * Dofs_write // dofs to write: dofs_write[0] - number of dofs
		//                dofs_write[i] - IDs of dofs to write
)
{
	int i, nrdofs, idofs;
	PVSolInfo * solInfos;
	int numNodes, max_node_id, numElems, nreq, ino, idof, el_id, nrnodes;
	int mesh_id, el_type;
	int el_nodes[10];
	//int el_nodes[APC_MAX_GEO_DOFS+1];        /* list of nodes of El */
	//int el_nodes_type[APC_MAX_GEO_DOFS+1];        /* list of nodes type of El */

	double el_dofs[APC_MAXELSD];  /* solution dofs in an element */
	utt_material_query_params qparams;
	utt_material_query_result qresult;

	mesh_id = apr_get_mesh_id(Field_id);
	max_node_id = mmr_get_max_node_id(mesh_id);
	numElems = mmr_get_nr_elem(mesh_id);
	nreq = apr_get_nreq(Field_id);
	if (nreq > PDC_MAXEQ) {
		printf("Number of equations %d > maximum allowed in uts_dump_paraview %d\n",
				nreq, PDC_MAXEQ);
		printf(" Recompile pdr_heat_write with greater PDC_MAXEQ. Exiting!\n");
		exit(0);
	}

	// dof_entities (nodes) are numbered from 1 to max_node_id
	solInfos = (PVSolInfo *)malloc((max_node_id + 1) * sizeof(PVSolInfo));
	for (ino = 0; ino <= max_node_id; ino++) {
		for (idof = 0; idof < nreq; idof++) {
			solInfos[ino].dofs[idof] = 0.0;
		}
	}

	/* loop over elements */
	el_id = 0;
	while ((el_id = mmr_get_next_act_elem(mesh_id, el_id)) != 0) {
		int dof_counter = 0;
		nrnodes = mmr_el_node_coor(mesh_id, el_id, el_nodes, NULL);
		//nrnodes = apr_get_el_geo_dofs(field_id,el_id,el_nodes,el_nodes_type,NULL);


		apr_get_el_dofs(Field_id, el_id, 1, el_dofs);
		for (ino = 1; ino <= nrnodes; ino++) {
			nrdofs = apr_get_ent_nrdofs(Field_id, APC_VERTEX, ino);
			for (idof = 0; idof < nreq; idof++) {
				int node_id = el_nodes[ino];
				solInfos[node_id].dofs[idof] = el_dofs[dof_counter];
				//Fk-for testing material
				//solInfos[node_id].dofs[idof]=mmr_el_groupID(mesh_id, el_id);
				/*kbw
					printf("node %d, value %.12lg\n",node_id, solInfos[node_id].dofs[0]);
				/*kew*/
				dof_counter++;
			}
		}
	}

	/****************************************************************************/
// everything below commented out since it does not take care of material number !!!

	/* pdt_heat_problem *problem_heat = (pdt_heat_problem *)pdr_get_problem_structure(pdv_heat_current_problem_id); */
	/* qparams.material_idx = 0; */
	/* qparams.name = ""; */

	/* fprintf(File, "SCALARS density double 1\n"); */
	/* fprintf(File, "LOOKUP_TABLE default\n"); */
	/* for(i=1;i<=max_node_id;i++){     */
	/*   for(idofs = 1; idofs <= Dofs_write[0]; idofs++){ */
	// CHECK WHETHER MATERIAL DATABASE USED AT ALL
	/*     qparams.temperature = solInfos[i].dofs[Dofs_write[idofs]]; */
	/*     pdr_heat_material_query(&problem_heat->materials, &qparams, &qresult); */
	/*     fprintf(File,"%.12lg\n",qresult.density); */
	/*   } */
	/*   fprintf(File,"\n"); */
	/* } */

	/* fprintf(File, "SCALARS thermal_conductivity double 1\n"); */
	/* fprintf(File, "LOOKUP_TABLE default\n"); */
	/* for(i=1;i<=max_node_id;i++){     */
	/*   for(idofs = 1; idofs <= Dofs_write[0]; idofs++){ */
	// CHECK WHETHER MATERIAL DATABASE USED AT ALL
	/*     qparams.temperature = solInfos[i].dofs[Dofs_write[idofs]]; */
	/*     pdr_heat_material_query(&problem_heat->materials, &qparams, &qresult); */
	/*     fprintf(File,"%.12lg\n",qresult.thermal_conductivity); */
	/*   } */
	/*   fprintf(File,"\n"); */
	/* } */

	/* fprintf(File, "SCALARS specific_heat double 1\n"); */
	/* fprintf(File, "LOOKUP_TABLE default\n"); */
	/* for(i=1;i<=max_node_id;i++){     */
	/*   for(idofs = 1; idofs <= Dofs_write[0]; idofs++){ */
	// CHECK WHETHER MATERIAL DATABASE USED AT ALL
	/*     qparams.temperature = solInfos[i].dofs[Dofs_write[idofs]]; */
	/*     pdr_heat_material_query(&problem_heat->materials, &qparams, &qresult); */
	/*     fprintf(File,"%.12lg\n",qresult.specific_heat); */
	/*   } */
	/*   fprintf(File,"\n"); */
	/* } */

	/* fprintf(File, "SCALARS enthalpy double 1\n"); */
	/* fprintf(File, "LOOKUP_TABLE default\n"); */
	/* for(i=1;i<=max_node_id;i++){     */
	/*   for(idofs = 1; idofs <= Dofs_write[0]; idofs++){ */
	// CHECK WHETHER MATERIAL DATABASE USED AT ALL
	/*     qparams.temperature = solInfos[i].dofs[Dofs_write[idofs]]; */
	/*     pdr_heat_material_query(&problem_heat->materials, &qparams, &qresult); */
	/*     fprintf(File,"%.12lg\n",qresult.enthalpy); */
	/*   } */
	/*   fprintf(File,"\n"); */
	/* } */

	/* fprintf(File, "SCALARS thermal_expansion_coefficient double 1\n"); */
	/* fprintf(File, "LOOKUP_TABLE default\n"); */
	/* for(i=1;i<=max_node_id;i++){     */
	/*   for(idofs = 1; idofs <= Dofs_write[0]; idofs++){ */
	// CHECK WHETHER MATERIAL DATABASE USED AT ALL
	/*     qparams.temperature = solInfos[i].dofs[Dofs_write[idofs]]; */
	/*     pdr_heat_material_query(&problem_heat->materials, &qparams, &qresult); */
	/*     fprintf(File,"%.12lg\n",qresult.thermal_expansion_coefficient); */
	/*   } */
	/*   fprintf(File,"\n"); */
	/* } */

	/* fprintf(File, "SCALARS electrical_resistivity double 1\n"); */
	/* fprintf(File, "LOOKUP_TABLE default\n"); */
	/* for(i=1;i<=max_node_id;i++){     */
	/*   for(idofs = 1; idofs <= Dofs_write[0]; idofs++){ */
	// CHECK WHETHER MATERIAL DATABASE USED AT ALL
	/*     qparams.temperature = solInfos[i].dofs[Dofs_write[idofs]]; */
	/*     pdr_heat_material_query(&problem_heat->materials, &qparams, &qresult); */
	/*     fprintf(File,"%.12lg\n",qresult.electrical_resistivity); */
	/*   } */
	/*   fprintf(File,"\n"); */
	/* } */


	/* fprintf(File, "SCALARS VOF double 1\n"); */
	/* fprintf(File, "LOOKUP_TABLE default\n"); */
	/* for(i=1;i<=max_node_id;i++){     */
	/*   for(idofs = 1; idofs <= Dofs_write[0]; idofs++){ */
	// CHECK WHETHER MATERIAL DATABASE USED AT ALL
	/*     qparams.temperature = solInfos[i].dofs[Dofs_write[idofs]]; */
	/*     pdr_heat_material_query(&problem_heat->materials, &qparams, &qresult); */
	/*     fprintf(File,"%.12lg\n",qresult.VOF); */
	/*   } */
	/*   fprintf(File,"\n"); */
	/* } */

	/****************************************************************************/


	free(solInfos);
	return 0;
}

/*---------------------------------------------------------
pdr_heat_write_paraview - write mesh in ParaView format
---------------------------------------------------------*/
int pdr_heat_write_paraview(
		char * Work_dir,
		FILE * Interactive_input,
		FILE * Interactive_output
)
{
	int i;
	int num_nodes;
	int num_elems;
	double node_coor[3];
	int num_conn = 0;
	int mesh_id, field_id;
	char filename[300];
	int dofs_write[PDC_MAXEQ + 1];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* open the output file */
#ifdef PARALLEL
	sprintf(filename, "heat_%d_proc%4d",  //pdv_heat_problem.ctrl.name,
			pdv_heat_problem.time.cur_step, pcr_my_proc_id());
#else
	sprintf(filename, "heat_%06d", //pdv_heat_problem.ctrl.name,
			pdv_heat_problem.time.cur_step);

#endif
	int problem_id = pdv_heat_current_problem_id;
	mesh_id = pdr_ctrl_i_params(problem_id, 2);

	const ute_paraview_flags VTK_FORMAT = UTE_VTK_XML;

	utr_io_result_add_value(RESULT_OUT_FILENAME, filename );

	utr_write_paraview_mesh(mesh_id, Work_dir, filename, VTK_FORMAT);

	// ---- Write data to paraview file ----
	utt_paraview_field_descriptor field_desc[99];
	int field_desc_nb = 0;

	// write pdd_heat temperature field data
	field_desc[0].field_id = pdv_heat_problem.ctrl.field_id;
	field_desc[0].dofs_write[0] = 1;
	field_desc[0].dofs_write[1] = 0; // temperature is the only component
	field_desc[0].entity_type = UTE_POINT_DATA;
	field_desc[0].quantity_type = UTE_SCALARS;
	field_desc[0].value_type = UTE_DOUBLE;
	field_desc[0].field_name = "temperature";
	field_desc_nb++;

#ifdef PHASE_TRANSFORMATION

	// write pdd_heat phase tranformation field data
	field_desc[1].field_id = pdv_heat_problem.ctrl.phase_transformation_field_id;
	field_desc[1].dofs_write[0] = 1;
	field_desc[1].dofs_write[1] = 0;
	field_desc[1].entity_type = UTE_POINT_DATA;
	field_desc[1].quantity_type = UTE_SCALARS;
	field_desc[1].value_type = UTE_DOUBLE;
	field_desc[1].field_name = "phase_transf";
	field_desc_nb++;

	// write pdd_heat austenite field data
	field_desc[2].field_id = pdv_heat_problem.ctrl.phases_field_id;
	field_desc[2].dofs_write[0] = 1;
	field_desc[2].dofs_write[1] = P_HEAT_AUSTENITE;
	field_desc[2].entity_type = UTE_POINT_DATA;
	field_desc[2].quantity_type = UTE_SCALARS;
	field_desc[2].value_type = UTE_DOUBLE;
	field_desc[2].field_name = "austenite";
	field_desc_nb++;

	// write pdd_heat martensite field data
	field_desc[3].field_id = pdv_heat_problem.ctrl.phases_field_id;
	field_desc[3].dofs_write[0] = 1;
	field_desc[3].dofs_write[1] = P_HEAT_MARTENSITE;
	field_desc[3].entity_type = UTE_POINT_DATA;
	field_desc[3].quantity_type = UTE_SCALARS;
	field_desc[3].value_type = UTE_DOUBLE;
	field_desc[3].field_name = "martensite";
	field_desc_nb++;

	// write pdd_heat bainite field data
	field_desc[4].field_id = pdv_heat_problem.ctrl.phases_field_id;
	field_desc[4].dofs_write[0] = 1;
	field_desc[4].dofs_write[1] = P_HEAT_BAINITE;
	field_desc[4].entity_type = UTE_POINT_DATA;
	field_desc[4].quantity_type = UTE_SCALARS;
	field_desc[4].value_type = UTE_DOUBLE;
	field_desc[4].field_name = "bainite";
	field_desc_nb++;

	// write pdd_heat ferrite field data
	field_desc[5].field_id = pdv_heat_problem.ctrl.phases_field_id;
	field_desc[5].dofs_write[0] = 1;
	field_desc[5].dofs_write[1] = P_HEAT_FERRITE;
	field_desc[5].entity_type = UTE_POINT_DATA;
	field_desc[5].quantity_type = UTE_SCALARS;
	field_desc[5].value_type = UTE_DOUBLE;
	field_desc[5].field_name = "bainite_Qv";  // for testing only here
	field_desc_nb++;

	// write pdd_heat pearlite field data
	field_desc[6].field_id = pdv_heat_problem.ctrl.phases_field_id;
	field_desc[6].dofs_write[0] = 1;
	field_desc[6].dofs_write[1] = P_HEAT_PEARLITE;
	field_desc[6].entity_type = UTE_POINT_DATA;
	field_desc[6].quantity_type = UTE_SCALARS;
	field_desc[6].value_type = UTE_DOUBLE;
	field_desc[6].field_name = "martensite_Qv";  // for testing only here
	field_desc_nb++;

#endif

	utr_write_paraview_fields(Work_dir, filename, pdv_heat_problem.time.cur_time,
			field_desc_nb, field_desc, VTK_FORMAT);


#ifdef PHASE_TRANSFORMATION_OLD

	// OLD - TO BE CHANGED

	utr_write_paraview_field(pdv_heat_problem.ctrl.field_id,
			fp, dofs_write);

	// phase tranformation field
	fprintf(fp, "SCALARS phase double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	dofs_write[0] = 1;
	dofs_write[1] = 0; //
	utr_write_paraview_field(pdv_heat_problem.ctrl.phase_transformation_field_id, fp, dofs_write);

	// austenite field
	fprintf(fp, "SCALARS austenite double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	dofs_write[0] = 1;
	dofs_write[1] = P_HEAT_AUSTENITE; //
	utr_write_paraview_field(pdv_heat_problem.ctrl.phases_field_id, fp, dofs_write);

//   // bainite field
//   fprintf(fp, "SCALARS bainite double 1\n");
//   fprintf(fp,"LOOKUP_TABLE default\n");
//
//   dofs_write[0]=1;
//   dofs_write[1]=P_HEAT_BAINITE; //
//   utr_write_paraview_field(pdv_heat_problem.ctrl.phases_field_id, fp, dofs_write);

	// martensite field
	fprintf(fp, "SCALARS martensite double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	dofs_write[0] = 1;
	dofs_write[1] = P_HEAT_MARTENSITE; // for testing only here
	utr_write_paraview_field(pdv_heat_problem.ctrl.phases_field_id, fp, dofs_write);

	// martensite heat source field
	fprintf(fp, "SCALARS martensite_Qv double 1\n");  // for testing only here
	fprintf(fp, "LOOKUP_TABLE default\n");

	dofs_write[0] = 1;
	dofs_write[1] = P_HEAT_PEARLITE; //
	utr_write_paraview_field(pdv_heat_problem.ctrl.phases_field_id, fp, dofs_write);

	// bainite field
	fprintf(fp, "SCALARS bainite double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	dofs_write[0] = 1;
	dofs_write[1] = P_HEAT_BAINITE;
	utr_write_paraview_field(pdv_heat_problem.ctrl.phases_field_id, fp, dofs_write);

	// bainite heat source field
	fprintf(fp, "SCALARS bainite_Qv double 1\n");  // for testing only here
	fprintf(fp, "LOOKUP_TABLE default\n");

	dofs_write[0] = 1;
	dofs_write[1] = P_HEAT_FERRITE; // for testing only here
	utr_write_paraview_field(pdv_heat_problem.ctrl.phases_field_id, fp, dofs_write);

	// temperature dependent coefficients
	dofs_write[0] = 1;
	dofs_write[1] = 0; // temperature is the only component
	pdr_heat_write_temperature_dpendent_fields(pdv_heat_problem.ctrl.field_id, fp, dofs_write);

#endif

	// write latent heat of fusion factor
	//fprintf(fp, "SCALARS latent_heat_factor double 1\n");
	//fprintf(fp, "LOOKUP_TABLE default\n");
	//dofs_write[0]=1;
	//dofs_write[1]=0; // fusion factor is the only component
	//utr_write_paraview_field(pdv_heat_dtdt_problem.ctrl.field_id,
	//			   fp, dofs_write);

	// write heating-cooling rate
	//fprintf(fp, "SCALARS heating_cooling_rate double 1\n");
	//fprintf(fp, "LOOKUP_TABLE default\n");
	//dofs_write[0]=1;
	//dofs_write[1]=1; // heating-cooling rate is the only component
	//utr_write_paraview_field(pdv_heat_dtdt_problem.ctrl.field_id,
	//			   fp, dofs_write);


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

