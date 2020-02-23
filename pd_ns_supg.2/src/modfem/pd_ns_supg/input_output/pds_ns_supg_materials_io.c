/************************************************************************
File pds_ns_supg_materials_io.c - definition of functions related to materials
			   read and management.

Contains definition of routines:
  pdr_ns_supg_material_read - read materials data from config file
  pdr_ns_supg_material_free - free materials resources

------------------------------
History:
	2011    - Przemyslaw Plaszewski (pplaszew@agh.edu.pl)
	2011    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libconfig.h>

#include <modfem/pd_ns_supg/pdh_ns_supg.h> /* USES */
/* types and functions related to materials handling */
#include <modfem/pd_ns_supg/pdh_ns_supg_materials.h> /* IMPLEMENTS */
#include <modfem/uth_mat.h>
#include <modfem/uth_log.h>

/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */

///*------------------------------------------------------------
//pdr_ns_supg_material_free - free materials resources
//------------------------------------------------------------*/
//int pdr_ns_supg_material_free(pdt_ns_supg_materials *Materials_db)
//{
//  return 0; //TODO
//}

/*------------------------------------------------------------
  pdr_material_read - read materials data from config file
------------------------------------------------------------*/
int pdr_ns_supg_material_read(char *Work_dir,
				  char *Filename,
				  FILE *Interactive_output)
{

  ute_mat_read_result result = utr_mat_read(Work_dir,Filename,Interactive_output);

  switch(result) {
  case UTE_FAIL:
	break;
  case UTE_REF_TEMP_MUST_BE_GEQ_ZERO:
	if(pdv_ns_supg_problem.ctrl.ref_temperature <= 0.0){
	mf_fatal_err("ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
	//fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
	  //exit(-1);
	  //return(EXIT_FAILURE);
	}
	break;
  case UTE_NOT_NEED_MATERIAL_DATABASE:
	// we do not need material database
	// just read some values from defualt material into problem structures
	if( utr_mat_get_n_materials() > 0 ) {

		const int n_mats = utr_mat_get_n_materials();
		int* matsIDs = (int*) calloc(n_mats,sizeof(int));
		if(utr_mat_get_materials_IDs(matsIDs) < 1) {
			mf_fatal_err("Error in reading materials IDs!");
		}

		pdv_ns_supg_problem.ctrl.ref_temperature = -1.0;
		pdv_ns_supg_problem.ctrl.dynamic_viscosity =
		  utr_mat_get_material_by_matID(matsIDs[0])->atT_dynamic_viscosity[0];
		pdv_ns_supg_problem.ctrl.density =
		  utr_mat_get_material_by_matID(matsIDs[0])->atT_density[0];
		free(matsIDs);
		result = UTE_SUCCESS;
	}
	else {
		mf_fatal_err("There is NOT any material, so values atT_dynamic_viscosity and atT_density are undefined!");
	}
	break;
  case UTE_SUCCESS:
  default:
	break;
  }

  return result;
}


