/************************************************************************
File pds_heat_materials.c- definition of functions related to materials
			   handling

Contains definition of routines:
  pdr_heat_material_query - gets material data

------------------------------
History:
	2011    - Przemyslaw Plaszewski (pplaszew@agh.edu.pl)
	2011    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
	2016    - Aleksander Siwek (Aleksander.Siwek@agh.edu.pl)
*************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* problem module's types and functions */
#include <modfem/pd_heat/pdh_heat.h>		/* USES */
/* types and functions related to materials handling */
#include <modfem/pd_heat/pdh_heat_materials.h>	/* IMPLEMENTS */
//#include <modfem/uth_log.h>
/* interface for all approximation modules */
#include <modfem/aph_intf.h>		/* USES */
/* interface for all mesh manipulation modules */
#include <modfem/mmh_intf.h>		/* USES */
/* from problem dependent module */
#include <modfem/pdh_control_intf.h>	/* IMPLEMENTS */


/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */

/*------------------------------------------------------------
pdr_heat_material_query - gets material data
------------------------------------------------------------*/
int pdr_heat_material_query(const utt_material_query_params * Params,
		utt_material_query_result * Result)
{
	utr_mat_init_query_result(Result);

	Result->density = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->thermal_conductivity = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->specific_heat = UTC_MAT_QUERY_RESULT_REQUIRED;

#ifdef PHASE_TRANSFORMATION

	Result->T_Fs = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->T_Ff = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->T_Ps = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->T_Pf = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->T_Bs = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->T_Bf = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->T_Ms = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->T_Mf = UTC_MAT_QUERY_RESULT_REQUIRED;

	Result->M_Vc = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->B_Vc = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->F_Vc = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->P_Vc = UTC_MAT_QUERY_RESULT_REQUIRED;

	Result->H_A_F = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->H_A_P = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->H_A_B = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->H_A_M = UTC_MAT_QUERY_RESULT_REQUIRED;

	Result->tB_lim = UTC_MAT_QUERY_RESULT_REQUIRED;
	Result->tM_lim = UTC_MAT_QUERY_RESULT_REQUIRED;

#endif

	/*cin
	  printf("\nParams.group_idx = %lf",Params->group_idx);
	/*cot*/

	return utr_material_query(Params, Result);
}

#ifdef PHASE_TRANSFORMATION
/**-----------------------------------------------------------
pdr_phases_field_init - initialize phases field
------------------------------------------------------------*/
int pdr_phases_field_init(int Problem_id)
{
	pdt_heat_problem * problem = (pdt_heat_problem *)pdr_get_problem_structure(Problem_id);
	int heat_field_id = problem->ctrl.field_id;
	int p_field_id = problem->ctrl.phases_field_id;
	int p_t_field_id = problem->ctrl.phase_transformation_field_id;
	int mesh_id = problem->ctrl.mesh_id;
	double p_current[APC_MAXELSD];
	double p_t_current[APC_MAXELSD];
	int p_numdofs, p_t_numdofs;
	int node_id, i;

	/*cin*/
	printf("\nheat_field_id = %d\tp_field_id = %d", heat_field_id, p_field_id);
	/*cot*/
	node_id = 0;
	while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
		if (apr_get_ent_pdeg(heat_field_id, APC_VERTEX, node_id) > 0) {
			p_current[P_HEAT_AUSTENITE] = 100.0;	// gear entirely austenite before cooling
			p_current[P_HEAT_BAINITE] = 0.0;
			p_current[P_HEAT_FERRITE] = 0.0;
			p_current[P_HEAT_MARTENSITE] = 0.0;
			p_current[P_HEAT_PEARLITE] = 0.0;
			p_numdofs = apr_get_ent_nrdofs(p_field_id, APC_VERTEX, node_id);
			apr_write_ent_dofs(p_field_id, APC_VERTEX, node_id, p_numdofs, Current_solution_ID, p_current);

			p_t_current[0] = PT_HEAT_UNKNOWN;
			p_t_current[1] = 0.0;
			p_t_current[2] = 0.0;
			p_t_current[3] = 0.0;
			p_t_current[4] = 0.0;	// dT -> cooling rate V_{8/5} = dT/dt
			p_t_current[5] = 0.0;	// dt -> cooling rate V_{8/5} = dT/dt
			p_t_current[6] = 0.0;	// V_{8/5} -> cooling rate V_{8/5} = dT/dt
			p_t_numdofs = apr_get_ent_nrdofs(p_t_field_id, APC_VERTEX, node_id);
			apr_write_ent_dofs(p_t_field_id, APC_VERTEX, node_id, p_t_numdofs, Current_solution_ID, p_t_current);
			/*cin
				  printf("\nnode = %d", node_id);
				  printf("\t[AUS] = %lf\t[FER] = %lf\t[PEA] = %lf\t[BAI] = %lf\t[MAR] = %lf",
				  p_current[P_HEAT_AUSTENITE], p_current[P_HEAT_FERRITE], p_current[P_HEAT_PEARLITE],
				  p_current[P_HEAT_BAINITE], p_current[P_HEAT_MARTENSITE]);
			/*cot*/
		}
	}
}

/**-----------------------------------------------------------
pdr_heat_material_phase_transformation - transforms phases
------------------------------------------------------------*/
int pdr_phases_transformation(int Problem_id)
{
	int node_id;
	int Current, Previous;
	double p_current[APC_MAXELSD], p_previous[APC_MAXELSD];	/* solution dofs */
	double p_t_current[APC_MAXELSD], p_t_previous[APC_MAXELSD];	/* solution dofs */
	double T_current[APC_MAXELSD], T_previous[APC_MAXELSD];	/* solution dofs */
	int T_numdofs, p_numdofs, p_t_numdofs, i;
	pdt_heat_problem * problem = (pdt_heat_problem *)pdr_get_problem_structure(Problem_id);
	int heat_field_id = problem->ctrl.field_id;
	int p_field_id = problem->ctrl.phases_field_id;
	int p_t_field_id = problem->ctrl.phase_transformation_field_id;
	int mesh_id = problem->ctrl.mesh_id;
	utt_material_query_params query_params;
	utt_material_query_result query_result;
	double Vc;
	double k_M = 0.01;	// martensitic transformation constant k_M=-ln(0.01)/(T_Ms-T_Mf) <- Vc=const
	double b_B, n_B;		// bainitic transformation constants n_B=6.12733/ln(t_Bs/t_Bf), b_B=0.01005/(t_Bf^n_B) <- Vc=const
	double ts, ts_v, tf, tf_v, dt, dts, dtf;

	query_params.group_idx = mmr_el_groupID(mesh_id, 1);	// ? node belongs to multiple elements
	query_params.name = "";	// ... not by material name
	query_params.time = problem->time.cur_time;

	apr_rewr_sol(p_field_id, Current_solution_ID, Previous_time_step_sol_ID);
	apr_rewr_sol(p_t_field_id, Current_solution_ID, Previous_time_step_sol_ID);

	node_id = 0;

	dt = problem->time.cur_dtime;
	node_id = 0;
	while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
		if (apr_get_ent_pdeg(heat_field_id, APC_VERTEX, node_id) > 0) {
			T_numdofs = apr_get_ent_nrdofs(heat_field_id, APC_VERTEX, node_id);
			apr_read_ent_dofs(heat_field_id, APC_VERTEX, node_id, T_numdofs, Current_solution_ID, T_current);
			query_params.temperature = T_current[0];
			apr_read_ent_dofs(heat_field_id, APC_VERTEX, node_id, T_numdofs, Previous_time_step_sol_ID, T_previous);
			Vc = (T_current[0] - T_previous[0]) / dt;

			p_t_numdofs = apr_get_ent_nrdofs(p_t_field_id, APC_VERTEX, node_id);
			apr_read_ent_dofs(p_t_field_id, APC_VERTEX, node_id, p_t_numdofs, Current_solution_ID, p_t_current);
			apr_read_ent_dofs(p_t_field_id, APC_VERTEX, node_id, p_t_numdofs, Previous_time_step_sol_ID, p_t_previous);

			query_params.Vc = fabs(Vc);
			// for high cooling rate, node jumps below {8/5} range, in one time step
			if ( (p_t_current[0] == PT_HEAT_UNKNOWN) && T_current[0] <= 1073.0 && T_current[0] >= 773.0 ) {
				p_t_current[4] += (T_current[0] - T_previous[0]);
				p_t_current[5] += dt;
				p_t_current[6] = p_t_current[4] / p_t_current[5];
				query_params.Vc = fabs(p_t_current[6]);
				printf("");
			}
			pdr_heat_material_query(&query_params, &query_result);

			p_numdofs = apr_get_ent_nrdofs(p_field_id, APC_VERTEX, node_id);
			apr_read_ent_dofs(p_field_id, APC_VERTEX, node_id, p_numdofs, Current_solution_ID, p_current);
			apr_read_ent_dofs(p_field_id, APC_VERTEX, node_id, p_numdofs, Previous_time_step_sol_ID, p_previous);

			if ( Vc < 0.0 ) {	// progress of transformation only when cooling, current cooling rate Vc<0.0
				if ( (p_t_current[0] == PT_HEAT_UNKNOWN) &&
						(query_params.time < query_result.tB_lim) &&	// only martensitic transf. possible before that time (<-CCT data)
						(query_params.temperature <= query_result.T_Ms) ) {
					p_t_current[0] = PT_HEAT_MARTENSITE_START;
					p_t_current[1] = query_result.M_Vc;
					p_t_current[2] = query_result.T_Ms;
					p_t_current[3] = p_current[P_HEAT_AUSTENITE];
				} else if ( (p_t_current[0] == PT_HEAT_UNKNOWN) &&
						(query_params.time >= query_result.tB_lim) &&	// bainitic transformation possible only after that time (<-CCT data)
						(query_params.temperature <= query_result.T_Bs) &&
						(query_params.temperature >= query_result.T_Bf) ) {
					p_t_current[0] = PT_HEAT_BAINITE_START;
					dts = (query_result.T_Bs - query_params.temperature) / fabs(Vc);
					ts_v = problem->time.cur_time - dts;	// approximation
					dtf = (query_params.temperature - query_result.T_Bf) / fabs(Vc);	// for current cooling rate
					tf_v = problem->time.cur_time + dtf;
// 		  tf_v = ts_v + (query_params.temperature - query_result.T_Bf )/fabs(p_t_current[6]);	// for averaged cooling rate
					if ( (ts_v > query_result.tB_lim) && (tf_v < query_result.tM_lim) ) {
						ts = ts_v;
						tf = tf_v;	// T_Bf ~const for time t<tM_lim[s]
						p_t_current[1] = log( log(0.99) / log(0.01) ) / log(ts / tf);	// n
						p_t_current[2] = -log(0.01) / ( pow(tf, p_t_current[1]) );	// b
						if ( p_current[P_HEAT_AUSTENITE] < query_result.B_Vc ) {
							p_t_current[3] = p_current[P_HEAT_AUSTENITE];
						} else {
							p_t_current[3] = query_result.B_Vc;
						}
						printf("");
					} else {
						p_t_current[1] = p_t_previous[1];
						p_t_current[2] = p_t_previous[2];
						p_t_current[3] = p_t_previous[3];
						printf("");
					}
				} else if ( (p_t_current[0] == PT_HEAT_MARTENSITE_START) &&
						(query_params.temperature <= query_result.T_Mf) ) {
					p_t_current[0] = PT_HEAT_FINISHED;
					// let last part of Austenite transform to Martensite befor finish
					p_current[P_HEAT_MARTENSITE] = p_t_previous[3] * ( 1.0 - exp(-k_M * ( p_t_previous[2] - query_result.T_Mf ) ) );
					p_current[P_HEAT_AUSTENITE] -= (p_current[P_HEAT_MARTENSITE] - p_previous[P_HEAT_MARTENSITE]);
// 		  p_current[P_HEAT_AUSTENITE] = p_t_current[3] - p_current[P_HEAT_MARTENSITE];
					dtf = (query_result.T_Mf - query_params.temperature) / fabs(Vc);
					p_current[P_HEAT_PEARLITE] = 0.01 * ( p_current[P_HEAT_MARTENSITE] - p_previous[P_HEAT_MARTENSITE] ) * query_result.H_A_M * query_result.density / (dt - dtf); //problem->time.cur_dtime;
					p_current[P_HEAT_FERRITE] = 0.0;
					assert(p_current[P_HEAT_AUSTENITE]  > -10e200 && p_current[P_HEAT_AUSTENITE]  < 10e200);
					assert(p_current[P_HEAT_MARTENSITE] > -10e200 && p_current[P_HEAT_MARTENSITE] < 10e200);
#ifdef DEBUG
					if (node_id == 686 || node_id == 2066 || node_id == 1) {
						printf("\n\t\tPT_HEAT_MARTENSITE_START >>>> node_id = %d\tV_c = %lf\tquery Vc =%lf", node_id, Vc, query_params.Vc);
						printf("\n\t\tPT_HEAT_MARTENSITE_START >>>> k = %lf\tMs = %lf\tMf = %lf\tT =%lf", k_M, query_result.T_Ms, query_result.T_Mf, query_params.temperature);
						printf("\n\t\tPT_HEAT_MARTENSITE_START >>>> A = %lf\tM = %lf\n", p_current[P_HEAT_AUSTENITE], p_current[P_HEAT_MARTENSITE]);
						printf("\n");
					}
#endif


				} else if ( (p_t_current[0] == PT_HEAT_BAINITE_START) &&
						(query_params.time < query_result.tM_lim) &&			// martensitic transformation possible only befor that time (<-CCT data)
						(query_params.temperature <= query_result.T_Ms) ) {		// !!! model is true when T_Bf = T_Ms
					p_t_current[0] = PT_HEAT_MARTENSITE_START;
					// let last part of Austenite transform to Bainite, befor Martensite start
					dtf = (query_result.T_Bf - query_params.temperature) / fabs(Vc);
					p_current[P_HEAT_BAINITE] = p_t_previous[3] * ( 1.0 - exp(-p_t_previous[2] * pow((query_params.time - dtf), p_t_previous[1]) ) );
					p_current[P_HEAT_AUSTENITE] -= (p_current[P_HEAT_BAINITE] - p_previous[P_HEAT_BAINITE]);
// 			p_current[P_HEAT_AUSTENITE] = p_t_current[3] - p_current[P_HEAT_BAINITE];
					p_current[P_HEAT_FERRITE] = 0.01 * ( p_current[P_HEAT_BAINITE] - p_previous[P_HEAT_BAINITE] ) * query_result.H_A_B * query_result.density / dt; //problem->time.cur_dtime;
					p_current[P_HEAT_PEARLITE] = 0.0;

					p_t_current[1] = query_result.M_Vc;
					p_t_current[2] = query_result.T_Ms;
					p_t_current[3] = p_current[P_HEAT_AUSTENITE];

				}
// 		else => node jumps below T_Mf in first time step

				switch ( (int)p_t_current[0] ) {
					case PT_HEAT_UNKNOWN:
						// check if ferritic, pearlitic, bainitic or martensitic transformation start
						break;
					case PT_HEAT_FERRITE_START:
						// check if ferritic transformation stop
						break;
					case PT_HEAT_FERRITE_STOP:
						// check if pearlitic transformation start
						break;
					case PT_HEAT_PEARLITE_START:
						// check if pearlitic transformation stop
						break;
					case PT_HEAT_PEARLITE_STOP:
						// check if bainitic transformation start
						break;
					case PT_HEAT_BAINITE_START:
						// check if bainitic transformation stop
						p_current[P_HEAT_BAINITE] = p_t_current[3] * ( 1.0 - exp(-p_t_current[2] * pow(query_params.time, p_t_current[1]) ) );
						p_current[P_HEAT_AUSTENITE] -= (p_current[P_HEAT_BAINITE] - p_previous[P_HEAT_BAINITE]);
// 			p_current[P_HEAT_AUSTENITE] = p_t_current[3] - p_current[P_HEAT_BAINITE];
						p_current[P_HEAT_FERRITE] = 0.01 * ( p_current[P_HEAT_BAINITE] - p_previous[P_HEAT_BAINITE] ) * query_result.H_A_B * query_result.density / dt; //problem->time.cur_dtime;
						p_current[P_HEAT_PEARLITE] = 0.0;
						assert(p_current[P_HEAT_AUSTENITE] > -10e200 && p_current[P_HEAT_AUSTENITE] < 10e200);
						assert(p_current[P_HEAT_BAINITE]   > -10e200 && p_current[P_HEAT_BAINITE]   < 10e200);
#ifdef DEBUG
						if (node_id == 686 || node_id == 2066 || node_id == 1) {
							printf("\n\t\tPT_HEAT_BAINITE_START >>>> node_id = %d\tV_c = %lf\tquery Vc =%lf", node_id, Vc, query_params.Vc);
							printf("\n\t\tPT_HEAT_BAINITE_START >>>> b = % lf\tn = %lf\tBs = %lf\tBf = %lf\tT =%lf", p_t_current[2], p_t_current[1], query_result.T_Bs, query_result.T_Bf, query_params.temperature);
							printf("\n\t\tPT_HEAT_BAINITE_START >>>> A = %lf\tB = %lf\n", p_current[P_HEAT_AUSTENITE], p_current[P_HEAT_BAINITE]);
							printf("\n");
						}
#endif
						break;
					case PT_HEAT_BAINITE_STOP:
						break;
					case PT_HEAT_MARTENSITE_START:
						// check if martensitic transformation start
						p_current[P_HEAT_MARTENSITE] = p_t_current[3] * ( 1.0 - exp(-k_M * ( p_t_current[2] - query_params.temperature ) ) );
						p_current[P_HEAT_AUSTENITE] -= (p_current[P_HEAT_MARTENSITE] - p_previous[P_HEAT_MARTENSITE]);
// 			p_current[P_HEAT_AUSTENITE] = p_t_current[3] - p_current[P_HEAT_MARTENSITE];
						p_current[P_HEAT_PEARLITE] = 0.01 * ( p_current[P_HEAT_MARTENSITE] - p_previous[P_HEAT_MARTENSITE] ) * query_result.H_A_M * query_result.density / dt; //problem->time.cur_dtime;
						p_current[P_HEAT_FERRITE] = 0.0;
						assert(p_current[P_HEAT_AUSTENITE]  > -10e200 && p_current[P_HEAT_AUSTENITE]  < 10e200);
						assert(p_current[P_HEAT_MARTENSITE] > -10e200 && p_current[P_HEAT_MARTENSITE] < 10e200);
#ifdef DEBUG
						if (node_id == 686 || node_id == 2066 || node_id == 1) {
							printf("\n\t\tPT_HEAT_MARTENSITE_START >>>> node_id = %d\tV_c = %lf\tquery Vc =%lf", node_id, Vc, query_params.Vc);
							printf("\n\t\tPT_HEAT_MARTENSITE_START >>>> k = %lf\tMs = %lf\tMf = %lf\tT =%lf", k_M, query_result.T_Ms, query_result.T_Mf, query_params.temperature);
							printf("\n\t\tPT_HEAT_MARTENSITE_START >>>> A = %lf\tM = %lf\n", p_current[P_HEAT_AUSTENITE], p_current[P_HEAT_MARTENSITE]);
							printf("\n");
						}
#endif
						break;
					case PT_HEAT_MARTENSITE_STOP:
						// check if subquenching start
						break;
					case PT_HEAT_FINISHED:
						p_current[P_HEAT_FERRITE] = 0.0;
						p_current[P_HEAT_PEARLITE] = 0.0;
						// check if subquenching start
						break;
					default:
						break;
				}
				apr_write_ent_dofs(p_t_field_id, APC_VERTEX, node_id, p_t_numdofs, Current_solution_ID, p_t_current);
				apr_write_ent_dofs(p_field_id, APC_VERTEX, node_id, p_numdofs, Current_solution_ID, p_current);
			}
		}
	}

}
#endif
