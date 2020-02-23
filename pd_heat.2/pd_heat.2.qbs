import qbs

import cutehmi

Project {
	name: "ModFEM.pd_heat.2"

	cutehmi.CppExtension {
		name: parent.name

		minor: 0

		micro: 0

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		friendlyName: "Heat Problem"

		description: "Heat problem."

		files: [
			"include/modfem/pd_heat/pdh_heat.h",
			"include/modfem/pd_heat/pdh_heat_bc.h",
			"include/modfem/pd_heat/pdh_heat_materials.h",
			"include/modfem/pd_heat/pdh_heat_problem.h",
			"include/modfem/pd_heat/pdh_heat_weakform.h",
			"src/modfem/pd_heat/adaptation/pds_heat_adapt.c",
			"src/modfem/pd_heat/input_output/pds_heat_bc_io.c",
			"src/modfem/pd_heat/input_output/pds_heat_io.c",
			"src/modfem/pd_heat/input_output/pds_heat_materials_io.c",
			"src/modfem/pd_heat/input_output/pds_heat_paraview_io.c",
			"src/modfem/pd_heat/input_output/pds_heat_problem_io.c",
			"src/modfem/pd_heat/linear_solver_interface/pds_heat_ls_intf.c",
			"src/modfem/pd_heat/materials/pds_heat_materials.c",
			"src/modfem/pd_heat/time_integration/pds_heat_time_integration.c",
			"src/modfem/pd_heat/weak_formulation/pdh_heat_weakform.history",
			"src/modfem/pd_heat/weak_formulation/pds_heat_bc.c",
			"src/modfem/pd_heat/weak_formulation/pds_heat_weakform.c",
		]

		Group {
			name: "allocation"

//			condition: modfem.config.pd === "heat"

			files: [
				"src/modfem/pd_heat/pds_heat.c",
			]
		}

		Depends { name: "ModFEM.2" }

		Export {
			Depends { name: "ModFEM.Interfaces.2" }
		}
	}
}
