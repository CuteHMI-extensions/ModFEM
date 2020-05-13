import qbs

import cutehmi

Project {
	name: "ModFEM.pd_ns_supg.2"

	cutehmi.CppExtension {
		name: parent.name

		type: ["staticlibrary"]

		minor: 0

		micro: 0

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		friendlyName: "Navier Stokes SUPG"

		description: "Navier Stokes flow simulation using streamline upwind Petrov Galerkin method."

		files: [
			"include/modfem/pd_ns_supg/pdh_ns_supg.h",
			"include/modfem/pd_ns_supg/pdh_ns_supg_bc.h",
			"include/modfem/pd_ns_supg/pdh_ns_supg_materials.h",
			"include/modfem/pd_ns_supg/pdh_ns_supg_problem.h",
			"include/modfem/pd_ns_supg/pdh_ns_supg_weakform.h",
			"src/modfem/pd_ns_supg/adaptation/pds_ns_supg_adapt.c",
			"src/modfem/pd_ns_supg/input_output/pds_ns_supg_bc_io.c",
			"src/modfem/pd_ns_supg/input_output/pds_ns_supg_io.c",
			"src/modfem/pd_ns_supg/input_output/pds_ns_supg_materials_io.c",
			"src/modfem/pd_ns_supg/input_output/pds_ns_supg_paraview_io.c",
			"src/modfem/pd_ns_supg/input_output/pds_ns_supg_problem_io.c",
			"src/modfem/pd_ns_supg/linear_solver_interface/pds_ns_supg_ls_intf.c",
			"src/modfem/pd_ns_supg/materials/pds_ns_supg_materials.c",
			"src/modfem/pd_ns_supg/time_integration/pds_ns_supg_time_integration.c",
			"src/modfem/pd_ns_supg/weak_formulation/pdh_ns_supg_weakform.history",
			"src/modfem/pd_ns_supg/weak_formulation/pds_ns_supg_bc.c",
			"src/modfem/pd_ns_supg/weak_formulation/pds_ns_supg_weakform.c",
			"src/modfem/pd_ns_supg/weak_formulation/pds_ns_supg_weakform.c_fully_implicit",
		]

		Group {
			name: "allocation"

//			condition: modfem.config.pd === "ns_supg"

			files: [
				"src/modfem/pd_ns_supg/pds_ns_supg.c",
			]
		}

		Depends { name: "ModFEM.2" }

		Export {
			Depends { name: "ModFEM.Interfaces.2" }
		}
	}
}
