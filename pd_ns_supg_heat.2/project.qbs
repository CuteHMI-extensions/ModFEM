import qbs

import cutehmi

Project {
	name: "ModFEM.pd_ns_supg_heat.2"

	cutehmi.CppExtension {
		name: parent.name

//		condition: false

		minor: 0

		micro: 0

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		friendlyName: "Navier Stokes SUPG with heat"

		description: "Navier Stokes flow simulation using streamline upwind Petrov Galerkin method with heat transfer."

		files: [
			"include/modfem/pd_ns_supg_heat/pdh_ns_supg_heat.h",
			"src/modfem/pd_ns_supg_heat/adaptation/pds_ns_supg_heat_adapt.c",
			"src/modfem/pd_ns_supg_heat/input_output/pds_ns_supg_heat_io.c",
			"src/modfem/pd_ns_supg_heat/input_output/pds_ns_supg_heat_paraview_io.c",
			"src/modfem/pd_ns_supg_heat/linear_solver_interface/pds_ns_supg_heat_ls_intf.c",
			"src/modfem/pd_ns_supg_heat/time_integration/pds_ns_supg_heat_time_integration.c",
		]

		Group {
			name: "allocation"

//			condition: modfem.config.pd === "ns_supg_heat"

			files: [
				"src/modfem/pd_ns_supg_heat/pds_ns_supg_heat.c",
			]
		}

		Depends { name: "ModFEM.2" }

		Depends { name: "ModFEM.pd_heat.2" }

		Depends { name: "ModFEM.pd_ns_supg.2" }

		Export {
			Depends { name: "ModFEM.Interfaces.2" }
		}
	}
}
