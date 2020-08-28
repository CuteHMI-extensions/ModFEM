import qbs

import cutehmi

Project {
	name: "ModFEM.NSSUPGHeat.1"

	cutehmi.CppExtension {
		name: parent.name

		friendlyName: "Heat Qt interface"

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		description: "Heat user interface."

		files: [
         "include/modfem/nssupgheat/AbstractColorMapper.hpp",
         "include/modfem/nssupgheat/AbstractProbe.hpp",
         "include/modfem/nssupgheat/BoundaryConditionsData.hpp",
         "include/modfem/nssupgheat/ElementData.hpp",
         "include/modfem/nssupgheat/IntegrationData.hpp",
         "include/modfem/nssupgheat/IntegrationThread.hpp",
         "include/modfem/nssupgheat/PaletteColorMapper.hpp",
         "include/modfem/nssupgheat/Problem.hpp",
         "include/modfem/nssupgheat/ScalarFieldNodes.hpp",
         "include/modfem/nssupgheat/HueColorMapper.hpp",
         "include/modfem/nssupgheat/ScalarProbe.hpp",
         "include/modfem/nssupgheat/Vector3Probe.hpp",
         "include/modfem/nssupgheat/internal/common.hpp",
         "include/modfem/nssupgheat/internal/platform.hpp",
         "include/modfem/nssupgheat/logging.hpp",
         "include/modfem/nssupgheat/metadata.hpp",
         "src/modfem/nssupgheat/AbstractColorMapper.cpp",
         "src/modfem/nssupgheat/AbstractProbe.cpp",
         "src/modfem/nssupgheat/BoundaryConditionsData.cpp",
         "src/modfem/nssupgheat/ElementData.cpp",
         "src/modfem/nssupgheat/IntegrationData.cpp",
         "src/modfem/nssupgheat/IntegrationThread.cpp",
         "src/modfem/nssupgheat/PaletteColorMapper.cpp",
         "src/modfem/nssupgheat/Problem.cpp",
         "src/modfem/nssupgheat/ScalarFieldNodes.cpp",
         "src/modfem/nssupgheat/HueColorMapper.cpp",
         "src/modfem/nssupgheat/ScalarProbe.cpp",
         "src/modfem/nssupgheat/Vector3Probe.cpp",
         "src/modfem/nssupgheat/internal/QMLPlugin.cpp",
         "src/modfem/nssupgheat/internal/QMLPlugin.hpp",
         "src/modfem/nssupgheat/logging.cpp",
     ]

		Depends { name: "Qt.3drender" }

		Depends { name: "Qt.concurrent" }

		Depends { name: "CuteHMI.2" }

		Depends { name: "cutehmi.doxygen" }
		cutehmi.doxygen.warnIfUndocumented: false
		cutehmi.doxygen.useDoxyqml: true
		cutehmi.doxygen.exclude: ['dev', 'tests']

		Depends { name: "cutehmi.metadata" }

		Depends { name: "cutehmi.qmldir" }

		Depends { name: "cutehmi.qmltypes" }

		Depends { name: "modfem.config" }

		Depends { name: "ModFEM.pd_ns_supg_heat.2" }

		Depends { name: "ModFEM.pd_ns_supg.2"; cpp.link: false }

		Depends { name: "ModFEM.pd_heat.2"; cpp.link: false }

		Export {
			Depends { name: "CuteHMI.2" }

			Depends { name: "ModFEM.pd_ns_supg_heat.2" }

			Depends { name: "ModFEM.pd_ns_supg.2"; cpp.link: false }

			Depends { name: "ModFEM.pd_heat.2"; cpp.link: false }
		}
	}
}

//(c)C: Copyright © 2020, Michał Policht <michal@policht.pl>. All rights reserved.
//(c)C: This file is a part of CuteHMI.
//(c)C: CuteHMI is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//(c)C: CuteHMI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//(c)C: You should have received a copy of the GNU Lesser General Public License along with CuteHMI.  If not, see <https://www.gnu.org/licenses/>.
