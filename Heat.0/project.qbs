import qbs

import cutehmi

Project {
	name: "ModFEM.Heat.0"

	cutehmi.CppExtension {
		name: parent.name

		friendlyName: "Heat Qt interface"

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		description: "Heat user interface."

		files: [
         "include/modfem/heat/BoundaryConditionsData.hpp",
         "include/modfem/heat/ElementData.hpp",
         "include/modfem/heat/Problem.hpp",
         "include/modfem/heat/VertexColorMapper.hpp",
         "include/modfem/heat/internal/common.hpp",
         "include/modfem/heat/internal/platform.hpp",
         "include/modfem/heat/logging.hpp",
         "include/modfem/heat/metadata.hpp",
         "src/modfem/heat/BoundaryConditionsData.cpp",
         "src/modfem/heat/ElementData.cpp",
         "src/modfem/heat/Problem.cpp",
         "src/modfem/heat/VertexColorMapper.cpp",
         "src/modfem/heat/internal/QMLPlugin.cpp",
         "src/modfem/heat/internal/QMLPlugin.hpp",
         "src/modfem/heat/logging.cpp",
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

		Depends { name: "ModFEM.pd_heat.2" }

		Export {
			Depends { name: "CuteHMI.2" }

			Depends { name: "ModFEM.pd_heat.2" }
		}
	}
}

//(c)C: Copyright © 2020, Michał Policht <michal@policht.pl>. All rights reserved.
//(c)C: This file is a part of CuteHMI.
//(c)C: CuteHMI is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//(c)C: CuteHMI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//(c)C: You should have received a copy of the GNU Lesser General Public License along with CuteHMI.  If not, see <https://www.gnu.org/licenses/>.
