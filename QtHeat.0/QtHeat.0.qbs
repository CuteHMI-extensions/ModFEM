import qbs

import cutehmi

Project {
	name: "ModFEM.QtHeat.0"

	cutehmi.CppExtension {
		name: parent.name

		friendlyName: "Heat Qt interface"

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		description: "Heat user interface."

		files: [
         "include/modfem/qtheat/ElementData.hpp",
         "include/modfem/qtheat/FaceData.hpp",
         "include/modfem/qtheat/Mesh.hpp",
         "include/modfem/qtheat/Problem.hpp",
         "include/modfem/qtheat/internal/common.hpp",
         "include/modfem/qtheat/internal/platform.hpp",
         "include/modfem/qtheat/logging.hpp",
         "include/modfem/qtheat/metadata.hpp",
         "src/modfem/qtheat/ElementData.cpp",
         "src/modfem/qtheat/FaceData.cpp",
         "src/modfem/qtheat/Mesh.cpp",
         "src/modfem/qtheat/Problem.cpp",
         "src/modfem/qtheat/internal/QMLPlugin.cpp",
         "src/modfem/qtheat/internal/QMLPlugin.hpp",
         "src/modfem/qtheat/logging.cpp",
     ]

		Depends { name: "Qt.3drender" }

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
