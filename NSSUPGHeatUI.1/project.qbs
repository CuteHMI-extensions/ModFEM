import qbs

import cutehmi

Project {
	name: "ModFEM.NSSUPGHeatUI.1"

	//	condition: false

	cutehmi.Extension {
		name: parent.name

		friendlyName: "Heat User Interface"

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		description: "Heat user interface."

		files: [
         "BCColorMappers.qml",
         "ElementsEntity.qml",
         "MeshInfo.qml",
         "ProblemInfo.qml",
         "Scene3D.qml",
         "SurfaceEntity.qml",
         "TemperatureColorMappers.qml",
         "TestEntity.qml",
         "View.qml",
     ]

		Depends { name: "CuteHMI.2" }

		Depends { name: "cutehmi.doxygen" }
		cutehmi.doxygen.warnIfUndocumented: false
		cutehmi.doxygen.useDoxyqml: true
		cutehmi.doxygen.exclude: ['dev', 'tests']

		Depends { name: "cutehmi.metadata" }

		Depends { name: "cutehmi.qmldir" }

		Depends { name: "cutehmi.qmltypes" }

		Export {
			Depends { name: "CuteHMI.2" }
		}

//		FileTagger {
//			patterns: "*.fbx"
//			fileTags: ["fbx"]
//		}

//		FileTagger {
//			patterns: "*.obj"
//			fileTags: ["obj"]
//		}

//		FileTagger {
//			patterns: "*.png"
//			fileTags: ["png"]
//		}

//		Group {
//			name: "3D assets"
//			fileTagsFilter: ["fbx", "obj", "png"]
//			qbs.install: true
//			qbs.installSourceBase: installSourceBase
//			qbs.installDir: dedicatedInstallSubdir
//		}
	}
}

//(c)C: Copyright © 2020, Michał Policht <michal@policht.pl>. All rights reserved.
//(c)C: This file is a part of CuteHMI.
//(c)C: CuteHMI is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//(c)C: CuteHMI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//(c)C: You should have received a copy of the GNU Lesser General Public License along with CuteHMI.  If not, see <https://www.gnu.org/licenses/>.
