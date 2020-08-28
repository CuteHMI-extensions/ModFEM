import qbs

import cutehmi

Project {
	name: "ModFEM.NSSUPGHeatUI.2"

	//	condition: false

	cutehmi.Extension {
		name: parent.name

		friendlyName: "Heat User Interface"

		vendor: "ModFEM"

		domain: "modfem.agh.edu.pl"

		description: "Heat user interface."

		files: [
         "ArrowEntity.qml",
         "Axes.qml",
         "BCColorMappers.qml",
         "CentrifugalFanEntity.qml",
         "ClipPlaneEntity.qml",
         "Console.qml",
         "ControlsGroup.qml",
         "CustomPhongAlphaMaterial.qml",
         "ElementSelectionGroup.qml",
         "ElementsEntity.qml",
         "ItemEntity.qml",
         "MeshInfo.qml",
		 "NumberDisplayEntity.qml",
         "PhongAlphaClipMaterial.qml",
         "PressureColorMappers.qml",
         "ProbeEntity.qml",
         "ProblemGroup.qml",
         "ProblemInfo.qml",
         "Project.qml",
         "RootEntity.qml",
         "RotationGroup.qml",
         "ScaleGroup.qml",
         "Scene3D.qml",
         "SimulationGroup.qml",
         "SimulationInfo.qml",
         "SolutionInfo.qml",
         "SurfaceEntity.qml",
         "TemperatureColorMappers.qml",
         "TestEntity.qml",
         "TranslationGroup.qml",
         "TransparentTextureMaterial.qml",
         "VelocityColorMappers.qml",
         "VertexAlphaClipMaterial.qml",
         "View.qml",
         "VisibilityGroup.qml",
         "dev/ModFEM.NSSUPGHeatUI-1.workaround.nvidia.bug.txt",
         "images/keyboard-black-18dp.svg",
         "images/keyboard-white-18dp.svg",
         "shaders/es2/transparenttexture.frag",
         "shaders/es2/transparenttexture.vert",
         "shaders/gl3/pervertexcolorclip.frag",
         "shaders/gl3/pervertexcolorclip.vert",
         "shaders/gl3/phongalpha.frag",
         "shaders/gl3/phongalpha.vert",
         "shaders/gl3/phongalphaclip.frag",
         "shaders/gl3/phongalphaclip.vert",
         "shaders/gl3/phongclip.frag",
         "shaders/gl3/phongclip.vert",
         "shaders/gl3/transparenttexture.frag",
         "shaders/gl3/transparenttexture.vert",
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

		FileTagger {
			patterns: "*.fbx"
			fileTags: ["fbx"]
		}

		FileTagger {
			patterns: "*.obj"
			fileTags: ["obj"]
		}

		FileTagger {
			patterns: "*.png"
			fileTags: ["png"]
		}

		FileTagger {
			patterns: "*.vert"
			fileTags: ["vert"]
		}

		FileTagger {
			patterns: "*.frag"
			fileTags: ["frag"]
		}

		Group {
			name: "3D assets"
			fileTagsFilter: ["fbx", "obj", "png", "vert", "frag"]
			qbs.install: true
			qbs.installSourceBase: installSourceBase
			qbs.installDir: dedicatedInstallSubdir
		}
	}
}

//(c)C: Copyright © 2020, Michał Policht <michal@policht.pl>. All rights reserved.
//(c)C: This file is a part of CuteHMI.
//(c)C: CuteHMI is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//(c)C: CuteHMI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//(c)C: You should have received a copy of the GNU Lesser General Public License along with CuteHMI.  If not, see <https://www.gnu.org/licenses/>.
