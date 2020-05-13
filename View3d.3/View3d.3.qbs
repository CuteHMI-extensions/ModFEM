import qbs

import cutehmi

Project {
	name: "ModFEM.View3d.3"

	//	condition: false

	cutehmi.Extension {
		name: parent.name

		friendlyName: "View 3d"

		vendor: "CuteHMI"

		domain: "cutehmi.kde.org"

		description: "3d View."

		files: [
         "AnimatedEntity.qml",
         "CustomMaterial.qml",
         "LogoControls.qml",
         "Main.qml",
         "PlaneEntity.qml",
         "Qt_logo.obj",
		 "WireframeEffect.qml",
         "WireframeMaterial.qml",
         "awesomeface.png",
         "cube.obj",
         "cube_ascii.fbx",
         "man.png",
         "shaders/es2/transparentTexture.frag",
         "shaders/es2/transparentTexture.vert",
         "shaders/gl3/transparentTexture.frag",
         "shaders/gl3/transparentTexture.vert",
         "shaders/texture.frag",
         "shaders/texture.vert",
         "shaders/transparentTexture.frag",
         "shaders/transparentTexture.vert",
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
