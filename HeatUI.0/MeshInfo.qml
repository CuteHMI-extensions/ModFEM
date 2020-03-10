import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import ModFEM.QtHeat 0.0

GroupBox {
	title: qsTr("Mesh")

	property Problem problem

	ColumnLayout {
//		Label {
//			text: qsTr("Nodes: ") + problem.mesh.nodeCount
//		}

//		Label {
//			text: qsTr("Triangle faces: ") + problem.mesh.faceData.triangleCount
//		}

//		Label {
//			text: qsTr("Quad faces: ") + problem.mesh.faceData.quadCount
//		}

		Label {
			text: qsTr("Nodes: ") + problem.mesh.nodes.count
		}

		Label {
			text: qsTr("Triangles: ") + problem.mesh.triangles.count
		}

		Label {
			text: qsTr("Quads: ") + problem.mesh.quads.count
		}

		Label {
			text: qsTr("Lines: ") + problem.mesh.lines.count
		}
	}
}
