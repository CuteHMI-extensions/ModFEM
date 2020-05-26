import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import ModFEM.NSSUPGHeat 1.0

ColumnLayout {
	property ElementData elementData


	GroupBox {
		Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

		title: qsTr("Total")

		ColumnLayout {
			Label {
				text: qsTr("Elements: ") + elementData.count.elements
			}

			Label {
				text: qsTr("Tetrahedrons: ") + elementData.count.tetrahedrons
			}

			Label {
				text: qsTr("Prisms: ") + elementData.count.prisms
			}

			Label {
				text: qsTr("Bricks: ") + elementData.count.bricks
			}

			Label {
				text: qsTr("Triangles: ") + elementData.count.triangles
			}

			Label {
				text: qsTr("Quads: ") + elementData.count.quads
			}

			Label {
				text: qsTr("Faces: ") + elementData.count.faces
			}

			Label {
				text: qsTr("Edges: ") + elementData.count.edges
			}

			Label {
				text: qsTr("Lines: ") + elementData.count.lines
			}

			Label {
				text: qsTr("Nodes: ") + elementData.count.nodes
			}
		}
	}

	GroupBox {
		Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

		title: qsTr("Display")

		ColumnLayout {
			Label {
				text: qsTr("Triangles: ") + elementData.triangles.count
			}

			Label {
				text: qsTr("Quads: ") + elementData.quads.count
			}

			Label {
				text: qsTr("Quad triangles: ") + elementData.quads.triangleCount
			}

			Label {
				text: qsTr("Lines: ") + elementData.lines.count
			}

			Label {
				text: qsTr("Nodes: ") + elementData.nodes.count
			}
		}
	}
}
