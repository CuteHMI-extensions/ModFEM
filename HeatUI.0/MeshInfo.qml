import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import ModFEM.Heat 0.0

ColumnLayout {
	property Problem problem


	GroupBox {
		Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

		title: qsTr("Total")

		ColumnLayout {
			Label {
				text: qsTr("Elements: ") + problem.elements.count.elements
			}

		}
	}

	GroupBox {
		Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

		title: qsTr("Display")

		ColumnLayout {
			Label {
				text: qsTr("Nodes: ") + problem.elements.nodes.count
			}

			Label {
				text: qsTr("Triangles: ") + problem.elements.triangles.count
			}

			Label {
				text: qsTr("Quads: ") + problem.elements.quads.count
			}

			Label {
				text: qsTr("Lines: ") + problem.elements.lines.count
			}
		}
	}
}
