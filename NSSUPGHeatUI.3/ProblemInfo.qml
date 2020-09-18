import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import ModFEM.NSSUPGHeat 1.0

GroupBox {
	title: qsTr("Problem")

	property Problem problem

	ColumnLayout {
		Label {
			text: qsTr("Problem id: ") + problem.problemId
		}

		Label {
			text: qsTr("Mesh id: ") + problem.meshId
		}

		Label {
			text: qsTr("Field id: ") + problem.fieldId
		}

		Label {
			text: qsTr("Solution count: ") + problem.solutionCount
		}

		Label {
			text: qsTr("Equation count: ") + problem.equationCount
		}
	}
}
