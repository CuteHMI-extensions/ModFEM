import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import ModFEM.NSSUPGHeat 1.0

GroupBox {
	title: qsTr("Solution")

	property Problem problem

	ColumnLayout {
		Label {
			text: qsTr("Max temperature: ") + problem.elementData.maxTemperature
		}

		Label {
			text: qsTr("Min temperature: ") + problem.elementData.minTemperature
		}

		Label {
			text: qsTr("Max pressure: ") + problem.elementData.maxPressure
		}

		Label {
			text: qsTr("Min pressure: ") + problem.elementData.minPressure
		}
	}
}
