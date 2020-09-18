import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3

import ModFEM.NSSUPGHeat 1.0

GroupBox {
	title: qsTr("Simulation")

	property Problem problem

	RowLayout {
		Label {
			text: qsTr("Time step:")
		}

		TextField {
			id: timeStepField

			validator: DoubleValidator{}
			enabled: !realTimeBox.checked
			selectByMouse: true

			onAccepted: problem.integrationData.requestedTimeStep = text

			Connections {
				target: problem.integrationData

				function onRequestedTimeStepChanged() {
					timeStepField.text = problem.integrationData.requestedTimeStep
				}
			}
		}

		CheckBox {
			id: realTimeBox

			text: qsTr("Real time")

			checked: problem.integrationData.realTimeSimulation

			onCheckedChanged: problem.integrationData.realTimeSimulation = checked
		}
	}
}
