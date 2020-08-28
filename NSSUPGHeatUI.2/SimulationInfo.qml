import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import ModFEM.NSSUPGHeat 1.0

GroupBox {
	title: qsTr("Simulation")

	property Problem problem

	ColumnLayout {
		Label {
			text: qsTr("Requested time step: ") + problem.integrationData.requestedTimeStep
		}

		Label {
			text: qsTr("Current time step: ") + problem.integrationData.currentTimeStep
		}

		Label {
			text: qsTr("Reall time: ") + problem.integrationData.realTime
		}

		Label {
			text: qsTr("Simulation time: ") + problem.integrationData.simulationTime
		}
	}
}
