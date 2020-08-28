import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3

import CuteHMI.Services 2.0

import ModFEM.NSSUPGHeat 1.0

GroupBox {
	title: qsTr("Controls")

	property Problem problem

	RowLayout {
		Button {
			text: qsTr("Init")

			onClicked: problem.init()
		}

		Button {
			text: qsTr("Start")

			onClicked: {
				problem.start()
				ServiceManager.start()
			}
		}

		Button {
			text: qsTr("Stop")

			onClicked: {
				problem.stop()
				ServiceManager.stop()
			}
		}

		Button {
			text: qsTr("Solve")

			onClicked: problem.solve()
		}

		Button {
			text: qsTr("Update fields")

			onClicked: problem.elementData.updateFields()
		}

		Button {
			text: qsTr("Integrate")

			onClicked: problem.integrate()
		}

		Button {
			text: qsTr("Write ParaView")

			onClicked: problem.writeParaview()
		}

		Button {
			text: qsTr("View entity")
		}
	}
}
