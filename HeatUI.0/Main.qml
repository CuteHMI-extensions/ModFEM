import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import Qt.labs.platform 1.1
import Qt.labs.settings 1.1

import ModFEM.QtHeat 0.0

Item {
	anchors.fill: parent
	anchors.margins: 10

	Problem {
		id: heatProblem
	}

	RowLayout {
		anchors.fill: parent

		Column {
			Layout.alignment: Qt.AlignTop

			Label {
				text: qsTr("Problem id: ") + heatProblem.problemId
			}

			Label {
				text: qsTr("Mesh id: ") + heatProblem.meshId
			}

			Label {
				text: qsTr("Field id: ") + heatProblem.fieldId
			}

			Label {
				text: qsTr("Solution count: ") + heatProblem.solutionCount
			}

			Label {
				text: qsTr("Equation count: ") + heatProblem.equationCount
			}
		}

		ColumnLayout {
			RowLayout {
				Label {
					text: qsTr("Problem directory: ")
				}

				TextField {
					id: directoryTextField
					text: heatProblem.directory
					selectByMouse: true

					onEditingFinished: heatProblem.directory = text
				}

				Button {
					text: qsTr("Browse...")
					onClicked: directoryDialog.open()
				}

				FolderDialog {
					id: directoryDialog
					currentFolder: heatProblem.directory

					onAccepted: heatProblem.setDirectoryFromURL(folder)
				}
			}

			Row {
				spacing: 5

				Button {
					text: "Init"

					onClicked: heatProblem.init()
				}

				Button {
					text: "Reset buffer"

					onClicked: heatProblem.resetBuffer()
				}
			}

			Rectangle {
				Layout.fillHeight: true
				Layout.fillWidth: true
				color: "black"

				Scene3D {
					id: scene3d

					anchors.fill: parent

//					controller: heat
				}
			}
		}
	}

}
