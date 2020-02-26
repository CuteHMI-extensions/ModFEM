import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import Qt.labs.platform 1.1
import Qt.labs.settings 1.1

import ModFEM.QtHeat 0.0

Item {
	anchors.fill: parent
	anchors.margins: 10

	Controller {
		id: heat
	}

	Row {
		spacing: 10

		Column {
			Label {
				text: qsTr("Problem id: ") + heat.problemId
			}

			Label {
				text: qsTr("Mesh id: ") + heat.meshId
			}

			Label {
				text: qsTr("Field id: ") + heat.fieldId
			}

			Label {
				text: qsTr("Solution count: ") + heat.solutionCount
			}

			Label {
				text: qsTr("Equation count: ") + heat.equationCount
			}
		}

		ColumnLayout {
			RowLayout {
				Label {
					text: qsTr("Problem directory: ")
				}

				TextField {
					id: problemDirectoryTextField
					text: heat.problemDirectory
					selectByMouse: true

					onEditingFinished: heat.problemDirectory = text
				}

				Button {
					text: qsTr("Browse...")
					onClicked: problemDirectoryDialog.open()
				}

				FolderDialog {
					id: problemDirectoryDialog
					currentFolder: heat.problemDirectory

					onAccepted: heat.setProblemDirectoryFromURL(folder)
				}

			}

			Button {
				text: "Init"

				onClicked: heat.init()
			}
		}
	}


}
