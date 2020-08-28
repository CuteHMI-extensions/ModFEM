import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3

import Qt.labs.platform 1.1
import Qt.labs.settings 1.1

import ModFEM.NSSUPGHeat 1.0

GroupBox {
	title: qsTr("Problem")

	property Problem problem

	RowLayout {
		Label {
			text: qsTr("Problem directory: ")
		}

		TextField {
			id: directoryTextField
			text: problem.directory
			selectByMouse: true

			onEditingFinished: problem.directory = text
		}

		Button {
			text: qsTr("Browse...")
			onClicked: directoryDialog.open()
		}

		FolderDialog {
			id: directoryDialog
			currentFolder: problem.directory

			onAccepted: problem.setDirectoryFromURL(folder)
		}
	}
}
