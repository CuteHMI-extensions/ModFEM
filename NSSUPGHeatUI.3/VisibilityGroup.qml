import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3

GroupBox {
	title: qsTr("Visibility")

	Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

	property alias nodesEnabled: nodesCheckbox.checked

	property alias linesEnabled: linesCheckbox.checked

	property alias facesEnabled: facesCheckbox.checked

	property alias capsEnabled: capsCheckbox.checked

	property alias alpha: alphaSlider.value

	Column {
		CheckBox {
			id: nodesCheckbox

			checked: false
			text: qsTr("Nodes")
		}

		CheckBox {
			id: linesCheckbox

			checked: true
			text: qsTr("Lines")
		}

		CheckBox {
			id: facesCheckbox

			checked: true
			text: qsTr("Faces")
		}

		CheckBox {
			id: capsCheckbox

			checked: false
			text: qsTr("Caps")
		}

		RowLayout {
			Label {
				text: qsTr("Alpha:")
			}

			Slider {
				id: alphaSlider

				from: 0.0
				to: 1.0

				value: 1.0
			}
		}
	}
}
