import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3

GroupBox {
	title: activeRange.name + " " + qsTr("color range")

	property QtObject temperature: QtObject {
		property string name: qsTr("Temperature")
		property real minimum: 18
		property real maximum: 40.0
	}

	property QtObject velocityMagnitude: QtObject {
		property string name: qsTr("Velocity magnitude")
		property real minimum: 0.0
		property real maximum: 5.0
	}

	property QtObject pressure: QtObject {
		property string name: qsTr("Pressure")
		property real minimum: 0
		property real maximum: 30.0
	}

	property QtObject activeRange: temperature

	onActiveRangeChanged: {
		minimumField.text = activeRange.minimum
		maximumField.text = activeRange.maximum
	}

	RowLayout {
		RowLayout {
			Label {
				text: qsTr("Min:")
			}

			TextField {
				id: minimumField

				validator: DoubleValidator {bottom: -1000000; top: 1000000}
				text: activeRange.minimum
				onAccepted: activeRange.minimum = text
			}

			Label {
				text: qsTr("Max:")
			}

			TextField {
				id: maximumField

				validator: DoubleValidator {bottom: -1000000; top: 1000000}
				text: activeRange.maximum
				onAccepted: activeRange.maximum = text
			}
		}
	}
}
