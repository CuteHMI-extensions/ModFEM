import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3
import Qt3D.Core 2.5 as Qt3D

GroupBox {
	title: qsTr("Translation")

	Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

	property var entity

	GridLayout {
		columns: 2

		Label {
			text: qsTr("X:")
		}

		Slider {
			id: xSlider

			from: -5
			to: 5

			onValueChanged: entity.transform.translation.x = value
		}

		Label {
			text: qsTr("Y:")
		}

		Slider {
			id: ySlider

			from: -5
			to: 5

			onValueChanged: entity.transform.translation.y = value
		}

		Label {
			text: qsTr("Z:")
		}

		Slider {
			id: zSlider

			from: -5
			to: 5

			onValueChanged: entity.transform.translation.z = value
		}
	}

	onEntityChanged: {
		xSlider.value = entity.transform.translation.x
		ySlider.value = entity.transform.translation.y
		zSlider.value = entity.transform.translation.z
	}
}
