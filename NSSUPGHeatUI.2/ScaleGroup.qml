import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3
import Qt3D.Core 2.5 as Qt3D

GroupBox {
	title: qsTr("Scale")

	Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

	property var entity

	GridLayout {
		columns: 2

		Label {
			text: qsTr("XYZ:")
		}

		Slider {
			id: scaleSlider

			from: -2
			to: 2

			onValueChanged: entity.transform.scale = Math.pow(10, scaleSlider.value)
		}
	}

	onEntityChanged: scaleSlider.value = Math.log10(entity.transform.scale)
}
