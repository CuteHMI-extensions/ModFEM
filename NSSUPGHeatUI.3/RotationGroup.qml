import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3
import Qt3D.Core 2.5 as Qt3D

GroupBox {
	title: qsTr("Rotation")

	Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

	property var entity

	GridLayout {
		columns: 2

		Label {
			text: qsTr("XY:")
		}

		Slider {
			id: xySlider

			from: -180
			to: 180

			onValueChanged: entity.transform.rotationZ = value
		}

		Label {
			text: qsTr("XZ:")
		}

		Slider {
			id: xzSlider

			from: -180
			to: 180

			onValueChanged: entity.transform.rotationY = value
		}

		Label {
			text: qsTr("YZ:")
		}

		Slider {
			id: yzSlider


			from: -180
			to: 180

			onValueChanged: entity.transform.rotationX = value
		}
	}

	onEntityChanged: {
		xySlider.value = entity.transform.rotationZ
		xzSlider.value = entity.transform.rotationY
		yzSlider.value = entity.transform.rotationX
	}
}
