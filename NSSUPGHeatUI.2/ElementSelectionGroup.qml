import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import ModFEM.NSSUPGHeat 1.0

GroupBox {
	title: qsTr("Element selection")

	Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

	property Problem problem

	ButtonGroup { id: elementSelectionButtonGroup }

	Column {
		Row {
			RadioButton {
				id: elementSelectionAllButton

				checked: true
				text: qsTr("All")

				ButtonGroup.group: elementSelectionButtonGroup
			}

			RadioButton {
				id: elementSelectionSingleButton

				checked: false
				text: qsTr("Single")

				ButtonGroup.group: elementSelectionButtonGroup
			}
		}

		SpinBox {
			enabled: elementSelectionSingleButton.checked

			from: 1
			to: problem.elementData.count.elements

			onEnabledChanged: enabled ? problem.elementData.selectElement(value) : problem.elementData.selectAll()

			onValueChanged: if (enabled) problem.elementData.selectElement(value)
		}
	}
}
