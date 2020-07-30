import QtQml 2.12
import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Input 2.5
import Qt3D.Extras 2.14

import Qt.labs.platform 1.1
import Qt.labs.settings 1.1

import ModFEM.NSSUPGHeat 1.0

Item {
	anchors.fill: parent
	anchors.margins: 10

	Problem {
		id: problem
	}

	Instantiator {
		id: colorMappersInstantiator

		property AbstractColorMapper triangleColorMapper

		property AbstractColorMapper quadTriangleColorMapper

		onObjectAdded: {
			object.elementData = Qt.binding(function() { return problem.elementData })
			triangleColorMapper = object.triangle
			quadTriangleColorMapper = object.quadTriangle
		}
	}

	RowLayout {
		id: settingsLayout

		anchors.fill: parent

		property var transformEntity: elementsEntity

		ColumnLayout {
			Layout.alignment: Qt.AlignTop

			ElementSelectionGroup {
				id: elementSelectionGroup

				problem: problem
			}

			VisibilityGroup {
				id: visibilityGroup
			}

			GridLayout {
				columns: 2

				Label {
					text: qsTr("Color:")
				}

				ComboBox {
					model: [qsTr("BC"), qsTr("Temperature")]

					onCurrentIndexChanged: {
						switch (currentIndex) {
							case 0:
								colorMappersInstantiator.delegate = Qt.createComponent("BCColorMappers.qml")
								break
							case 1:
								colorMappersInstantiator.delegate = Qt.createComponent("TemperatureColorMappers.qml")
						}
					}
				}

				Label {
					text: qsTr("Point size:")
				}

				Slider {
					id: pointSizeSlider

					from: 1
					to: 20
					stepSize: 1
					value: rootEntity.pointSize

					onValueChanged: rootEntity.pointSize = value
				}

				Label {
					text: qsTr("Line width:")
				}

				Slider {
					id: lineWidthSlider

					from: 1
					to: 20
					stepSize: 1
					value: rootEntity.lineWidth

					onValueChanged: rootEntity.lineWidth = value
				}

				Label {
					text: qsTr("Transform:")
				}

				ComboBox {
					model: ["Elements", "Clip plane"]

					property var entities: [elementsEntity, clipPlane0Entity]

					onCurrentIndexChanged: settingsLayout.transformEntity = entities[currentIndex]
				}
			}

			ScaleGroup {
				entity: settingsLayout.transformEntity
			}

			RotationGroup {
				entity: settingsLayout.transformEntity
			}

			TranslationGroup {
				entity: settingsLayout.transformEntity
			}
		}

		ColumnLayout {
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

			Row {
				spacing: 5

				Button {
					text: qsTr("Init")

					onClicked: problem.init()
				}

				Button {
					text: qsTr("Start")

					onClicked: problem.start()
				}

				Button {
					text: qsTr("Stop")

					onClicked: problem.stop()
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
			}

			Rectangle {
				Layout.fillHeight: true
				Layout.fillWidth: true
				color: "black"

				Scene3D {
					id: scene3d

					anchors.fill: parent

					RootEntity {
						id: rootEntity

						property ShaderData clipPlanesData: ShaderData {
							property int count: 1
							property ShaderDataArray planes: ShaderDataArray {
								ShaderData {
									property vector4d equation: clipPlane0Entity.equation
								}
							}
						}

						ElementsEntity {
							id: elementsEntity

							elementData: problem.elementData
							triangleColorMapper: colorMappersInstantiator.triangleColorMapper
							quadTriangleColorMapper: colorMappersInstantiator.quadTriangleColorMapper
							nodesEnabled: visibilityGroup.nodesEnabled
							linesEnabled: visibilityGroup.linesEnabled
							facesEnabled: visibilityGroup.facesEnabled
							clipPlanesData: rootEntity.clipPlanesData
						}

						ClipPlaneEntity {
							id: clipPlane0Entity
						}
					}
				}

				Row {
					x: parent.width - width - 10
					y: 10
					spacing: 10

					/// @todo turn rectangle to icon with arrow keys symbol.
					Rectangle {
						opacity: scene3d.focus

						width: 10
						height: 10

						color: "transparent"
						border.color: "white"
						border.width: 1

						Behavior on opacity { NumberAnimation {} }
					}
				}
			}
		}

		ColumnLayout {
			Layout.alignment: Qt.AlignTop

			ProblemInfo {
				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

				problem: problem
			}

			MeshInfo {
				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

				elementData: problem.elementData
			}

		}
	}

}
