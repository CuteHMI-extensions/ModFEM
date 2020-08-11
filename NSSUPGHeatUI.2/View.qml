import QtQml 2.12
import QtQuick 2.14
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import Qt3D.Core 2.14
import Qt3D.Render 2.14
import Qt3D.Input 2.5
import Qt3D.Extras 2.14

import Qt.labs.platform 1.1
import Qt.labs.settings 1.1

import CuteHMI.DataAcquisition 0.0
import CuteHMI.Services 2.0

import ModFEM.NSSUPGHeat 1.0

import QtMultimedia 5.12	//temp?

Item {
	anchors.fill: parent
	anchors.margins: 10

	Problem {
		id: problem

		elementData.probes: [
			ScalarProbe {
				id: temperatureProbe1

				position: temperatureProbe1Entity.transform.translation
				field: "temperature"
			},

			Vector3Probe {
				id: velocityProbe1

				position: velocityProbe1Entity.transform.translation
				field: "velocity"
			}
		]
	}

	Project {
		id: project
	}

	Service {
		id: historyService

		name: "ModFEM History Service"

		HistoryWriter {
			schema: project.schema

			TagValue {
				name: "temperatureProbe1"
				value: temperatureProbe1.value
			}
		}
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

			// Temporary label to show probe value
			Label {
				text: "Value: " + temperatureProbe1.value
			}

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
					model: [qsTr("BC"), qsTr("Temperature"), qsTr("Pressure"), qsTr("Velocity magnitude"), qsTr("Velocity X"), qsTr("Velocity Y"), qsTr("Velocity Z")]

					onCurrentIndexChanged: {
						switch (currentIndex) {
						case 0:
							colorMappersInstantiator.delegate = Qt.createComponent("BCColorMappers.qml")
							break
						case 1:
							colorMappersInstantiator.delegate = Qt.createComponent("TemperatureColorMappers.qml")
							break
						case 2:
							colorMappersInstantiator.delegate = Qt.createComponent("PressureColorMappers.qml")
							break
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
					model: ["Elements", "Clip plane", "Temperature probe 1", "Velocity probe 1"]

					property var entities: [elementsEntity, clipPlane0Entity, temperatureProbe1Entity, velocityProbe1Entity]

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

					onClicked: {
						problem.start()
						ServiceManager.start()
					}
				}

				Button {
					text: qsTr("Stop")

					onClicked: {
						problem.stop()
						ServiceManager.stop()
					}
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

				Button {
					text: qsTr("View entity")
				}
			}

			Rectangle {
				Layout.fillHeight: true
				Layout.fillWidth: true
				color: "black"

//				VideoOutput {
//					anchors.fill: parent

//					source: camera

//					Camera {
//						id: camera
//						// You can adjust various settings in here

//						deviceId: QtMultimedia.availableCameras[0].deviceId
//					}

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

							ClipPlaneEntity {
								id: clipPlane0Entity
							}

							ElementsEntity {
								id: elementsEntity

								elementData: problem.elementData
								triangleColorMapper: colorMappersInstantiator.triangleColorMapper
								quadTriangleColorMapper: colorMappersInstantiator.quadTriangleColorMapper
								nodesEnabled: visibilityGroup.nodesEnabled
								linesEnabled: visibilityGroup.linesEnabled
								facesEnabled: visibilityGroup.facesEnabled
								alpha: visibilityGroup.alpha
								clipPlanesData: rootEntity.clipPlanesData

								ProbeEntity {
									id: temperatureProbe1Entity

									transform.translation.x: -5
								}

								ProbeEntity {
									id: velocityProbe1Entity

									transform.translation.x: 5
								}

								ArrowEntity {
									id: velocity1ArrowEntity

									transform.translation: velocityProbe1Entity.transform.translation
									vector: velocityProbe1.value
								}
							}

							CentrifugalFanEntity {

							}

							NumberDisplayEntity {
								display.value: temperatureProbe1.value

								transform.translation: Qt.vector3d(temperatureProbe1Entity.transform.worldMatrix.column(3).x, elementsEntity.maxExtent.y + 2.5, temperatureProbe1Entity.transform.worldMatrix.column(3).z)
							}
							NumberDisplayEntity {

								display.value: velocityProbe1.value.length()
								display.unit: "m/s"

								transform.translation: Qt.vector3d(velocityProbe1Entity.transform.worldMatrix.column(3).x, elementsEntity.maxExtent.y + 2.5, velocityProbe1Entity.transform.worldMatrix.column(3).z)
							}
						}
//					}
				}

				Row {
					id: header3d

					x: parent.width - width - 10
					y: 10
					spacing: 10

					Image {
						opacity: scene3d.focus

						fillMode: Image.PreserveAspectFit
						sourceSize.width: 25
						sourceSize.height: 25
						source: "images/keyboard-white-18dp.svg"

						Behavior on opacity { NumberAnimation {} }
					}
				}

				Row {
					id: footer3d

					x: 10
					y: parent.height - height
					spacing: 10

					Axes {
						width: 150
						height: 150

						rotationMatrix: rootEntity.cameraRotationMatrix
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

			SolutionInfo {
				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

				problem: problem
			}
		}
	}

}
