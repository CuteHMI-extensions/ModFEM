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

//	Service {
//		id: historyService

//		name: "ModFEM History Service"

//		HistoryWriter {
//			schema: project.schema

//			TagValue {
//				name: "temperatureProbe1"
//				value: temperatureProbe1.value
//			}
//		}
//	}

	Instantiator {
		id: colorMappersInstantiator

		property AbstractColorMapper triangleColorMapper

		property AbstractColorMapper quadTriangleColorMapper

		property AbstractColorMapper capColorMapper

		onObjectAdded: {
			object.elementData = Qt.binding(function() { return problem.elementData })
			if (object.minimum !== undefined)
				object.minimum = Qt.binding(function() { return colorRangesGroup.activeRange.minimum })
			if (object.maximum !== undefined)
				object.maximum = Qt.binding(function() { return colorRangesGroup.activeRange.maximum })

			triangleColorMapper = object.triangle
			quadTriangleColorMapper = object.quadTriangle
			capColorMapper = object.cap
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

				onCapsEnabledChanged: if (problem.meshId && capsEnabled) problem.elementData.clip(clipPlane0Entity.equation)
			}

			GridLayout {
				columns: 2

				Label {
					text: qsTr("Color:")
				}

				ComboBox {
					model: [qsTr("BC"), qsTr("Temperature"), qsTr("Pressure"), qsTr("Velocity magnitude")]

					onCurrentIndexChanged: {
						switch (currentIndex) {
						case 0:
							colorMappersInstantiator.delegate = Qt.createComponent("BCColorMappers.qml")
							break
						case 1:
							colorMappersInstantiator.delegate = Qt.createComponent("TemperatureColorMappers.qml")
							colorRangesGroup.activeRange = colorRangesGroup.temperature
							break
						case 2:
							colorMappersInstantiator.delegate = Qt.createComponent("PressureColorMappers.qml")
							colorRangesGroup.activeRange = colorRangesGroup.pressure
							break
						case 3:
							colorMappersInstantiator.delegate = Qt.createComponent("VelocityColorMappers.qml")
							colorRangesGroup.activeRange = colorRangesGroup.velocityMagnitude
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
					model: ["Elements", "Clip plane", "Temperature probe 1", "Velocity probe 1", "Fan"]

					property var entities: [elementsEntity, clipPlane0Entity, temperatureProbe1Entity, velocityProbe1Entity, fanModelEntity]

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
			Flow {
				Layout.fillWidth: true

				spacing: 5

				SimulationGroup {
					problem: problem
				}

				ProblemGroup {
					problem: problem
				}

				ControlsGroup {
					id: controlsGroup

					problem: problem
				}

				ColorRangesGroup {
					id: colorRangesGroup
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

								onEquationChanged: if (problem.meshId && visibilityGroup.capsEnabled) problem.elementData.clip(equation)
							}

							ElementsEntity {
								id: elementsEntity

								elementData: problem.elementData
								triangleColorMapper: colorMappersInstantiator.triangleColorMapper
								quadTriangleColorMapper: colorMappersInstantiator.quadTriangleColorMapper
								capColorMapper: colorMappersInstantiator.capColorMapper
								nodesEnabled: visibilityGroup.nodesEnabled
								linesEnabled: visibilityGroup.linesEnabled
								facesEnabled: visibilityGroup.facesEnabled
								capsEnabled:  visibilityGroup.capsEnabled
								alpha: visibilityGroup.alpha
								clipPlanesData: rootEntity.clipPlanesData

								ProbeEntity {
									id: temperatureProbe1Entity

									transform.translation.x: 1
								}

								ProbeEntity {
									id: velocityProbe1Entity

									transform.translation.x: -1
								}

								ArrowEntity {
									id: velocity1ArrowEntity

									transform.translation: velocityProbe1Entity.transform.translation
									vector: velocityProbe1.value
								}
							}

							DuctModelEntity {
								id: ductModelEntity

								clipPlanesData: rootEntity.clipPlanesData
							}

							FanModelEntity {
								id: fanModelEntity


								running: controlsGroup.started
							}

							CentrifugalFanEntity {
								fan.active: controlsGroup.started

								transform.translation: Qt.vector3d(fanModelEntity.transform.worldMatrix.column(3).x, elementsEntity.maxExtent.y + 1.0, elementsEntity.transform.worldMatrix.column(3).z)
							}

							NumberDisplayEntity {
								display.value: temperatureProbe1.value

								transform.translation: Qt.vector3d(temperatureProbe1Entity.transform.worldMatrix.column(3).x, elementsEntity.maxExtent.y + 1.0, temperatureProbe1Entity.transform.worldMatrix.column(3).z)
							}
							NumberDisplayEntity {

								display.value: velocityProbe1.value.length()
								display.unit: "m/s"

								transform.translation: Qt.vector3d(velocityProbe1Entity.transform.worldMatrix.column(3).x, elementsEntity.maxExtent.y + 1.0, velocityProbe1Entity.transform.worldMatrix.column(3).z)
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

			SimulationInfo {
				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

				problem: problem
			}
		}
	}

}
