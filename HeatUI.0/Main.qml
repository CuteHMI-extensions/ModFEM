import QtQuick 2.0
import QtQuick.Controls 2.5
import QtQuick.Layouts 1.3

import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Input 2.5
import Qt3D.Extras 2.14

import Qt.labs.platform 1.1
import Qt.labs.settings 1.1

import ModFEM.Heat 0.0

Item {
	anchors.fill: parent
	anchors.margins: 10

	Problem {
		id: heatProblem
	}

	RowLayout {
		anchors.fill: parent

		ColumnLayout {
			Layout.alignment: Qt.AlignTop

			GroupBox {
				title: qsTr("Element selection")

				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

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
						to: heatProblem.elementData.count.elements

					}
				}
			}

			GroupBox {
				title: qsTr("Visibility")

				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

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
				}
			}

			//			GroupBox {
			//				title: qsTr("Surface visibility")

			//				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

			//				Column {
			//					CheckBox {
			//						id: surfaceCheckbox

			//						checked: false
			//						text: qsTr("Nodes")
			//					}

			//					CheckBox {
			//						id: surfaceLinesCheckbox

			//						checked: false
			//						text: qsTr("Lines")
			//					}

			//					CheckBox {
			//						id: surfaceFacesCheckbox

			//						checked: true
			//						text: qsTr("Faces")
			//					}
			//				}
			//			}


			GridLayout {
				columns: 2

				Label {
					text: qsTr("Point size:")
				}

				Slider {
					from: 1
					to: 20
					stepSize: 1
					value: renderSettings.pointSize

					onValueChanged: renderSettings.pointSize = value
				}

				Label {
					text: qsTr("Line width:")
				}

				Slider {
					from: 1
					to: 20
					stepSize: 1
					value: renderSettings.lineWidth

					onValueChanged: renderSettings.lineWidth = value
				}
			}

		}

		ColumnLayout {
			RowLayout {
				Label {
					text: qsTr("Problem directory: ")
				}

				TextField {
					id: directoryTextField
					text: heatProblem.directory
					selectByMouse: true

					onEditingFinished: heatProblem.directory = text
				}

				Button {
					text: qsTr("Browse...")
					onClicked: directoryDialog.open()
				}

				FolderDialog {
					id: directoryDialog
					currentFolder: heatProblem.directory

					onAccepted: heatProblem.setDirectoryFromURL(folder)
				}
			}

			Row {
				spacing: 5

				Button {
					text: qsTr("Init")

					onClicked: heatProblem.init()
				}

				Button {
					text: qsTr("Solve")

					onClicked: heatProblem.solve()
				}

				Button {
					text: qsTr("Integrate")

					onClicked: heatProblem.integrate()
				}

				Button {
					text: qsTr("Write ParaView")

					onClicked: heatProblem.writeParaview()
				}
			}

			Rectangle {
				Layout.fillHeight: true
				Layout.fillWidth: true
				color: "black"

				Scene3D {
					id: scene3d

					anchors.fill: parent

					aspects: ["input", "logic"]

					cameraAspectRatioMode: Scene3D.AutomaticAspectRatio

					Entity {
						id: sceneRoot

						components: [renderSettings, inputSettings]

						RenderSettings {
							id: renderSettings

							property real pointSize: 10

							property real lineWidth: 3

							//							ForwardRenderer {
							//								clearColor: "transparent"
							//								camera: camera1
							//							}

							RenderSurfaceSelector {
								ClearBuffers {
									buffers : ClearBuffers.ColorDepthBuffer

									CameraSelector {
										camera: camera1

										RenderStateSet {
											renderStates: [
												DepthTest { depthFunction: DepthTest.Less },
												CullFace { mode: CullFace.NoCulling },
												LineWidth { value: renderSettings.lineWidth },
												PointSize { value: renderSettings.pointSize; sizeMode: PointSize.Fixed }
											]
										}
									}
								}
							}
						}

						InputSettings {
							id: inputSettings
						}

						MouseDevice {
							id: mouseDevice
						}

						MouseHandler {
							sourceDevice: mouseDevice

							onPressed: scene3d.focus = true
						}

						Camera {
							id: camera1
							projectionType: CameraLens.PerspectiveProjection
							fieldOfView: 45
							nearPlane : 0.1
							farPlane : 1000.0
							position: Qt.vector3d( 0.0, 0.0, 40.0 )
							//							position: Qt.vector3d( 0.0, 0.0, 1.0 )
							upVector: Qt.vector3d( 0.0, 1.0, 0.0 )
							viewCenter: Qt.vector3d( 0.0, 0.0, 0.0 )

							//							onViewVectorChanged: console.log(viewVector)
						}

						FirstPersonCameraController {
							camera: camera1
							linearSpeed: 100
						}

						//						TestEntity {
						//							nodesEnabled: nodesCheckbox.checked

						//							rotationX: rotationXSlider.value
						//							rotationY: rotationYSlider.value
						//							rotationZ: rotationZSlider.value
						//							scale: scaleSlider.value
						//						}

						ElementsEntity {
							elementData: heatProblem.elementData

							rotationX: rotationXSlider.value
							rotationY: rotationYSlider.value
							rotationZ: rotationZSlider.value
							scale: Math.pow(10, scaleSlider.value)

							nodesEnabled: nodesCheckbox.checked
							linesEnabled: linesCheckbox.checked
							facesEnabled: facesCheckbox.checked
						}

						//						SurfaceEntity {
						//							surface: heatProblem.faces

						//							rotationX: rotationXSlider.value
						//							rotationY: rotationYSlider.value
						//							rotationZ: rotationZSlider.value
						//							scale: Math.pow(10, scaleSlider.value)

						//							nodesEnabled: surfaceNodesCheckbox.checked
						//							linesEnabled: surfaceLinesCheckbox.checked
						//							facesEnabled: surfaceFacesCheckbox.checked
						//						}
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

				Column {
					Slider {
						id: rotationXSlider

						from: -180
						to: 180
					}

					Slider {
						id: rotationYSlider

						from: -180
						to: 180
					}

					Slider {
						id: rotationZSlider

						from: -180
						to: 180
					}

					Slider {
						id: scaleSlider
						//						from: 0.1
						//						to: 2
						from: -10
						to: 10
						value: 0
					}
				}
			}
		}

		ColumnLayout {
			Layout.alignment: Qt.AlignTop

			ProblemInfo {
				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

				problem: heatProblem
			}

			MeshInfo {
				Layout.minimumWidth: Math.max(parent.width, Layout.preferredWidth)

				elementData: heatProblem.elementData
			}

		}
	}

}
