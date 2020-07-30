import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Scene3D 2.0
import QtQuick.Scene2D 2.9
import Qt3D.Render 2.13
import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Input 2.0
import Qt3D.Extras 2.0

import CuteHMI.Symbols.HVAC 1.0

Item {
	Rectangle {
		id: scene
		anchors.fill: parent
		anchors.margins: 50
		color: "black"

		Pump {
			x: scene3d.width * 0.5 + xTranslateSlider.value
			y: scene3d.height * 0.5 - yTranslateSlider.value
			z: 1
			active: true
			scale: scaleSlider.value

		}

		Scene3D {
			id: scene3d
			anchors.fill: parent
			anchors.margins: 10
			focus: true
			aspects: ["input", "logic"]
			cameraAspectRatioMode: Scene3D.AutomaticAspectRatio

			Entity {
				id: sceneRoot

				components: [
					RenderSettings {
						activeFrameGraph: ForwardRenderer {
							camera: camera
							clearColor: "transparent"
						}
					},
					InputSettings {}
				]

				Camera {
					id: camera
					projectionType: CameraLens.PerspectiveProjection
					fieldOfView: 90
					nearPlane : 0.1
					farPlane : 1000.0
					position: Qt.vector3d( 0.0, 0.0, 200.0 )
					upVector: Qt.vector3d( 0.0, 1.0, 0.0 )
					viewCenter: Qt.vector3d( 0.0, 0.0, 0.0 )

					onViewCenterChanged: {
						console.log(viewCenter, viewVector, position, viewMatrix)
					}
				}

				FirstPersonCameraController {
					camera: camera
				}

//				Entity {
//					id: test

//					GeometryRenderer {
//						instanceCount: 1
//						indexOffset: 0
//						firstInstance: 0
//						primitiveType: GeometryRenderer.Triangles
//						geometry: Geometry { ... }
//					}
//				}

				Entity {
					id: planeEntity

					Transform {
						id: planeTransform
						translation: Qt.vector3d( -90, 0, 110)
//						scale3D: Qt.vector3d(1000, 1000, 100)
						rotation: fromAxisAndAngle(Qt.vector3d(1,0,0), 90)
					}

					PlaneMesh {
						id: planeMesh

						width: 100
						height: 100
						meshResolution: Qt.size(20, 20)
					}

					PhongMaterial {
						id: planeMaterial
						diffuse: Qt.rgba(0.0, 1.0, 0.0, 1.0)
						ambient: Qt.rgba(0.0, 1.0, 0.0, 1.0)
						shininess: 1.0
					}

					components: [ planeTransform, planeMesh, planeMaterial ]
				}


				Entity {
					id: cubeEntity

					Transform {
						id: cubeTransform

						scale: scaleSlider.value
						rotationX: xRotationSlider.value
						rotationY: yRotationSlider.value
						rotationZ: zRotationSlider.value
						translation: Qt.vector3d(xTranslateSlider.value, yTranslateSlider.value, 0.0)
					}

//					Mesh {
//						id: cubeMesh
//						source: Qt.resolvedUrl("table.obj")
//					}

					PhongMaterial {
						id: cubeMaterial
						ambient: "green"
					}

					SceneLoader {
						id: cubeScene

//						source: Qt.resolvedUrl("cube.fbx")
						source: Qt.resolvedUrl("cube_ascii.fbx")
//						source: Qt.resolvedUrl("handgun_binary.fbx")
//						source: Qt.resolvedUrl("table.fbx")
//						source: Qt.resolvedUrl("cube_r14.dxf")
//						source: Qt.resolvedUrl("cube_2018_ascii.dxf")
//						source: Qt.resolvedUrl("cube_2018_binary.dxf")
					}

					Mesh {
						id: cubeMesh

						source: Qt.resolvedUrl("cube.obj")
					}

//					components: [cubeTransform, cubeScene, cubeMaterial]
					components: [cubeTransform, cubeScene, cubeMaterial]
//					components: [cubeTransform, cubeMesh]
				}
			}
		}

		Column {
			Slider {
				id: xRotationSlider

				from: -180
				to: 180
			}

			Slider {
				id: yRotationSlider

				from: -180
				to: 180
			}

			Slider {
				id: zRotationSlider

				from: -180
				to: 180
			}

			Slider {
				id: scaleSlider
				from: 0.1
				to: 2
				value: 1
			}

			Slider {
				id: xTranslateSlider
				from: -100
				to: 100
			}

			Slider {
				id: yTranslateSlider
				from: -100
				to: 100
			}

			Button {
				text: "dump camera projection matrix"
				onClicked: {
					console.log(camera.viewMatrix)
					console.log(cubeTransform.matrix)
					console.log(camera.projectionMatrix.times(cubeTransform.matrix))
					var mat = camera.projectionMatrix.times(camera.viewMatrix).times(cubeTransform.matrix)
					var vec = Qt.vector4d(xTranslateSlider.value, yTranslateSlider.value, 100.0, 0.0)
					console.log(vec.times(mat), mat.times(vec))
				}
			}
		}
	}
}
