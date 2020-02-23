import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Scene3D 2.0
import Qt3D.Render 2.13
import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Input 2.0
import Qt3D.Extras 2.0

Item {
	Text {
		text: "Click me!"
		anchors.top: parent.top
		anchors.topMargin: 10
		anchors.horizontalCenter: parent.horizontalCenter

		MouseArea {
			anchors.fill: parent
			onClicked: animation.start()
		}
	}

	Text {
		text: "Multisample: " + scene3d.multisample
		anchors.bottom: parent.bottom
		anchors.bottomMargin: 10
		anchors.horizontalCenter: parent.horizontalCenter

		MouseArea {
			anchors.fill: parent
			onClicked: scene3d.multisample = !scene3d.multisample
		}
	}

	Rectangle {
		id: scene
		anchors.fill: parent
		anchors.margins: 50
		color: "black"

		transform: Rotation {
			id: sceneRotation
			axis.x: 1
			axis.y: 0
			axis.z: 0
			origin.x: scene.width / 2
			origin.y: scene.height / 2
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
					fieldOfView: 45
					nearPlane : 0.1
					farPlane : 1000.0
					position: Qt.vector3d( 0.0, 0.0, 40.0 )
					upVector: Qt.vector3d( 0.0, 1.0, 0.0 )
					viewCenter: Qt.vector3d( 0.0, 0.0, 0.0 )
				}

				FirstPersonCameraController {
					camera: camera
				}

				Entity {
					id: cubeEntity

					Transform {
						id: cubeTransform

						scale: scaleSlider.value
						rotationX: xRotationSlider.value
						rotationY: yRotationSlider.value
						rotationZ: zRotationSlider.value
					}

					PhongMaterial {
						id: cubeMaterial
						ambient: "green"
					}

//					Mesh {
//						id: cubeMesh
//						source: Qt.resolvedUrl("table.obj")
//					}

					SceneLoader {
						id: cubeScene

//						source: Qt.resolvedUrl("cube.fbx")
						source: Qt.resolvedUrl("cube_ascii.fbx")
//						source: Qt.resolvedUrl("handgun_binary.fbx")
//						source: Qt.resolvedUrl("table.fbx")
//						source: Qt.resolvedUrl("cube_r14.dxf")
//						source: Qt.resolvedUrl("cube_2018_ascii.dxf")
//						source: Qt.resolvedUrl("cube_2018_binary.dxf")

						Component.onCompleted: {
//							console.log("cubeScene")
//							console.log("entities: " + cubeScene.entities)
//							console.log("childNodes: " + childNodes.length)
//							console.log("childNodes: " + childNodes[0])
//							console.log("childNodes: " + childNodes[1])
//							console.log("entity names: " + cubeScene.entityNames())
						}
					}

					components: [cubeTransform, cubeScene, cubeMaterial]

					Component.onCompleted: {
						console.log("cubeEntity")
						console.log("components: " + components.length)
						console.log("childNodes: " + childNodes.length)
						for (var i = 0; i < childNodes.length; i++)
							console.log("childNode[" + i + "]: " + childNodes[i])

						console.log("components[0].length: " + components[0])
						console.log("components[1].length: " + components[1])
						console.log("components[2].length: " + components[2])
////						console.log("components[1]: " + components[1].components)
//						console.log("components[1].childNodes: " + components[1].childNodes.length)
//						console.log("components[1].childNodes: " + components[1].childNodes[0])
//						console.log("components[1].childNodes: " + components[1].childNodes[1])
//						console.log("components[1].entityNames: " + components[1].entityNames().length)
//						console.log("childNodes: " + childNodes.length)
//						console.log("childNodes: " + childNodes[0])
//						console.log("childNodes: " + childNodes[1])
//						console.log("childNodes: " + childNodes[2])
//						console.log("childNodes: " + childNodes[3])
//						console.log("childNodes: " + childNodes[4])
//						console.log(cubeScene.entityNames())
//						console.log(cubeScene.entity(""))
					}

				}

				Component.onCompleted: {
					console.log("sceneRoot")
					console.log("components: " + components.length)
					console.log("childNodes: " + childNodes.length)
					for (var i = 0; i < childNodes.length; i++)
						console.log("childNode[" + i + "]: " + childNodes[i])
//					console.log("components[0].length: " + components[0])
//					console.log("components[1].length: " + components[1])
//					console.log("components[2].length: " + components[2])
////						console.log("components[1]: " + components[1].components)
//						console.log("components[1].childNodes: " + components[1].childNodes.length)
//						console.log("components[1].childNodes: " + components[1].childNodes[0])
//						console.log("components[1].childNodes: " + components[1].childNodes[1])
//						console.log("components[1].entityNames: " + components[1].entityNames().length)
//						console.log("childNodes: " + childNodes.length)
//						console.log("childNodes: " + childNodes[0])
//						console.log("childNodes: " + childNodes[1])
//						console.log("childNodes: " + childNodes[2])
//						console.log("childNodes: " + childNodes[3])
//						console.log("childNodes: " + childNodes[4])
//						console.log(cubeScene.entityNames())
//						console.log(cubeScene.entity(""))
				}

			}

//			SceneLoader {

//			}

//			AnimatedEntity {}
		}

		Column {
			Button {
				text: "test"
				onClicked: {
//					console.log("components: " + cubeScene.components.length)
//					console.log("childNodes: " + cubeScene.childNodes.length)
					console.log("cubeEntity")
					for (var i = 0; i < cubeEntity.childNodes.length; i++)
						console.log("childNode[" + i + "]: " + cubeEntity.childNodes[i])
					console.log("sceneRoot")
					for (i = 0; i < sceneRoot.childNodes.length; i++)
						console.log("childNode[" + i + "]: " + sceneRoot.childNodes[i])
					console.log(cubeScene.entityNames())

					console.log("cubeEntity.chilNodes[6]")
					var entity = cubeEntity.childNodes[6]
					for (i = 0; i < entity.childNodes.length; i++)
						console.log("entity[" + i + "]: " + entity.childNodes[i])

					console.log("cubeEntity.chilNodes[6][0]")
					for (i = 0; i < entity.childNodes[0].childNodes.length; i++)
						console.log("entity childNodes[" + i + "]: " + entity.childNodes[0].childNodes[i])
					console.log("cubeEntity.chilNodes[6][0]")
//					entity.childNodes[0].childNodes[0] = null
//					entity.childNodes[0].removeComponent(entity.childNodes[0].components[0])
					entity.childNodes[0].components = [cubeMaterial]
					for (i = 0; i < entity.childNodes[0].components.length; i++)
						console.log("entity components[" + i + "]: " + entity.childNodes[0].components[i])
					entity.childNodes[0].components.push(cubeMaterial)
//					entity.childNodes[0].components[0] = cubeMaterial

//					console.log(entity.childNodes[0].childNodes.length)
//					console.log(entity.childNodes[0].components.length)
//					console.log(entity[0].components.length)
//					entity = cubeEntity.childNodes[6].childNodes[0].components
//					for (i = 0; i < entity.childNodes.length; i++)
//						console.log("entity[" + i + "]: " + entity.childNodes[i])
				}
			}

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
				to: 1
				value: 0.1
			}
		}
	}

	SequentialAnimation {
		id: animation

		RotationAnimation {
			to: 45
			duration: 1000
			target: sceneRotation
			property: "angle"
			easing.type: Easing.InOutQuad
		}
		PauseAnimation { duration: 500 }
		NumberAnimation {
			to: 0.5
			duration: 1000
			target: scene
			property: "scale"
			easing.type: Easing.OutElastic
		}
		PauseAnimation { duration: 500 }
		NumberAnimation {
			to: 1.0
			duration: 1000
			target: scene
			property: "scale"
			easing.type: Easing.OutElastic
		}
		PauseAnimation { duration: 500 }
		RotationAnimation {
			to: 0
			duration: 1000
			target: sceneRotation
			property: "angle"
			easing.type: Easing.InOutQuad
		}
	}
}
