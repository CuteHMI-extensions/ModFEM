import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Scene3D 2.0
import QtQuick.Scene2D 2.9
import Qt3D.Render 2.13
import Qt3D.Core 2.0
import Qt3D.Render 2.14
import Qt3D.Input 2.0
import Qt3D.Extras 2.14

import CuteHMI.Symbols.HVAC 0.0

Item {
	Rectangle {
		id: scene
		anchors.fill: parent
		//		anchors.margins: 50
		color: "black"

		Scene3D {
			id: scene3d
			anchors.fill: parent
			anchors.margins: 10
			focus: true
			aspects: ["input", "logic"]
			cameraAspectRatioMode: Scene3D.AutomaticAspectRatio

			Entity {
				id: sceneRoot

				Camera {
					id: camera
					projectionType: CameraLens.PerspectiveProjection
					position: Qt.vector3d( 0.0, 0.0, 20 )
				}

				FirstPersonCameraController {
					camera: camera
				}

				components: [
					RenderSettings {
						activeFrameGraph: ForwardRenderer {
							camera: camera
							clearColor: "white"
						}
						pickingSettings.pickMethod: PickingSettings.TrianglePicking
					},
					InputSettings {}
				]

				Entity {
					id: logoEntity

					Mesh {
						id: logoMesh
						source: "Qt_logo.obj"
					}

//					PhongAlphaMaterial{
//						id: logoAlphaMaterial

//						ambient: Qt.rgba( 1, 0, 0, 1.0 )
//						diffuse: Qt.rgba( 1, 0, 0, 1.0 )
////						specular: Qt.rgba(1, 0, 0, 1.0 )
//						shininess: 1.0
//						alpha: 0.1
//					}

					PhongMaterial {
						id: logoMaterial
						ambient: Qt.rgba( 0.1, 0.1, 0.1, 1.0 )
						shininess: 1.0
					}

					components: [ logoMesh, logoMaterial ]
				}

				Entity {
					id: planeEntity

					Transform {
						id: planeTransform
						translation: Qt.vector3d( -25, 0, -100)
//						scale3D: Qt.vector3d(1000, 1000, 100)
						rotation: fromAxisAndAngle(Qt.vector3d(1,0,0), -90)
					}

					PlaneMesh {
						id: planeMesh

						width: 10
						height: 10
						meshResolution: Qt.size(20, 20)
					}

					PhongMaterial {
						id: planeMaterialFallback
						diffuse: Qt.rgba(1.0, 0.0, 0.0, 1.0)
						ambient: Qt.rgba(1.0, 0.0, 0.0, 1.0)
						shininess: 1.0
					}

					TextureMaterial {
						id: planeMaterial
//						texture: plane2dTexture
						alphaBlending: true
						texture: TextureLoader {
							id: textureLoader
							source: "man.png"
							mirrored: true
						}
					}

					CustomMaterial {
						id: customMaterial

						texture: pumpTexture

						Texture2D {
							id: manTexture

							generateMipMaps: true
							minificationFilter: Texture.Linear
							magnificationFilter: Texture.Linear
							wrapMode { // #TRYIT: change wrap mode
								x: WrapMode.Repeat
								y: WrapMode.Repeat
							}

							TextureImage {
								mipLevel: 0
								source: Qt.resolvedUrl("man.png")
							}
						}


						Scene2D {
							id: pumpScene2D
							output: RenderTargetOutput {
								attachmentPoint: RenderTargetOutput.Color0
								texture: Texture2D {
									id: pumpTexture
									width: pump2.width
									height: pump2.height
									format: Texture.RGBA8_UNorm
									generateMipMaps: true
									magnificationFilter: Texture.Linear
									minificationFilter: Texture.LinearMipMapLinear
									wrapMode {
										x: WrapMode.ClampToEdge
										y: WrapMode.ClampToEdge
									}
								}
							}

							mouseEnabled: false

							Pump {
								id: pump2

								active: true
							}
						}
					}

					ObjectPicker {
						id: planePicker
						hoverEnabled: true
						dragEnabled: true
					}

					Scene2D {
						id: planeScene2d
						output: RenderTargetOutput {
							attachmentPoint: RenderTargetOutput.Color0
							texture: Texture2D {
								id: plane2dTexture
								width: 100
								height: 80
								format: Texture.RGBA8_UNorm
								generateMipMaps: true
								magnificationFilter: Texture.Linear
								minificationFilter: Texture.LinearMipMapLinear
								wrapMode {
									x: WrapMode.ClampToEdge
									y: WrapMode.ClampToEdge
								}
							}
						}

						entities: [ planeEntity ]
						mouseEnabled: false

						CentrifugalFan {
							active: true
						}

					}

//					components: [ planeTransform, planeMesh, planeMaterial, planePicker ]
					components: [ planeTransform, planeMesh, customMaterial, planePicker ]
				}

				PlaneEntity {
					id: customPlaneEntity
				}


				Entity {
					id: cube

					components: [cubeTransform, cubeMaterial, cubeMesh, cubePicker]

					property real rotationAngle: 0

					Behavior on rotationAngle {
						RotationAnimation {
							direction: RotationAnimation.Shortest
							duration: 450
						}
					}

					RotationAnimation on rotationAngle {
						loops: Animation.Infinite
						from: 0; to: 360
						duration: 4000
						onStopped: cube.rotationAngle = 0
					}

					Transform {
						id: cubeTransform
						translation: Qt.vector3d(2, 0, 10)
						scale3D: Qt.vector3d(1, 1, 1)
						rotation: fromAxisAndAngle(Qt.vector3d(0,1,0), cube.rotationAngle)
					}

					CuboidMesh {
						id: cubeMesh
					}

					ObjectPicker {
						id: cubePicker
						hoverEnabled: true
						dragEnabled: true

						// Explicitly require a middle click to have the Scene2D grab the mouse
						// events from the picker
						onPressed: {
							if (pick.button === PickEvent.MiddleButton) {
								qmlTexture.mouseEnabled = !qmlTexture.mouseEnabled
								logoControls.enabled = !logoControls.enabled
							}
						}
					}

					TextureMaterial {
						id: cubeMaterial
						texture: offscreenTexture
					}

					Scene2D {
						id: qmlTexture
						output: RenderTargetOutput {
							attachmentPoint: RenderTargetOutput.Color0
							texture: Texture2D {
								id: offscreenTexture
//								width: 256
//								height: 1024
								width: pump.width
								height: pump.height
								format: Texture.RGBA8_UNorm
								generateMipMaps: true
								magnificationFilter: Texture.Linear
								minificationFilter: Texture.LinearMipMapLinear
								wrapMode {
									x: WrapMode.ClampToEdge
									y: WrapMode.ClampToEdge
								}
							}
						}

						entities: [ cube ]
						mouseEnabled: false

//						Rectangle {
//							width: offscreenTexture.width
//							height: offscreenTexture.height
//							color: "blue"
//						}

						Pump {
							id: pump

							active: true
						}

//						LogoControls {
//							id: logoControls
//							width: offscreenTexture.width
//							height: offscreenTexture.height
//						}
					}
				}
			}
		}
	}
}
