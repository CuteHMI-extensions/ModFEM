import QtQuick 2.2 as QQ2
import Qt3D.Core 2.0
import Qt3D.Render 2.14
import Qt3D.Input 2.0
import Qt3D.Extras 2.0
import QtQuick.Scene3D 2.0

Scene3D {
	anchors.fill: parent

	aspects: ["input", "logic"]

	cameraAspectRatioMode: Scene3D.AutomaticAspectRatio

	Entity {
		id: sceneRoot

		Camera {
			id: camera
			projectionType: CameraLens.PerspectiveProjection
			fieldOfView: 45
			aspectRatio: 16/9
			nearPlane : 0.1
			farPlane : 1000.0
			position: Qt.vector3d( 0.0, 0.0, -40.0 )
			upVector: Qt.vector3d( 0.0, 1.0, 0.0 )
			viewCenter: Qt.vector3d( 0.0, 0.0, 0.0 )
		}

		OrbitCameraController {
			camera: camera
		}

		components: [
			RenderSettings {
				id: renderSettings

				property real pointSize: 10

				property real lineWidth: 3

				RenderSurfaceSelector {
					ClearBuffers {
						buffers : ClearBuffers.ColorDepthBuffer

						CameraSelector {
							camera: camera

							RenderStateSet {
								renderStates: [
//									ScissorTest {
//										bottom: 100
//										height: 1000
//										left: -100
//										width: 2000
//									},
									ClipPlane {
										distance: 100
										normal: Qt.vector3d( 0.0, 1.0, 0.0 )
									},

									DepthTest { depthFunction: DepthTest.Less },
									CullFace { mode: CullFace.NoCulling },
									LineWidth { value: renderSettings.lineWidth },
									PointSize { value: renderSettings.pointSize; sizeMode: PointSize.Fixed }
								]
							}
						}
					}
				}
			},
			// Event Source will be set by the Qt3DQuickWindow
			InputSettings { }
		]

		PhongMaterial {
			id: material

			diffuse: "pink"
			ambient: "darkred"
		}

		TorusMesh {
			id: torusMesh
			radius: 5
			minorRadius: 1
			rings: 100
			slices: 20
		}

		Transform {
			id: torusTransform
			scale3D: Qt.vector3d(1.5, 1, 0.5)
			rotation: fromAxisAndAngle(Qt.vector3d(1, 0, 0), 45)
		}

		Entity {
			id: torusEntity
			components: [ torusMesh, material, torusTransform ]
		}

		SphereMesh {
			id: sphereMesh
			radius: 3
		}

		Transform {
			id: sphereTransform
			property real userAngle: 0.0
			matrix: {
				var m = Qt.matrix4x4();
				m.rotate(userAngle, Qt.vector3d(0, 1, 0));
				m.translate(Qt.vector3d(20, 0, 0));
				return m;
			}
		}

		QQ2.NumberAnimation {
			target: sphereTransform
			property: "userAngle"
			duration: 10000
			from: 0
			to: 360

			loops: QQ2.Animation.Infinite
			running: true
		}

		Entity {
			id: sphereEntity
			components: [ sphereMesh, material, sphereTransform ]
		}
	}
}

