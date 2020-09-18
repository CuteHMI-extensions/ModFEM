import QtQuick 2.0 as QQ2

import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Input 2.5
import Qt3D.Extras 2.5

Entity {
	id: root

	components: [renderSettings, inputSettings]

	property alias pointSize: renderSettings.pointSize

	property alias lineWidth: renderSettings.lineWidth

	property alias cameraRotationMatrix: cameraRotationTransform.matrix

	RenderSettings {
		id: renderSettings

		property real pointSize: 10

		property real lineWidth: 3

//		ForwardRenderer {
//			clearColor: "grey"
//			camera: camera1
//		}

		RenderSurfaceSelector {
			ClearBuffers {
				buffers : ClearBuffers.ColorDepthBuffer
				clearColor: "transparent"

				CameraSelector {
					camera: camera1

					RenderStateSet {
						renderStates: [
							ClipPlane { planeIndex: 0 },
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
		nearPlane : 0.025
		farPlane : 1000.0
		position: Qt.vector3d( 0.0, 0.0, 10.0 )
		upVector: Qt.vector3d( 0.0, 1.0, 0.0 )
		viewCenter: Qt.vector3d( 0.0, 0.0, 0.0 )
	}

	OrbitCameraController {
		camera: camera1

		linearSpeed: 100
	}

	Transform {
		id: cameraRotationTransform

		rotation: camera1.transform.rotation
	}
}
