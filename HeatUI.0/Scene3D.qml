import QtQuick 2.0 as QQ2
import QtQuick.Scene3D 2.0

import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Input 2.5
import Qt3D.Extras 2.5

import ModFEM.QtHeat 0.0 as QtHeat

Scene3D {
	id: scene3d

	aspects: ["input", "logic"]
	cameraAspectRatioMode: Scene3D.AutomaticAspectRatio

	//	property QtHeat.Controller controller

	Entity {
		id: sceneRoot

		Camera {
			id: camera1
			projectionType: CameraLens.PerspectiveProjection
			fieldOfView: 45
			nearPlane : 0.1
			farPlane : 1000.0
			position: Qt.vector3d( 0.0, 0.0, 40.0 )
			upVector: Qt.vector3d( 0.0, 1.0, 0.0 )
			viewCenter: Qt.vector3d( 0.0, 0.0, 0.0 )
		}

		FirstPersonCameraController {
			camera: camera1
			linearSpeed: 100
		}

		components: [
			RenderSettings {
				activeFrameGraph: RenderSurfaceSelector {
					ClearBuffers {
						buffers : ClearBuffers.ColorDepthBuffer

						CameraSelector {
							camera: camera1

							RenderStateSet {
								renderStates: [
									LineWidth { value: 3 },
									PointSize { value: 10; sizeMode: PointSize.Fixed }
								]
							}
						}
					}
				}
			},

			InputSettings { }
		]


		//				components: [
		//					RenderSettings {
		//						activeFrameGraph: ForwardRenderer {
		//							clearColor: "transparent"

		//							camera: camera

		//							TechniqueFilter {
		//								RenderStateSet {
		//									renderStates: [
		//										LineWidth { value: 12 },	//hehe
		//										PointSize { value: 10 }
		//									]
		//								}}
		//							}
		//					},

		//					InputSettings { }
		//				]

//		components: [
//			RenderSettings {
//				activeFrameGraph: RenderSurfaceSelector {
//					Viewport {
//						normalizedRect: Qt.rect(0.0, 0.0, 1.0, 1.0)

//						ClearBuffers {
//							buffers: ClearBuffers.ColorDepthBuffer

//							RenderStateSet {
//								renderStates: [
//									LineWidth { value: 3 },	//hehe
//									PointSize { value: 10; sizeMode: PointSize.Fixed }
//								]
//							}
//						}
//					}
//				}
//			},
//			InputSettings {}
//		]


		MeshEntity {
			//			Controller: controller.buffer
		}
	}
}
