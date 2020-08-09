import QtQuick 2.0

import Qt3D.Core 2.0
import Qt3D.Render 2.0

Material {
	id: root

	property ShaderData clipPlanesData: ShaderData {
		property int count: 0
	}

	property color ambient:  Qt.rgba(0.05, 0.05, 0.05, 1.0)

	property color diffuse:  Qt.rgba(0.7, 0.7, 0.7, 1.0)

	property color specular: Qt.rgba(0.95, 0.95, 0.95, 1.0)

	property real shininess: 1.0

	property real alpha: 1.0

	property RenderState activeDepthState: DepthTest { depthFunction: DepthTest.Less }

	property BlendEquationArguments blendEquationArguments: BlendEquationArguments {
		sourceRgb: BlendEquationArguments.SourceAlpha
		sourceAlpha: BlendEquationArguments.SourceAlpha
		destinationRgb: BlendEquationArguments.OneMinusSourceAlpha
		destinationAlpha: BlendEquationArguments.OneMinusSourceAlpha
	}

	property BlendEquation blendEquation: BlendEquation {
		blendFunction: BlendEquation.Add
	}

	property CullFace cullFace: CullFace {
		mode: CullFace.Back
	}

	effect: Effect {
		parameters: [
			Parameter { name: "ka"; value: Qt.vector3d(root.ambient.r, root.ambient.g, root.ambient.b) },
			Parameter { name: "kd"; value: Qt.vector3d(root.diffuse.r, root.diffuse.g, root.diffuse.b) },
			Parameter { name: "ks"; value: Qt.vector3d(root.specular.r, root.specular.g, root.specular.b) },
			Parameter { name: "alpha"; value: root.alpha },
			Parameter { name: "shininess"; value: root.shininess },
			Parameter { name: "clipPlanesData"; value: root.clipPlanesData }
		]

		FilterKey {
			id: forward

			name: "renderingStyle"
			value: "forward"
		}

		techniques: [
			Technique {
				filterKeys: [forward]

				graphicsApiFilter {
					api: GraphicsApiFilter.OpenGL
					profile: GraphicsApiFilter.CoreProfile
					majorVersion: 3
					minorVersion: 2
				}

				renderPasses: [
					RenderPass {
						renderStates: [
							cullFace,
							blendEquationArguments,
							blendEquation,
							activeDepthState
						]

						shaderProgram: ShaderProgram {
							vertexShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/phongalphaclip.vert"))
							fragmentShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/phongalphaclip.frag"))
						}
					}
				]
			}
		]
	}
}
