import QtQuick 2.0

import Qt3D.Core 2.0
import Qt3D.Render 2.0

Material {
	id: root

	property ShaderData clipPlanesData: ShaderData {
		property int count: 0
	}

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
			Parameter { name: "clipPlanesData"; value: root.clipPlanesData },
			Parameter { name: "alpha"; value: root.alpha }
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
							vertexShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/pervertexcolorclip.vert"))
							fragmentShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/pervertexcolorclip.frag"))
						}
					}
				]
			}
		]
	}
}
