import QtQuick 2.0

import Qt3D.Core 2.0
import Qt3D.Render 2.0

Material {
	id: root

	property alias texture: diffuseTextureParameter.value

	parameters: [
		Parameter {
			id: diffuseTextureParameter

			name: "diffuseTexture"
		}
	]

	effect: Effect {
		FilterKey {
			id: forward
			name: "renderingStyle"
			value: "forward"
		}

		ShaderProgram {
			id: gl3Shader

			vertexShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/transparenttexture.vert"))
			fragmentShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/transparenttexture.frag"))
		}

		ShaderProgram {
			id: es2Shader

			vertexShaderCode: loadSource(Qt.resolvedUrl("shaders/es2/transparenttexture.vert"))
			fragmentShaderCode: loadSource(Qt.resolvedUrl("shaders/es2/transparenttexture.frag"))
		}

		techniques: [
			// OpenGL 3.1
			Technique {
				filterKeys: [forward]

				graphicsApiFilter {
					api: GraphicsApiFilter.OpenGL
					profile: GraphicsApiFilter.CoreProfile
					majorVersion: 3
					minorVersion: 1
				}

				renderPasses: RenderPass {
					renderStates: [
						CullFace {
							mode: CullFace.Back
						},
						DepthTest {
							depthFunction: DepthTest.Less
						},
//						NoDepthMask {},							//??
						BlendEquationArguments {
							sourceRgb: BlendEquationArguments.SourceAlpha
							sourceAlpha: BlendEquationArguments.SourceAlpha
							destinationRgb: BlendEquationArguments.OneMinusSourceAlpha
							destinationAlpha: BlendEquationArguments.OneMinusSourceAlpha
						},
						BlendEquation {
							blendFunction: BlendEquation.Add
						}
					]

					shaderProgram: gl3Shader
				}
			},

			// OpenGL 2.0
			Technique {
				filterKeys: [forward]

				graphicsApiFilter {
					api: GraphicsApiFilter.OpenGL
					profile: GraphicsApiFilter.NoProfile
					majorVersion: 2
					minorVersion: 0
				}

				renderPasses: RenderPass {
					renderStates: [
						CullFace {
							mode: CullFace.Back
						},
						DepthTest {
							depthFunction: DepthTest.Less
						},
						NoDepthMask {},
						BlendEquationArguments {
							sourceRgb: BlendEquationArguments.SourceAlpha
							sourceAlpha: BlendEquationArguments.SourceAlpha
							destinationRgb: BlendEquationArguments.OneMinusSourceAlpha
							destinationAlpha: BlendEquationArguments.OneMinusSourceAlpha
						},
						BlendEquation {
							blendFunction: BlendEquation.Add
						}
					]
					shaderProgram: es2Shader
				}
			},

			// ES 2.0
			Technique {
				filterKeys: [forward]

				graphicsApiFilter {
					api: GraphicsApiFilter.OpenGLES
					profile: GraphicsApiFilter.CoreProfile
					majorVersion: 2
					minorVersion: 0
				}

				renderPasses: RenderPass {
					renderStates: [
						CullFace {
							mode: CullFace.Back
						},
						DepthTest {
							depthFunction: DepthTest.Less
						},
						BlendEquationArguments {
							sourceRgb: BlendEquationArguments.SourceAlpha
							sourceAlpha: BlendEquationArguments.SourceAlpha
							destinationRgb: BlendEquationArguments.OneMinusSourceAlpha
							destinationAlpha: BlendEquationArguments.OneMinusSourceAlpha
						},
						BlendEquation {
							blendFunction: BlendEquation.Add
						}
					]
					shaderProgram: es2Shader
				}
			}
		]
	}
}
