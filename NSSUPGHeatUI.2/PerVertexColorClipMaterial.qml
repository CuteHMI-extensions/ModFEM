import QtQuick 2.0

import Qt3D.Core 2.0
import Qt3D.Render 2.0

Material {
	property ShaderData clipPlanesData

	effect: Effect {
		parameters: [
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
