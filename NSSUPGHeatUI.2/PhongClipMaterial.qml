import QtQuick 2.0

import Qt3D.Core 2.0
import Qt3D.Render 2.0

Material {
	id: root

	property ShaderData clipPlanesData

	effect: Effect {
		parameters: [
			Parameter { name: "ka"; value: "black" },
			Parameter { name: "kd"; value: "blue" },
			Parameter { name: "ks"; value: "white" },
			Parameter { name: "shininess"; value: 100 },
			Parameter { name: "lightPosition"; value: Qt.vector4d(1.0, 1.0, 0.0, 1.0) },
			Parameter { name: "lightIntensity"; value: Qt.vector3d(1.0, 1.0, 1.0) },
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
							vertexShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/phongclip.vert"))
							fragmentShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/phongclip.frag"))
						}
					}
				]
			}
		]
	}
}
