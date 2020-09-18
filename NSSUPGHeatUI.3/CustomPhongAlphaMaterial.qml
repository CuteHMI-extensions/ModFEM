import QtQuick 2.0

import Qt3D.Core 2.0
import Qt3D.Render 2.0

/// @todo remove.
Material {
	effect: Effect {

		parameters: [
			Parameter { name: "alpha";  value: 0.5 },
			Parameter { name: "ka"; value: "black" },
			Parameter { name: "kd"; value: "blue" },
			Parameter { name: "ks"; value: "white" },
			Parameter { name: "shininess"; value: 100 },
			Parameter { name: "lightPosition"; value: Qt.vector4d(100.0, 100.0, 0.0, 1.0) },
			Parameter { name: "lightIntensity"; value: Qt.vector3d(1.0, 1.0, 1.0) }
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
							vertexShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/phongalpha.vert"))
							fragmentShaderCode: loadSource(Qt.resolvedUrl("shaders/gl3/phongalpha.frag"))
						}
					}
				]
			}
		]
	}
}
