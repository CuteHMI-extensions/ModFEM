import QtQuick 2.0
import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Input 2.0
import Qt3D.Extras 2.0

Entity {
	id: plane

	GeometryRenderer {
		id: geometry
		geometry: Geometry {
			boundingVolumePositionAttribute: position

			Attribute {
				id: position
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Float
				vertexSize: 3
				count: 4
				byteOffset: 0
				byteStride: 8 * 4
//				name: "position"
				name: defaultPositionAttributeName
				buffer: vertexBuffer

//				Component.onCompleted: console.log(name, defaultColorAttributeName, defaultTextureCoordinateAttributeName)
			}

			Attribute {
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Float
				vertexSize: 3
				count: 4
				byteOffset: 3 * 4
				byteStride: 8 * 4
//				name: "color"
				name: defaultColorAttributeName
				buffer: vertexBuffer
			}

			Attribute {
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Float
				vertexSize: 2
				count: 4
				byteOffset: 6 * 4
				byteStride: 8 * 4
//				name: "texCoord"
				name: defaultTextureCoordinateAttributeName
				buffer: vertexBuffer
			}

			Attribute {
				attributeType: Attribute.IndexAttribute
				vertexBaseType: Attribute.UnsignedShort
				vertexSize: 1
				count: 6
				buffer: Buffer {
					type: Buffer.IndexBuffer
					data: new Uint16Array([
						0, 1, 3,  // First Triangle
						1, 2, 3,  // Second Triangle
					])
				}
			}
		}

		Buffer {
			id: vertexBuffer
			type: Buffer.VertexBuffer
			data: new Float32Array([
				// Positions	  // Colors		// Texture Coords
				 0.5,  0.5, -1.0,   1.0, 0.0, 0.0,   1.0, 1.0, // Top Right
				 0.5, -0.5, -1.0,   0.0, 1.0, 0.0,   1.0, 0.0, // Bottom Right
				-0.5, -0.5, -1.0,   0.0, 0.0, 1.0,   0.0, 0.0, // Bottom Left
				-0.5,  0.5, -1.0,   1.0, 1.0, 0.0,   0.0, 1.0, // Top Left
			])
		}
	}

	Material {
		id: material
		effect: Effect {
			FilterKey {
				id: forward
				name: "renderingStyle"
				value: "forward"
			}

			techniques: Technique {
				filterKeys: [forward]

				graphicsApiFilter {
					api: GraphicsApiFilter.OpenGL
					profile: GraphicsApiFilter.CoreProfile
					majorVersion: 3
					minorVersion: 1
				}

				renderPasses: RenderPass {
					renderStates: CullFace { mode: CullFace.NoCulling }

					shaderProgram: ShaderProgram {
						vertexShaderCode: loadSource(Qt.resolvedUrl("shaders/texture.vert"))
						fragmentShaderCode: loadSource(Qt.resolvedUrl("shaders/texture.frag"))
					}
				}
			}
		}

		parameters: Parameter {
			name: "ourTexture"
			value: Texture2D { // #TRYIT: change filtering mode
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
		}
	}

	CustomMaterial {
		id: customMaterial
	}

	PhongMaterial {
		id: phongMaterial
		diffuse: Qt.rgba(0.7, 0.4, 0.4, 1.0)
		ambient: Qt.rgba(0.3, 0.3, 0.3, 1.0)
		shininess: 1.0
	}

	SphereMesh {
		id: sphereMesh
	}

//	components: [geometry, material]
	components: [sphereMesh, customMaterial]
}
