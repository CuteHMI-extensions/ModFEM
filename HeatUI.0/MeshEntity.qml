import QtQuick 2.0 as QQ2
import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Extras 2.5

Entity {
	id: root

	PhongMaterial {
		id: material

		ambient: "red"
	}

	GeometryRenderer {
		id: renderer

		primitiveType: GeometryRenderer.LineLoop
		instanceCount: 1
		geometry: Geometry {
			Attribute {
				name: defaultPositionAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Double
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
//				count: buffer.floatArr.length / vertexSize
				count: 18 / vertexSize //temp
//				buffer: heat.buffer
				buffer: Buffer {
//					type: Buffer.VertexBuffer
//					data: heatProblem.meshData
					data: floatArr
					property var floatArr: new Float64Array([
																// #1 triangle
																-5, -5, 0,
																0, -5, 0,
																-5, 5, 0,
																// #2 triangle
																0, -5, 0,
																10, -5, 0,
																5, 5, 0
															])
				}
				QQ2.Component.onCompleted: console.log("buffer size " + buffer.floatArr.length)

//				QQ2.Component.onCompleted: console.log(buffer.data)

			}
		}
	}
	components: [renderer, material]
}
