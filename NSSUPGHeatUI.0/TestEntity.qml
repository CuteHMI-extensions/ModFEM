import QtQuick 2.0 as QQ2
import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Extras 2.5

Entity {
	id: root

	property alias nodesEnabled: nodeRenderer.enabled

	property bool facesEnabled: true

	property real rotationX: 0

	property real rotationY: 0

	property real rotationZ: 0

	property real scale: 1

	PhongMaterial {
		id: lineMaterial

		diffuse: Qt.rgba(0.5, 0.2, 0.1, 1.0)
		ambient: Qt.rgba(0.6, 0.2, 0.2, 1.0)
		specular: Qt.rgba(0.8, 0.2, 0.2, 1.0)
		shininess: 1.0
	}

	PhongMaterial {
		id: faceMaterial

		//		ambient: "red"
		diffuse: Qt.rgba(0.1, 0.5, 0.1, 1.0)
		ambient: Qt.rgba(0.2, 0.6, 0.2, 1.0)
		specular: Qt.rgba(0.2, 0.8, 0.2, 1.0)
		shininess: 1.0
	}

	GeometryRenderer {
		id: nodeRenderer

		primitiveType: GeometryRenderer.Points
		instanceCount: 1
		geometry: Geometry {
			Attribute {
				name: defaultPositionAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Double
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
				count: heatProblem.mesh.nodeCount
				buffer: Buffer {
					data: heatProblem.mesh.nodeCoords
				}
			}
		}
	}

	GeometryRenderer {
		id: testRenderer

		primitiveType: GeometryRenderer.Triangles
		instanceCount: 1
		geometry: Geometry {
			Attribute {
				name: defaultNormalAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Double
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
				count: buffer.floatArr.length / vertexSize
				buffer: Buffer {
					//					type: Buffer.VertexBuffer
					data: floatArr
					property var floatArr: new Float64Array([
																// #1 triangle
																1, 0, 0,
																1, 0, 0,
																1, 0, 0,
																// #2 triangle
																0, 1, 0,
																0, 1, 0,
																0, 1, 0
															])
				}
				//				QQ2.Component.onCompleted: console.log("buffer size " + buffer.floatArr.length)

				//				QQ2.Component.onCompleted: console.log(buffer.data)
			}

			Attribute {
				name: defaultPositionAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Double
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
				count: buffer.floatArr.length / vertexSize
				buffer: Buffer {
					//					type: Buffer.VertexBuffer
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
				//				QQ2.Component.onCompleted: console.log("buffer size " + buffer.floatArr.length)

				//				QQ2.Component.onCompleted: console.log(buffer.data)
			}
		}
	}

	GeometryRenderer {
		id: triangleFaceRenderer

		enabled: facesEnabled

		primitiveType: GeometryRenderer.Triangles
		instanceCount: 1
		geometry: Geometry {
			Attribute {
				name: defaultNormalAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Double
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
				count: heatProblem.mesh.faceData.triangleCount * 3
				buffer: Buffer {
					data: heatProblem.mesh.faceData.triangleNormals
				}
			}

			Attribute {
				name: defaultPositionAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Double
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
				count: heatProblem.mesh.faceData.triangleCount * 3
				buffer: Buffer {
					data: heatProblem.mesh.faceData.triangleCoords
				}
			}
		}
	}

	GeometryRenderer {
		id: trianglesRenderer

		primitiveType: GeometryRenderer.Triangles
		instanceCount: 1
		geometry: Geometry {
			Attribute {
				name: defaultNormalAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Double
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
				count: heatProblem.mesh.triangles.count * 3			// 3 vertices per triangle.
				buffer: Buffer {
					data: heatProblem.mesh.triangles.normals
				}
			}

			Attribute {
				name: defaultPositionAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Double
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
				count: heatProblem.mesh.triangles.count * 3		// 3 vertices per triangle.
				buffer: Buffer {
					data: heatProblem.mesh.triangles.coords
				}
			}
		}
	}

	TorusMesh {
		id: torusMesh

		radius: 5
		minorRadius: 1
		rings: 100
		slices: 20

		//		primitiveType: GeometryRenderer.LineLoop
	}

	Transform {
		id: torusTransform
		scale3D: Qt.vector3d(1.5, 1, 0.5)
		rotation: fromAxisAndAngle(Qt.vector3d(1, 0, 0), 45)
		translation: Qt.vector3d(10, 2, -10)
	}

	Transform {
		id: trianglesTransforms
		scale: root.scale
		rotationX: root.rotationX
		rotationY: root.rotationY
		rotationZ: root.rotationZ
	}

	//		components: [torusMesh, faceMaterial, torusTransform]
	//		components: [triangleFaceRenderer, faceMaterial]
	components: [trianglesRenderer, faceMaterial, trianglesTransforms]
	//		components: [nodeRenderer, faceMaterial]
	//	components: [nodeRenderer, triangleFaceRenderer, faceMaterial]
	//		components: [testRenderer, faceMaterial]
}
