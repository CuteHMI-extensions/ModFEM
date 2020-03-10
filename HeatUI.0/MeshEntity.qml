import QtQuick 2.0 as QQ2
import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Extras 2.5

Entity {
	id: root

	property bool nodesEnabled: false

	property bool linesEnabled: true

	property bool facesEnabled: true

	property real rotationX: 0

	property real rotationY: 0

	property real rotationZ: 0

	property real scale: 1

	PhongMaterial {
		id: nodeMaterial

		ambient: "red"
		diffuse: "red"
		shininess: 0.0
	}

	PhongMaterial {
		id: lineMaterial

		ambient: "white"
		diffuse: "white"
		shininess: 0.0
	}

	PhongMaterial {
		id: faceMaterial

		diffuse: Qt.rgba(0.1, 0.5, 0.1, 1.0)
		ambient: Qt.rgba(0.2, 0.6, 0.2, 1.0)
		specular: Qt.rgba(0.2, 0.8, 0.2, 1.0)
		shininess: 1.0
	}

	Transform {
		id: transform

		scale: root.scale
		rotationX: root.rotationX
		rotationY: root.rotationY
		rotationZ: root.rotationZ
	}

	Entity {
		id: nodesEntity

		GeometryRenderer {
			id: nodeRenderer

			enabled: root.nodesEnabled

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
					count: heatProblem.mesh.nodes.count
					buffer: Buffer {
						data: heatProblem.mesh.nodes.coords
					}
				}
			}
		}

		components: [nodeRenderer, nodeMaterial]
	}

	Entity {
		id: linesEntity

		GeometryRenderer {
			id: lineRenderer

			enabled: root.linesEnabled

			primitiveType: GeometryRenderer.Lines
			instanceCount: 1
			geometry: Geometry {
				Attribute {
					name: defaultPositionAttributeName
					attributeType: Attribute.VertexAttribute
					vertexBaseType: Attribute.Double
					vertexSize: 3
					byteOffset: 0
					byteStride: vertexSize * Float64Array.BYTES_PER_ELEMENT
					count: heatProblem.mesh.lines.count *  2		// Each line is composed of two vertices.
					buffer: Buffer {
						data: heatProblem.mesh.lines.coords
					}
				}
			}
		}

		components: [lineRenderer, lineMaterial]
	}

	Entity {
		id: triangleEntity

		GeometryRenderer {
			id: triangleRenderer

			enabled: root.facesEnabled

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

		components: [triangleRenderer, faceMaterial]
	}

	components: [transform]
}
