import QtQuick 2.0 as QQ2
import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Extras 2.5

import ModFEM.Heat 0.0 as Heat

Entity {
	id: root

	// Move Mesh from QtHeat to ModFEM.Mesh (also refactor Qt Heat to just Heat)
	property Heat.ElementData elementData

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

//	PhongMaterial {
//		id: faceMaterial

//		diffuse: Qt.rgba(0.1, 0.5, 0.1, 1.0)
//		ambient: Qt.rgba(0.2, 0.6, 0.2, 1.0)
//		specular: Qt.rgba(0.2, 0.8, 0.2, 1.0)
//		shininess: 1.0
//	}

	PerVertexColorMaterial {
		id: faceMaterial
	}

	Transform {
		id: transform

		scale: root.scale
		rotationX: root.rotationX
		rotationY: root.rotationY
		rotationZ: root.rotationZ
	}

	property var bcColorMap: [
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.hsla(0.0, 1.0, 0.35, 1.0),
//		Qt.hsla(0.125, 1.0, 0.35, 1.0),
//		Qt.hsla(0.25, 1.0, 0.35, 1.0),
//		Qt.hsla(0.375, 1.0, 0.35, 1.0),
//		Qt.hsla(0.5, 1.0, 0.35, 1.0),
//		Qt.hsla(0.625, 1.0, 0.35, 1.0),
//		Qt.hsla(0.75, 1.0, 0.35, 1.0),
//		Qt.hsla(0.875, 1.0, 0.35, 1.0)
	]

	Heat.VertexColorMapper {
		id: trianglesVertexColorMapper

		map: bcColorMap
		colorIndices: elementData.triangles.boundaryConditions ? elementData.triangles.boundaryConditions : new ArrayBuffer
	}

	Heat.VertexColorMapper {
		id: quadsTriangleVertexColorMapper

		map: bcColorMap
		colorIndices: elementData.quads.triangleBoundaryConditions ? elementData.quads.triangleBoundaryConditions : new ArrayBuffer
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
					vertexBaseType: Attribute.Float
					vertexSize: 3
					byteOffset: 0
					byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
					count: elementData.nodes.count
					buffer: Buffer {
						data: elementData.nodes.coords
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
					vertexBaseType: Attribute.Float
					vertexSize: 3
					byteOffset: 0
					byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
					count: elementData.lines.count *  2		// Each line is composed of two vertices.
					buffer: Buffer {
						data: elementData.lines.coords
					}
				}
			}
		}

		components: [lineRenderer, lineMaterial]
	}

	Entity {
		id: facesEntity

		enabled: root.facesEnabled

		Entity {
			GeometryRenderer {
				id: quadRenderer

				primitiveType: GeometryRenderer.Triangles
				instanceCount: 1
				geometry: Geometry {
					Attribute {
						name: defaultNormalAttributeName
						attributeType: Attribute.VertexAttribute
						vertexBaseType: Attribute.Float
						vertexSize: 3
						byteOffset: 0
						byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
						count: elementData.quads.triangleCount * 3			// 3 vertices per triangle.
						buffer: Buffer {
							data: elementData.quads.triangleNormals
						}
					}

					Attribute {
						name: defaultColorAttributeName
						attributeType: Attribute.VertexAttribute
						vertexBaseType: Attribute.Float
						vertexSize: 3
						byteOffset: 0
						byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
						count: quadsTriangleVertexColorMapper.count
						buffer: Buffer {
							data: quadsTriangleVertexColorMapper.colorVertices
						}
					}

					Attribute {
						name: defaultPositionAttributeName
						attributeType: Attribute.VertexAttribute
						vertexBaseType: Attribute.Float
						vertexSize: 3
						byteOffset: 0
						byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
						count: elementData.quads.triangleCount * 3			// 3 vertices per triangle.
						buffer: Buffer {
							data: elementData.quads.triangleCoords
						}
					}
				}
			}

			components: [quadRenderer, faceMaterial]
		}

		Entity {
			GeometryRenderer {
				id: triangleRenderer

				primitiveType: GeometryRenderer.Triangles
				instanceCount: 1
				geometry: Geometry {
					Attribute {
						name: defaultNormalAttributeName
						attributeType: Attribute.VertexAttribute
						vertexBaseType: Attribute.Float
						vertexSize: 3
						byteOffset: 0
						byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
						count: elementData.triangles.count * 3			// 3 vertices per triangle.
						buffer: Buffer {
							data: elementData.triangles.normals
						}
					}

					Attribute {
						name: defaultColorAttributeName
						attributeType: Attribute.VertexAttribute
						vertexBaseType: Attribute.Float
						vertexSize: 3
						byteOffset: 0
						byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
						count: trianglesVertexColorMapper.count
						buffer: Buffer {
							data: trianglesVertexColorMapper.colorVertices
						}
					}

					Attribute {
						name: defaultPositionAttributeName
						attributeType: Attribute.VertexAttribute
						vertexBaseType: Attribute.Float
						vertexSize: 3
						byteOffset: 0
						byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
						count: elementData.triangles.count * 3		// 3 vertices per triangle.
						buffer: Buffer {
							data: elementData.triangles.coords
						}
					}
				}
			}
		}

		components: [triangleRenderer, faceMaterial]
	}

	components: [transform]
}
