import QtQuick 2.0 as QQ2
import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Extras 2.5

import ModFEM.NSSUPGHeat 1.0 as NSSUPGHeat

Entity {
	id: root

	property NSSUPGHeat.ElementData elementData

	property NSSUPGHeat.AbstractColorMapper triangleColorMapper

	property NSSUPGHeat.AbstractColorMapper quadTriangleColorMapper

	property NSSUPGHeat.AbstractColorMapper capColorMapper

	property bool nodesEnabled: false

	property bool linesEnabled: true

	property bool facesEnabled: true

	property bool capsEnabled: true

	property real alpha: 0.3

	property alias transform: transform

	property vector3d maxExtent: lineRenderer.geometry.maxExtent

	property vector3d minExtent: lineRenderer.geometry.minExtent

	property ShaderData clipPlanesData

	PhongAlphaClipMaterial {
		id: nodeMaterial

		ambient: "red"
		diffuse: "red"
		shininess: 0.0

		clipPlanesData: root.clipPlanesData
	}

	PhongAlphaClipMaterial {
		id: lineMaterial

		ambient: "white"
		diffuse: "white"
		shininess: 0.0
		alpha: 1.0

		clipPlanesData: root.clipPlanesData
	}

	//	PhongAlphaClipMaterial {
	//		ambient: "green"
	//		diffuse: "green"

	VertexAlphaClipMaterial {
		id: faceMaterial

		clipPlanesData: root.clipPlanesData

		alpha: root.alpha
		cullFace: CullFace {
			mode: CullFace.NoCulling
		}
	}

//	PhongMaterial {
	VertexAlphaClipMaterial {
		id: capsMaterial

//		ambient: "red"
//		diffuse: "red"
		alpha: root.alpha

		cullFace: CullFace {
			mode: CullFace.NoCulling
		}
	}

	Transform {
		id: transform
	}

	property var bcPalette: [
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.hsla(0.0, 1.0, 0.35, 1.0),
		Qt.hsla(0.125, 1.0, 0.35, 1.0),
		Qt.hsla(0.25, 1.0, 0.35, 1.0),
		Qt.hsla(0.375, 1.0, 0.35, 1.0),
		Qt.hsla(0.5, 1.0, 0.35, 1.0),
		Qt.hsla(0.625, 1.0, 0.35, 1.0),
		Qt.hsla(0.75, 1.0, 0.35, 1.0),
		Qt.hsla(0.875, 1.0, 0.35, 1.0)
	]

	NSSUPGHeat.HueColorMapper {
		id: triangleTemperaturesColorMapper

		valueBegin: -65.0
		valueEnd: 300
		input: elementData.triangleFields.temperatures ? elementData.triangleFields.temperatures : new ArrayBuffer
		inputType: NSSUPGHeat.HueColorMapper.INPUT_DOUBLE
	}

	NSSUPGHeat.HueColorMapper {
		id: quadTriangleTemperaturesColorMapper

		valueBegin: -65.0
		valueEnd: 300
		input: elementData.quadFields.triangleTemperatures ? elementData.quadFields.triangleTemperatures : new ArrayBuffer
		inputType: NSSUPGHeat.HueColorMapper.INPUT_DOUBLE
	}

	NSSUPGHeat.PaletteColorMapper {
		id: trianglesPaletteColorMapper

		palette: bcPalette
		input: elementData.triangles.boundaryConditions ? elementData.triangles.boundaryConditions : new ArrayBuffer
	}

	NSSUPGHeat.PaletteColorMapper {
		id: quadsTrianglePaletteColorMapper

		palette: bcPalette
		input: elementData.quads.triangleBoundaryConditions ? elementData.quads.triangleBoundaryConditions : new ArrayBuffer
	}

	NSSUPGHeat.PaletteColorMapper {
		id: quadsPaletteColorMapper

		palette: bcPalette
		input: elementData.quads.quadBoundaryConditions ? elementData.quads.quadBoundaryConditions : new ArrayBuffer
	}

	Entity {
		id: nodesEntity

//		GeometryRenderer {
//			id: nodeRenderer

//			enabled: root.nodesEnabled

//			primitiveType: GeometryRenderer.Points
//			instanceCount: 1
//			geometry: Geometry {
//				Attribute {
//					name: defaultPositionAttributeName
//					attributeType: Attribute.VertexAttribute
//					vertexBaseType: Attribute.Float
//					vertexSize: 3
//					byteOffset: 0
//					byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
//					count: elementData.nodes.count
//					buffer: Buffer {
//						data: elementData.nodes.coords
//					}
//				}
//			}
//		}

		GeometryRenderer {
			id: triangleNodeRenderer

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
					count: elementData.triangles.count * 3		// 3 vertices per triangle.
					buffer: Buffer {
						data: elementData.triangles.coords
					}
				}

				Attribute {
					name: defaultColorAttributeName
					attributeType: Attribute.VertexAttribute
					vertexBaseType: Attribute.Float
					vertexSize: 3
					byteOffset: 0
					byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
					count: trianglesPaletteColorMapper.count
					buffer: Buffer {
						data: trianglesPaletteColorMapper.output
					}
				}
			}
		}

		GeometryRenderer {
			id: quadNodeRenderer

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
					count: elementData.quads.count * 4		// 4 vertices per quad.
					buffer: Buffer {
						data: elementData.quads.coords
					}
				}

				Attribute {
					name: defaultColorAttributeName
					attributeType: Attribute.VertexAttribute
					vertexBaseType: Attribute.Float
					vertexSize: 3
					byteOffset: 0
					byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
					count: elementData.quads.count * 4		// 4 vertices per quad.
					buffer: Buffer {
						data: quadsPaletteColorMapper.output
					}
				}
			}
		}

//		components: [nodeRenderer, nodeMaterial]
		components: [quadNodeRenderer, triangleNodeRenderer, nodeMaterial]
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
						count: quadTriangleColorMapper.count
//						count: quadTriangleTemperaturesColorMapper.count
//						count: quadsTrianglePaletteColorMapper.count
						buffer: Buffer {
//							data: quadTriangleTemperaturesColorMapper.output
							data: quadTriangleColorMapper.output
//							data: quadsTrianglePaletteColorMapper.output
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
//						count: trianglesPaletteColorMapper.count
//						count: triangleTemperaturesColorMapper.count
						count: triangleColorMapper.count
						buffer: Buffer {
							data: triangleColorMapper.output
//							data: triangleTemperaturesColorMapper.output
//							data: trianglesPaletteColorMapper.output
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

	Entity {
		id: capsEntity

		enabled: root.capsEnabled

		Entity {
			GeometryRenderer {
				id: capTrianglesRenderer

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
						count: elementData.caps.count * 3			// 3 vertices per triangle.
						buffer: Buffer {
							data: elementData.caps.normals
						}
					}

					Attribute {
						name: defaultColorAttributeName
						attributeType: Attribute.VertexAttribute
						vertexBaseType: Attribute.Float
						vertexSize: 3
						byteOffset: 0
						byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
//						count: trianglesPaletteColorMapper.count
//						count: triangleTemperaturesColorMapper.count
						count: capColorMapper.count
						buffer: Buffer {
							data: capColorMapper.output
//							data: triangleTemperaturesColorMapper.output
//							data: trianglesPaletteColorMapper.output
						}
					}

					Attribute {
						name: defaultPositionAttributeName
						attributeType: Attribute.VertexAttribute
						vertexBaseType: Attribute.Float
						vertexSize: 3
						byteOffset: 0
						byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
						count: elementData.caps.count * 3		// 3 vertices per triangle.
						buffer: Buffer {
							data: elementData.caps.coords
						}
					}
				}
			}
		}

		components: [capTrianglesRenderer, capsMaterial]
	}

	components: [transform]
}
