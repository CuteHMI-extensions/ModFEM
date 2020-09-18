import Qt3D.Core 2.0
import Qt3D.Render 2.14
import Qt3D.Extras 2.14

Entity {
    id: root

	property alias transform: transform

	property alias width: mesh.width

	property alias height: mesh.height

	property vector3d normal: Qt.vector3d(0.0, -1.0, 0.0)

	readonly property vector4d equation: Qt.vector4d(normal.x, normal.y, normal.z, -(normal.x * transform.translation.x + normal.y * transform.translation.y + normal.z * transform.translation.z))

	PhongAlphaClipMaterial {
		id: material

		ambient: "yellow"
		diffuse: "yellow"
		specular: "yellow"
		shininess: 0.0
	}

	GeometryRenderer {
		id: mesh

		property real width: 10.0

		property real height: 5.0

		primitiveType: GeometryRenderer.LineLoop
		instanceCount: 1
		geometry: Geometry {
			Attribute {
				name: defaultPositionAttributeName
				attributeType: Attribute.VertexAttribute
				vertexBaseType: Attribute.Float
				vertexSize: 3
				byteOffset: 0
				byteStride: vertexSize * Float32Array.BYTES_PER_ELEMENT
				count: buffer.floatArr.length / vertexSize
				buffer: Buffer {
					data: floatArr
					property var floatArr: new Float32Array([-width * 0.5, 0, height * 0.5,
															 width * 0.5, 0, height * 0.5,
															 width * 0.5, 0, -height * 0.5,
															 -width * 0.5, 0, -height * 0.5
															])
				}
			}
		}
	}

	// Tranform component used to obtain rotation-only matrix.
	Transform {
		rotation: transform.rotation

		onRotationChanged: normal = matrix.times(Qt.vector3d(0.0, -1.0, 0.0))
	}

    Transform {
        id: transform
    }

	components: [material, mesh, transform]
}

