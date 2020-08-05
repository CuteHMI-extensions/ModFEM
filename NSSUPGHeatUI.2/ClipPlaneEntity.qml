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

	PhongAlphaMaterial {
		id: material

		ambient: "red"
		diffuse: "red"
	}

	PlaneMesh {
        id: mesh

		width: 20.0
        height: 20.0
        meshResolution: Qt.size(2, 2)
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

