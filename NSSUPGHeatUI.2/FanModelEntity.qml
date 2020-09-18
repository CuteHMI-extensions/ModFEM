import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Extras 2.14

Entity {
	id: root

	property alias transform: transform

	components: [transform, mesh, material]

	PhongMaterial {
		id: material

		ambient: "#5c5c5c"
	}

	Transform {
		id: transform

		translation: Qt.vector3d(-10, -3, 4.2)
		scale: 0.0175
		rotationY: 180
	}

	Mesh {
		id: mesh

		source: "models/fan.obj"
	}
}
