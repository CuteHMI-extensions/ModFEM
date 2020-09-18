import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Extras 2.14

Entity {
	id: root

	property alias transform: transform

	components: [transform, mesh, material]

	PhongAlphaClipMaterial {
		id: material

		ambient: "silver"
	}

	Transform {
		id: transform
	}

	CuboidMesh {
		id: mesh

		xExtent: 20
		yExtent: 6
		zExtent: 6
	}
}
