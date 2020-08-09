import Qt3D.Core 2.0
import Qt3D.Render 2.14
import Qt3D.Extras 2.14

Entity {
	id: root

	property alias transform: transform

	PhongAlphaClipMaterial {
		id: material

		ambient: "red"
		diffuse: "red"
		alpha: 1.0
	}

	SphereMesh {
		id: mesh

		radius: 0.25
	}

	Transform {
		id: transform
	}

	components: [material, mesh, transform]
}

