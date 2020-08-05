import Qt3D.Core 2.0
import Qt3D.Render 2.14
import Qt3D.Extras 2.14

/** @todo remove */
Entity {
	id: root

	property matrix4x4 mat

	PhongMaterial {
		id: material

		ambient: "red"
		diffuse: "red"
	}

	CylinderMesh {
		id: mesh

		length: 10
		radius: 1
	}

	Transform {
		id: transform

		matrix: root.mat
//		translation: (0, 0, 10)
	}

	components: [material, mesh]
}

