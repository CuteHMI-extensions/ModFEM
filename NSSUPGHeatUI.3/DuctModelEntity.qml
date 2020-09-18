import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Extras 2.14

Entity {
	id: root

	property alias transform: transform

	property ShaderData clipPlanesData

	components: [transform, mesh, material]

	PhongAlphaClipMaterial {
		id: material

		ambient: "silver"

		clipPlanesData: root.clipPlanesData
		cullFace: CullFace {
			mode: CullFace.NoCulling
		}
	}

	Transform {
		id: transform

		translation: Qt.vector3d(0, 0.0025, 0.0025)
	}

	CuboidMesh {
		id: mesh

		xExtent: 10
		yExtent: 1.25
		zExtent: 1.25
	}
}
