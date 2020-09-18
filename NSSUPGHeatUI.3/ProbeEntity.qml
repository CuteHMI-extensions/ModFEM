import Qt3D.Core 2.14
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

		activeDepthState: DepthTest { depthFunction: DepthTest.Always }
	}

	SphereMesh {
		id: mesh

		radius: 0.05
	}

	Transform {
		id: transform
	}

	components: [material, mesh, transform]
}

