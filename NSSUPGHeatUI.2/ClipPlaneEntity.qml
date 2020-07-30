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

		// Multiply normal vector by rotation matrix (translation components are subtracted, uniform scale does not matter).
//		onRotationChanged: normal = matrix.minus(Qt.matrix4x4(0, 0, 0, translation.x,
//															  0, 0, 0, translation.y,
//															  0, 0, 0, translation.z,
//															  0, 0, 0, 0)).times(Qt.vector3d(0.0, -1.0, 0.0))

//		{
//			var sinAlpha = Math.sin(rotationZ * Math.PI / 180)
//			var cosAlpha = Math.cos(rotationZ * Math.PI / 180)
//			var sinBeta = Math.sin(rotationY * Math.PI / 180)
//			var cosBeta = Math.cos(rotationY * Math.PI / 180)
//			var sinGamma = Math.sin(rotationX * Math.PI / 180)
//			var cosGamma = Math.cos(rotationX * Math.PI / 180)
//			var rotMatrix = Qt.matrix4x4(cosAlpha * cosBeta, cosAlpha * sinBeta * sinGamma - sinAlpha * cosGamma, cosAlpha * sinBeta * cosGamma + sinAlpha * sinGamma, 0,
//										 sinAlpha * cosBeta, sinAlpha * sinBeta * sinGamma + cosAlpha * cosGamma, sinAlpha * sinBeta * cosGamma - cosAlpha * sinGamma, 0,
//										 - sinBeta, cosBeta * sinGamma, cosBeta * cosGamma, 0,
//										 0, 0, 0, 1)

//			var a = Qt.matrix4x4(10,1,0,0,
//								 0,0,0,0,
//								 0,0,0,0,
//								 0,0,0,1);
//			var b = Qt.vector3d(1,1,1);
//			var bb = Qt.vector4d(2,2,2, 0);
//			var c = a.times(b);
//			var cc = a.times(bb);
//			console.log(c.toString()); // QVector3D(0.155556, 0.437037, 0.718518)
//			console.log(cc.toString()); // QVector3D(0.155556, 0.437037, 0.718518)
//			console.log(rotation.x, rotation.y, rotation.z, normal)
//			console.log(Math.sin(90))
//			console.log(rotMatrix.times(normal).toString())
//			normal = rotMatrix.times(Qt.vector3d(0.0, -1.0, 0.0))
//			normal = matrix.times(Qt.vector3d(0.0, -1.0, 0.0))
//			console.log(rotationX, rotationY, rotationZ, normal)
//			console.log(matrix)
//			console.log(rotMatrix)
//		}
    }

	components: [material, mesh, transform]
}

