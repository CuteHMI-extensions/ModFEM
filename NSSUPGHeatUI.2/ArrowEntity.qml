import Qt3D.Core 2.0
import Qt3D.Render 2.14
import Qt3D.Extras 2.14

Entity {
	id: root

	property real radius: 0.025
//	property real radius: 0.25

	property real headRadiusRatio: 5

	property real headAngle: 22.5

	property real lengthMultiplier: 1.0

	property vector3d vector: Qt.vector3d(0, 1, 0)
//	property vector3d vector: Qt.vector3d(1, 1, 0)

	readonly property vector3d initialVector: Qt.vector3d(0, 1, 0)	// Vector conforming with initial mesh orientation.

	property alias transform: transform

	onVectorChanged: {
		var v = initialVector.crossProduct(vector);

		if (vector.normalized().fuzzyEquals(initialVector.times(-1))) {
			arrowTransform.rotationX = 180
			arrowTransform.rotationY = 0
			arrowTransform.rotationZ = 0
		} else {
			// q.w = vector.length() * initialVector.length() + vector.dotProduct(initialVector)
			var q = Qt.vector4d(v.x, v.y, v.z, vector.length() + vector.y).normalized()
			arrowTransform.rotation = Qt.quaternion(q.w, q.x, q.y, q.z)
		}

		arrowTransform.translation = vector.times(0.5 * lengthMultiplier)
	}

	Entity {
		id: arrowEntity

		PhongAlphaClipMaterial {
			id: material

			ambient: "red"
			diffuse: "red"
			alpha: 1.0
		}

		Entity {
			id: cylinderEntity

			CylinderMesh {
				id: cylinder

				length: vector.length() * lengthMultiplier
				radius: root.radius

//				onLengthChanged: console.log(length)
			}

			components: [cylinder, material]
		}

		Entity {
			id: coneEntity

			ConeMesh {
				id: cone

				length: bottomRadius / Math.tan(headAngle * Math.PI /180)
				bottomRadius: root.radius * root.headRadiusRatio
			}

			Transform {
				id: coneTransform

				translation.y: (cone.length + cylinder.length) * 0.5
			}

			components: [cone, coneTransform, material]
		}

		Transform {
			id: arrowTransform
		}

		components: [cylinderEntity, coneEntity, arrowTransform]
	}

	Transform {
		id: transform
	}

	components: [arrowEntity, transform]
}

