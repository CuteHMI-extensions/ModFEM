import Qt3D.Core 2.0
import Qt3D.Render 2.0
import Qt3D.Extras 2.14


import QtQuick 2.0 as QQ2

Entity {
	id: root

	property alias transform: transform

	property alias running: rotorAnimation.running

	components: [transform, rotor, motor]

	Transform {
		id: transform

		scale: 0.25
		rotationY: 180
		translation: Qt.vector3d(-3, 0, 0)
	}

	Entity {
		id: rotor

		PhongMaterial {
			id: rotorMaterial

			ambient: "silver"
			diffuse: "green"
			specular: "lime"
			shininess: 1.0
		}

		Transform {
			id: rotorTransform
		}

		Mesh {
			id: rotorMesh

			source: "models/rotor.obj"
		}

		QQ2.NumberAnimation {
			id: rotorAnimation

			target: rotorTransform
			property: "rotationX"
			duration: 2000
			from: 0
			to: 360

			loops: QQ2.Animation.Infinite
			running: true
		}

		components: [rotorTransform, rotorMesh, rotorMaterial]
	}

	Entity {
		id: motor

		PhongMaterial {
			id: motorMaterial

			ambient: "#3c5c5c"
		}

		Mesh {
			id: motorMesh

			source: "models/motor.obj"
		}

		components: [motorMesh, motorMaterial]
	}
}
