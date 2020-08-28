import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property AbstractColorMapper triangle: HueColorMapper {
		valueBegin: 0.0
		valueEnd: 10
		input: elementData.triangleFields.velocities ? elementData.triangleFields.velocities : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper quadTriangle: HueColorMapper {
		valueBegin: 0.0
		valueEnd: 10
		input: elementData.quadFields.triangleVelocities ? elementData.quadFields.triangleVelocities : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}
}
