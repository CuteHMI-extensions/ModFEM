import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property AbstractColorMapper triangle: HueColorMapper {
		valueBegin: 0.0
		valueEnd: 0.5
		input: elementData.triangleFields.velocityMagnitudes ? elementData.triangleFields.velocityMagnitudes : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper quadTriangle: HueColorMapper {
		valueBegin: 0.0
		valueEnd: 0.5
		input: elementData.quadFields.triangleVelocityMagnitudes ? elementData.quadFields.triangleVelocityMagnitudes : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}
}
