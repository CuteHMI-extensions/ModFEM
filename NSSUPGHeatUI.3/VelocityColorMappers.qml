import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property real minimum: 0.0

	property real maximum: 0.5

	property AbstractColorMapper triangle: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.triangleFields.velocityMagnitudes ? elementData.triangleFields.velocityMagnitudes : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper quadTriangle: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.quadFields.triangleVelocityMagnitudes ? elementData.quadFields.triangleVelocityMagnitudes : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper cap: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.capFields.velocityMagnitudes ? elementData.capFields.velocityMagnitudes : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}
}
