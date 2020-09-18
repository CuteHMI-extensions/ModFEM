import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property real minimum: -65.0

	property real maximum: 300.0

	property AbstractColorMapper triangle: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.triangleFields.pressures ? elementData.triangleFields.pressures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper quadTriangle: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.quadFields.trianglePressures ? elementData.quadFields.trianglePressures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper cap: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.capFields.pressures ? elementData.capFields.pressures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}
}
