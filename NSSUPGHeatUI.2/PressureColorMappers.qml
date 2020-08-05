import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property AbstractColorMapper triangle: HueColorMapper {
		valueBegin: -65.0
		valueEnd: 300
		input: elementData.triangleFields.pressures ? elementData.triangleFields.pressures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper quadTriangle: HueColorMapper {
		valueBegin: -65.0
		valueEnd: 300
		input: elementData.quadFields.trianglePressures ? elementData.quadFields.trianglePressures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}
}
