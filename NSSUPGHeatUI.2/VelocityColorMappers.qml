import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property AbstractColorMapper triangle: HueColorMapper {
		valueBegin: -65.0
		valueEnd: 300
		input: elementData.triangleFields.temperatures ? elementData.triangleFields.temperatures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper quadTriangle: HueColorMapper {
		valueBegin: -65.0
		valueEnd: 300
		input: elementData.quadFields.triangleTemperatures ? elementData.quadFields.triangleTemperatures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}
}
