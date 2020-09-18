import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property real minimum: -65.0

	property real maximum: 300.0

	property AbstractColorMapper triangle: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.triangleFields.temperatures ? elementData.triangleFields.temperatures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper quadTriangle: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.quadFields.triangleTemperatures ? elementData.quadFields.triangleTemperatures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}

	property AbstractColorMapper cap: HueColorMapper {
		valueBegin: minimum
		valueEnd: maximum
		input: elementData.capFields.temperatures ? elementData.capFields.temperatures : new ArrayBuffer
		inputType: HueColorMapper.INPUT_DOUBLE
	}
}
