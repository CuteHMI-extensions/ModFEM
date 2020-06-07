import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property var bcPalette: [
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
//		Qt.rgba(0.8, 0.1, 0.4, 1.0),
		Qt.hsla(0.0, 1.0, 0.35, 1.0),
		Qt.hsla(0.125, 1.0, 0.35, 1.0),
		Qt.hsla(0.25, 1.0, 0.35, 1.0),
		Qt.hsla(0.375, 1.0, 0.35, 1.0),
		Qt.hsla(0.5, 1.0, 0.35, 1.0),
		Qt.hsla(0.625, 1.0, 0.35, 1.0),
		Qt.hsla(0.75, 1.0, 0.35, 1.0),
		Qt.hsla(0.875, 1.0, 0.35, 1.0)
	]

	property AbstractColorMapper triangle: PaletteColorMapper {
		palette: bcPalette
		input: elementData.triangles.boundaryConditions ? elementData.triangles.boundaryConditions : new ArrayBuffer
	}

	property AbstractColorMapper quadTriangle: PaletteColorMapper {
		palette: bcPalette
		input: elementData.quads.triangleBoundaryConditions ? elementData.quads.triangleBoundaryConditions : new ArrayBuffer
	}
}
