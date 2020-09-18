import QtQuick 2.0

import ModFEM.NSSUPGHeat 1.0

QtObject {
	property ElementData elementData

	property var bcPalette: [
		"pink",
		"red",
		"yellow",
		"blue",
		"green",
		"cyan",
		"silver",
		"orange",
		"purple",

//		Qt.hsla(0.0, 1.0, 0.35, 1.0),
//		Qt.hsla(0.125, 1.0, 0.35, 1.0),
//		Qt.hsla(0.25, 1.0, 0.35, 1.0),
//		Qt.hsla(0.5, 1.0, 0.35, 1.0),
//		Qt.hsla(0.625, 1.0, 0.35, 1.0),
//		Qt.hsla(0.875, 1.0, 0.35, 1.0),

		Qt.hsla(0.0, 0.25, 0.75, 1.0),
		Qt.hsla(0.125, 0.25, 0.75, 1.0),
		Qt.hsla(0.25, 0.25, 0.75, 1.0),
		Qt.hsla(0.5, 0.25, 0.75, 1.0),
		Qt.hsla(0.625, 0.25, 0.75, 1.0),
		Qt.hsla(0.875, 0.25, 0.75, 1.0),
	]

	property AbstractColorMapper triangle: PaletteColorMapper {
		palette: bcPalette
		input: elementData.triangles.boundaryConditions ? elementData.triangles.boundaryConditions : new ArrayBuffer
	}

	property AbstractColorMapper quadTriangle: PaletteColorMapper {
		palette: bcPalette
		input: elementData.quads.triangleBoundaryConditions ? elementData.quads.triangleBoundaryConditions : new ArrayBuffer
	}

	property AbstractColorMapper cap: PaletteColorMapper {
		palette: bcPalette
		input: new ArrayBuffer
	}
}
