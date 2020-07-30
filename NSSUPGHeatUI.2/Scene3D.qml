import QtQuick 2.0 as QQ2
import QtQuick.Scene3D 2.0

import Qt3D.Core 2.5
import Qt3D.Render 2.14
import Qt3D.Input 2.5
import Qt3D.Extras 2.5

import ModFEM.NSSUPGHeat 1.0

Scene3D {
	id: scene3d

	aspects: ["input", "logic"]
	cameraAspectRatioMode: Scene3D.AutomaticAspectRatio
}
