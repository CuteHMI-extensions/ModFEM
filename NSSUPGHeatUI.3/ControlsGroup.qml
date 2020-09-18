import QtQuick 2.0
import QtQuick.Controls 2.0
import QtQuick.Layouts 1.3

import CuteHMI.Services 2.0

import ModFEM.NSSUPGHeat 1.0

GroupBox {
	title: qsTr("Controls")

	property Problem problem

	property bool started: false

	property bool initialized: false

	RowLayout {
		Button {
			text: qsTr("Init")

			onClicked: {
				if (!started) {
					problem.init()
					initialized = true
				}
			}
		}

		Button {
			text: qsTr("Start")

			onClicked: {
				if (initialized) {
					problem.start()
					ServiceManager.start()
					started = true
				}
			}
		}

		Button {
			text: qsTr("Stop")

			onClicked: {
				if (started) {
					problem.stop()
					ServiceManager.stop()
					started = false
				}
			}
		}
	}
}
