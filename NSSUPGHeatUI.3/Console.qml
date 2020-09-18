import QtQml 2.0

import CuteHMI.Services 2.0
import CuteHMI.SharedDatabase 0.0
import CuteHMI.DataAcquisition 0.0

QtObject {
	objectName: "console"

	property Project project: Project {
		db.threaded: false
	}
}
