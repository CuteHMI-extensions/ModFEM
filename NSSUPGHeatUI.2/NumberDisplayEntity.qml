import CuteHMI.GUI 1.0

ItemEntity {
	property alias display: display

	mirrored: true

	item: NumberDisplay {
		id: display
	}
}
