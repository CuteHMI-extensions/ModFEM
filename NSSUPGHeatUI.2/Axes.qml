import QtQuick 2.0

/**
  Overlay axes item.
  */
Canvas {
	id: root

	implicitWidth: 100
	implicitHeight: 100

	property matrix4x4 rotationMatrix: Qt.matrix4x4()

	property real arrowLength: 0.4 * width

	property real lineWidth: 2

	property real headRatio: 0.25

	property real headAngle: 22.5

	property color xColor: "red"

	property color yColor: "green"

	property color zColor: "blue"

	property color labelColor: "white"

	property vector3d implicitXVector: Qt.vector3d(arrowLength, 0, 0)

	property vector3d implicitYVector: Qt.vector3d(0, arrowLength, 0)

	property vector3d implicitZVector: Qt.vector3d(0, 0, arrowLength)

	property vector3d xVector: implicitXVector

	property vector3d yVector: implicitYVector

	property vector3d zVector: implicitZVector

	onPaint: {
		var ctx = getContext("2d");
		ctx.save()
		ctx.reset()

		var vectorArr = [{vector: xVector, color: xColor, label: "X"},
						 {vector: yVector, color: yColor, label: "Y"},
						 {vector: zVector, color: zColor, label: "Z"}]
		vectorArr.sort((a, b) => a.vector.z - b.vector.z);

		for (var i = 0; i < vectorArr.length; i++)
			drawArrow(ctx, vectorArr[i].vector, vectorArr[i].color, vectorArr[i].label)

		ctx.restore()
	}

	onRotationMatrixChanged: {
		xVector = implicitXVector.times(rotationMatrix)
		yVector = implicitYVector.times(rotationMatrix)
		zVector = implicitZVector.times(rotationMatrix)
		requestPaint()
	}

	function rotateVector2d(vector, angle)
	{
		var cosAngle = Math.cos(angle)
		var sinAngle = Math.sin(angle)
		return Qt.vector2d(vector.x * cosAngle - vector.y * sinAngle, vector.x * sinAngle + vector.y * cosAngle)
	}

	function drawArrow(ctx, vector, color, label)
	{
		var center = root.width / 2
		var headSize = root.arrowLength * root.headRatio

		ctx.lineWidth = lineWidth
		ctx.strokeStyle = Qt.darker(color, Math.max(1.0, 1.0 - vector.z / root.arrowLength))
		ctx.fillStyle = Qt.darker(color, Math.max(1.0, 1.0 - vector.z / root.arrowLength))
		ctx.beginPath()
		ctx.moveTo(center, center)
		ctx.lineTo(center + vector.x, center - vector.y)
		ctx.stroke()
		ctx.beginPath()
		ctx.moveTo(center + vector.x, center - vector.y)

		var projectionVector = Qt.vector2d(vector.x, vector.y)
		var projectionRatio = projectionVector.length() / vector.length()
		var normalizedProjectionVector = projectionVector.normalized()
		var angle = root.headAngle * Math.PI /180
		var cosAngle = Math.cos(angle)
		var headVector = rotateVector2d(vector, angle).normalized().times(headSize / cosAngle)
		var offsetVector = normalizedProjectionVector.times(headSize * (1 - projectionRatio))
		ctx.lineTo(center + vector.x - headVector.x + offsetVector.x, center - vector.y + headVector.y - offsetVector.y)
		headVector = rotateVector2d(vector, -angle).normalized().times(headSize / cosAngle)
		ctx.lineTo(center + vector.x - headVector.x + offsetVector.x, center - vector.y + headVector.y - offsetVector.y)
		ctx.closePath()
		ctx.fill()

		ctx.fillStyle = Qt.darker(labelColor, Math.max(1.0, 1.0 - vector.z / root.arrowLength))
		ctx.textAlign = "center"
		ctx.textBaseline = "middle"
		ctx.font = headSize.toString() + "px monospace"
		offsetVector = normalizedProjectionVector.times(0.5 * headSize * projectionRatio)
		ctx.fillText(label, center + vector.x + offsetVector.x, center - vector.y - offsetVector.y)
	}
}
