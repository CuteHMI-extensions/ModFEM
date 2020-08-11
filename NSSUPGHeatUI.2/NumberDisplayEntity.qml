import QtQuick.Scene2D 2.9
import Qt3D.Render 2.13
import Qt3D.Core 2.0
import Qt3D.Render 2.14
import Qt3D.Input 2.0
import Qt3D.Extras 2.14

import CuteHMI.GUI 1.0

Entity {
	id: root

	property alias transform: transform

	property alias display: display

	PlaneMesh {
		id: mesh

		width: display.width
		height: display.height
		mirrored: true
		meshResolution: Qt.size(2, 2)
	}

	Entity {
		id: frontEntity

		Transform {
			id: frontTransform

			rotationX: 90
			scale: 0.05
		}

		components: [mesh, frontTransform, material]
	}

	Entity {
		id: backEntity

		Transform {
			id: backTransform

			rotationX: 90
			rotationY: 180
			scale: 0.05
		}

		components: [mesh, backTransform, material]
	}

	TransparentTextureMaterial {
		id: material

		texture: texture2d

		Scene2D {
			output: RenderTargetOutput {
				attachmentPoint: RenderTargetOutput.Color0
				texture: Texture2D {
					id: texture2d

					width: display.width
					height: display.height
					format: Texture.RGBA8_UNorm
					generateMipMaps: true
					magnificationFilter: Texture.Linear
					minificationFilter: Texture.LinearMipMapLinear
					wrapMode {
						x: WrapMode.ClampToEdge
						y: WrapMode.ClampToEdge
					}
				}
			}

			mouseEnabled: false

			NumberDisplay {
				id: display

				active: true
			}
		}
	}

	Transform {
		id: transform
	}

	components: [frontEntity, backEntity, transform]
}
