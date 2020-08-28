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

	property real scale: 0.05

	property bool mirrored: false

	property alias item: scene2d.item

	property alias mouseEnabled: scene2d.mouseEnabled

	Entity {
		id: frontEntity

		PlaneMesh {
			id: frontMesh

			width: scene2d.item.width
			height: scene2d.item.height
			mirrored: true
			meshResolution: Qt.size(2, 2)
		}

		Transform {
			id: frontTransform

			rotationX: 90
			scale: root.scale
		}

		components: [frontMesh, frontTransform, material]
	}

	Entity {
		id: backEntity


		PlaneMesh {
			id: backMesh

			width: scene2d.item.width
			height: scene2d.item.height
			mirrored: root.mirrored
			meshResolution: Qt.size(2, 2)
		}

		Transform {
			id: backTransform

			rotationX: 90
			rotationY: 180
			scale: root.mirrored ? root.scale : -root.scale
		}

		components: [backMesh, backTransform, material]
	}

	TransparentTextureMaterial {
		id: material

		texture: texture2d

		Scene2D {
			id: scene2d

			output: RenderTargetOutput {
				attachmentPoint: RenderTargetOutput.Color0
				texture: Texture2D {
					id: texture2d

					width: item.width
					height: item.height
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
		}
	}

	Transform {
		id: transform
	}

	components: [frontEntity, backEntity, transform]
}
