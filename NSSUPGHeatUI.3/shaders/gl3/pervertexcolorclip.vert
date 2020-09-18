#version 150 core

in vec3 vertexPosition;
in vec3 vertexNormal;
in vec4 vertexColor;

out vec3 worldPosition;
out vec3 worldNormal;
out vec4 color;

uniform mat4 modelMatrix;
uniform mat3 modelNormalMatrix;
uniform mat4 mvp;

struct Plane {
    vec4 equation;
};

struct ClipPlanesData {
    int count;
    Plane planes[8];
};

uniform ClipPlanesData clipPlanesData;

void main()
{
    worldNormal = normalize(modelNormalMatrix * vertexNormal);
    worldPosition = vec3(modelMatrix * vec4( vertexPosition, 1.0 ));
    color = vertexColor;

    // <ModFEM.NSSUPGHeatUI-1.workaround target="nvidia" cause="bug">
    gl_ClipDistance[0] = 1.0;
    // </ModFEM.NSSUPGHeatUI-1.workaround>

    for (int i = 0; i < clipPlanesData.count; ++i) {
	gl_ClipDistance[i] = dot(vec4(worldPosition, 1.0), clipPlanesData.planes[i].equation);
    }

    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
