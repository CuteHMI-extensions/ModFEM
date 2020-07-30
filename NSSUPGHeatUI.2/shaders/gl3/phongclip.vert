#version 150 core

in vec3 vertexPosition;
in vec3 vertexNormal;

out VertexData {
    vec3 position;
    vec3 normal;
} vd;

struct Plane {
    vec4 equation;
};

struct ClipPlanesData {
    int count;
    Plane planes[8];
};

uniform ClipPlanesData clipPlanesData;
uniform mat4 modelView;
uniform mat4 modelMatrix;
uniform mat3 modelViewNormal;
uniform mat4 mvp;

void main()
{

    vd.normal = normalize( modelViewNormal * vertexNormal );
    vd.position = vec4(modelView * vec4( vertexPosition, 1.0 )).xyz;
    vec3 worldPos = vec4(modelMatrix * vec4(vertexPosition, 1.0)).xyz;

    for (int i = 0; i < clipPlanesData.count; ++i) {
	gl_ClipDistance[i] = dot(vec4(worldPos, 1.0), clipPlanesData.planes[i].equation);
    }

    gl_Position = mvp * vec4( vertexPosition, 1.0 );
}
