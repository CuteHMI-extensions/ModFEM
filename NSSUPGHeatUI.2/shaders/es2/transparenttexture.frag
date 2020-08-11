#version 150 core
#define FP highp

uniform sampler2D diffuseTexture;

varying FP vec3 position;
varying FP vec2 texCoord;

void main()
{
    FP vec4 texColor = texture2D(diffuseTexture, texCoord);
    if(texColor.a < 0.75)
	  discard;
    gl_FragColor = texColor;
}
