#version 150 core

uniform sampler2D diffuseTexture;

in vec3 position;
in vec2 texCoord;

out vec4 fragColor;

void main()
{
    vec4 texColor = texture(diffuseTexture, texCoord);
    if(texColor.a < 0.75)
	  discard;
    fragColor = texColor;

}
