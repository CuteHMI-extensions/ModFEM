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

//#version 330 core

//in vec3 ourColor;
//in vec2 texCoord;

//out vec4 color;

//uniform sampler2D diffuseTexture;

//void main()
//{
//    vec4 fragColor = texture(diffuseTexture, texCoord);
//    if(fragColor.a < 0.7)
//	  discard;
//    color = fragColor;
//}
