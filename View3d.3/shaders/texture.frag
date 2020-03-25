#version 330 core

in vec3 ourColor;
in vec2 TexCoord;

out vec4 color;

uniform sampler2D ourTexture;

void main()
{
    vec4 fragColor = texture(ourTexture, TexCoord) * vec4(ourColor, 1.0f);
    if(fragColor.a < 0.7)
	  discard;
    color = fragColor;
}
