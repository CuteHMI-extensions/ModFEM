#version 330 core

in vec3 vertexPosition;
in vec3 vertexColor;
in vec2 vertexTexCoord;

out vec3 ourColor;
out vec2 TexCoord;

void main()
{
    gl_Position = vec4(vertexPosition, 1.0f);
    ourColor = vertexColor;
    TexCoord = vertexTexCoord;
}
