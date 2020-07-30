#version 150 core

out vec4 fragColor;

in FragData {
    vec3 fragPosition;
    vec3 fragNormal;
} fd;

void main()
{
    fragColor = vec4(0.0, 1.0, 0.0, 1.0) + (1.0 * vec4(fd.fragPosition, 1.0));
}
