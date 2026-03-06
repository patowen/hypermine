#version 450

layout(location = 0) in vec3 texcoords;
layout(location = 1) in vec4 normal;

layout(location = 0) out vec4 color_out;

layout(set = 1, binding = 0) uniform sampler2DArray textures;

void main() {
    color_out = texture(textures, texcoords);
}
