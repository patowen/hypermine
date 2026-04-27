#version 450

#include "common.h"

layout(location = 0) in vec4 position;
layout(location = 1) in vec3 texcoords;
layout(location = 2) in vec4 normal;

layout(location = 0) out vec3 texcoords_out;
layout(location = 1) out vec4 normal_out;

layout(push_constant) uniform PushConstants {
    mat4 transform;
};

void main() {
    gl_Position = view_projection * transform * position;
    texcoords_out = texcoords;
    normal_out = transform * normal;
}
