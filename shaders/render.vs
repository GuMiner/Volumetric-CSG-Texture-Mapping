#version 430 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec2 texCoord;

out VS_OUT
{
    vec4 color;
    float dStat;
    vec2 texCoord;
} vs_out;

uniform mat4 mv_matrix;
uniform mat4 proj_matrix;
uniform float dStatus;

void main(void)
{
    vec4 pos = vec4(position, 1);

    // Compute where we are for this shape
    int count = 10;
    float spacing = 1.5;
    int y = gl_InstanceID/count;
    int x = gl_InstanceID - y*count;
    gl_Position = proj_matrix * mv_matrix * (pos + x*vec4(spacing, 0, 0, 0.0) + y*vec4(0, spacing, 0, 0.0));
    
    // Output stuff to the fragment shader
    vs_out.dStat = dStatus;
    
    if (dStatus > 0.75)
    {
        // Draw the triangles with texture mapping
        vs_out.color = vec4(1.0, 0.0, 0.0, 0.0);
        vs_out.texCoord = texCoord;
    }
    else
    {
        // Draw the triangles with fragment-shader mapping
        vs_out.texCoord = vec2(0, 0);
        vs_out.color = vec4(pos.x/1, pos.y/1,pos.z/1, 0);
    }   
}