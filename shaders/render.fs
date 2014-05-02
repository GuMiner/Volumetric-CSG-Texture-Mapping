#version 430 core 

uniform sampler2D sampTex;

out vec4 color;

in VS_OUT
{
    vec4 color;
    float dStat;
    vec2 texCoord;
} fs_in;

void main(void)
{
    if (fs_in.dStat < 0.75)
    {
        color = fs_in.color;
    }
    else
    {
        ivec2 size = textureSize(sampTex, 0);
        color = texelFetch(sampTex, ivec2(size.x*fs_in.texCoord.x, size.y*fs_in.texCoord.y), 0);
    }
}