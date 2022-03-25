#version 150

// see as reference: https://github.com/mhalber/Lines/blob/master/geometry_shader_lines.h

uniform mat4 view_model_matrix;
uniform mat4 projection_matrix;

in vec3 v_position;

void main()
{
    gl_Position = projection_matrix * view_model_matrix * vec4(v_position, 1.0);
}