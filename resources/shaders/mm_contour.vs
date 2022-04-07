#version 110

void main()
{
    // Add small epsilon to z to solve z-fighting between painted triangles and contour lines.
	vec4 clip_position = gl_ModelViewProjectionMatrix * gl_Vertex;
	clip_position.z -= 0.00001 * abs(clip_position.w);
    gl_Position = clip_position;
}
