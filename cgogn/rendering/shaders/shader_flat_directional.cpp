/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/rendering/shaders/shader_flat_directional.h>

namespace cgogn
{

namespace rendering
{

ShaderFlatDirectional* ShaderFlatDirectional::instance_ = nullptr;

ShaderFlatDirectional::ShaderFlatDirectional()
{
	const char* vertex_shader_source = R"(
		#version 150
		
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;
	
		in vec3 vertex_position;
		in vec2 vertex_tc;
		
		out vec3 position;
		out vec2 tc;
		
		void main()
		{
			tc = vertex_tc;
			vec4 position4 = model_view_matrix * vec4(vertex_position, 1.0);
			position = position4.xyz;
			gl_Position = projection_matrix * position4;
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform vec3 light_position;
		uniform bool ghost_mode;
		
		in vec3 position;
		in vec2 tc; 

		out vec4 frag_out;

		void main()
		{
			vec3 N = normalize(cross(dFdx(position), dFdy(position)));

			vec3 LightDirection = light_position - position; 
			float distance = length(LightDirection); 
			distance = distance*distance; 
			LightDirection = normalize(LightDirection);
			
			frag_out = vec4(abs(LightDirection.x-N.x), abs(LightDirection.y-N.y), 
									abs(LightDirection.z - N.z), 1.0); 
		}
	)";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_position","vertex_tc");
	
}

void ShaderParamFlatDirectional::set_uniforms()
{

}

} // namespace rendering

} // namespace cgogn
