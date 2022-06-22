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

#include <cgogn/rendering/shaders/shader_flat.h>

namespace cgogn
{

namespace rendering
{

ShaderFlat* ShaderFlat::instance_ = nullptr;

ShaderFlat::ShaderFlat()
{
	const char* vertex_shader_source = R"(
		#version 150
		uniform mat4 projection_matrix;
		uniform mat4 model_view_matrix;

		in vec3 vertex_position;
		
		out vec3 position;
		
		void main()
		{
			vec4 position4 = model_view_matrix * vec4(vertex_position, 1.0);
			position = position4.xyz;
			gl_Position = projection_matrix * position4;
		}
	)";

	const char* fragment_shader_source = R"(
		#version 150
		uniform vec4 front_color;
		uniform vec4 back_color;
		uniform vec4 ambiant_color;
		uniform vec3 light_position;
		uniform bool double_side;
		uniform bool ghost_mode;
		
		in vec3 position;

		out vec4 frag_out;

		void main()
		{
			vec3 N = normalize(cross(dFdx(position), dFdy(position)));
			vec3 L = normalize(light_position - position);
			float lambert = dot(N, L);
			if (ghost_mode)
				lambert = 0.4 * pow(1.0 - lambert, 2);
			if (gl_FrontFacing)
				//frag_out = vec4(ambiant_color.rgb + lambert * front_color.rgb, front_color.a);
				frag_out = vec4(abs(L.x-N.x), abs(L.y-N.y), abs(L.z - N.z), 1.0); 
			else
				if (!double_side)
					discard;
				else frag_out = vec4(ambiant_color.rgb + lambert * back_color.rgb, back_color.a);
		}
	)";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_position");
	get_uniforms("front_color", "back_color", "ambiant_color", "light_position", "double_side", "ghost_mode");
}

void ShaderParamFlat::set_uniforms()
{
	shader_->set_uniforms_values(front_color_, back_color_, ambiant_color_, light_position_, double_side_, ghost_mode_);
}

} // namespace rendering

} // namespace cgogn
