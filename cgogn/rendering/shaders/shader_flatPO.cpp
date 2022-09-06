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

#include <cgogn/rendering/shaders/shader_flatPO.h>

namespace cgogn
{

namespace rendering
{

ShaderFlatPO* ShaderFlatPO::instance_ = nullptr;

ShaderFlatPO::ShaderFlatPO()
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
		//uniform vec4 ambient_color;
		//uniform vec4 diffuse_color;
		//uniform vec4 specular_color;
		//uniform float shininess;
		//uniform float alpha;
		uniform vec3 light_position;
		uniform bool ghost_mode;
		
		in vec3 position;
		in vec2 tc; 

		out vec4 frag_out;

		

		void main()
		{
			vec3 N = normalize(cross(dFdx(position), dFdy(position)));

			vec3 LightDir = light_position - position; 
			float distance = length(LightDir); 
			distance = distance*distance; 
			LightDir = normalize(LightDir);
			float lambert = max(dot(N, LightDir), 0.0);
			if (ghost_mode)
				lambert = 0.4 * pow(1.0 - lambert, 2);
			
			float specular = 0.0;
			//if (lambert > 0.0){
				//vec3 viewDir = normalize(-position);

				//vec3 halfDir = normalize(LightDir + viewDir);
    			//float specAngle = max(dot(halfDir, N), 0.0);
    			//specular = pow(specAngle, shininess);

      			//vec3 reflectDir = reflect(-LightDir, N);
      			//float specAngle = max(dot(reflectDir, viewDir), 0.0);
      			// note that the exponent is different here
      			//specular = pow(specAngle, shininess/4.0);
    
  			//}		

			//vec3 newColor = normalize(vec3(abs(position.y), abs(position.z), abs(position.x))); 
			//vec3 colorLinear = ambient_color.xyz + diffuse_color.xyz*lambert* vec3(1.0, 1.0, 1.0) * 40.0 / distance +
                    //specular_color.xyz * specular * vec3(1.0, 1.0, 1.0) * 10.0 / distance; 


			//vec3 colorGammaCorrected = pow(colorLinear, vec3(1.0 / 2.2));
			
			frag_out = vec4(abs(LightDir.x-N.x), abs(LightDir.y-N.y), abs(LightDir.z - N.z), 1.0); 
			//frag_out = vec4(colorLinear, alpha); 
		}
	)";

	load2_bind(vertex_shader_source, fragment_shader_source, "vertex_position","vertex_tc");
	//get_uniforms();
}

void ShaderParamFlatPO::set_uniforms()
{
	//shader_->set_uniforms_values();
		/*ambient_color_, diffuse_color_, specular_color_, shininess_, alpha_, light_position_,
								 ghost_mode_*/ 
}

} // namespace rendering

} // namespace cgogn
