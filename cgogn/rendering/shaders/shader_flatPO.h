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

#ifndef CGOGN_RENDERING_SHADERS_FLAT_PO_H_
#define CGOGN_RENDERING_SHADERS_FLAT_PO_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(FlatPO, false, CGOGN_STR(FlatPO))

class CGOGN_RENDERING_EXPORT ShaderParamFlatPO : public ShaderParam
{

	void set_uniforms() override;

public:
	GLColor ambient_color_;
	GLColor diffuse_color_;
	GLColor specular_color_;
	GLfloat shininess_;
	GLfloat alpha_;  
	GLVec3 light_position_;
	bool double_side_;
	bool ghost_mode_;

	using ShaderType = ShaderFlatPO;

	ShaderParamFlatPO(ShaderType* sh)
		: ShaderParam(sh), ambient_color_(0.05f, 0.05f, 0.05f, 1), diffuse_color_(0.05f, 0.05f, 0.05f, 1), specular_color_(0.05f, 0.05f, 0.05f, 1),
		shininess_(40.0), alpha_(1.0), light_position_(10, 100, 1000), ghost_mode_(false)
	{
	}

	inline ~ShaderParamFlatPO() override
	{
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_FLAT_PO_H_
