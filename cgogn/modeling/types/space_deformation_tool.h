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

#ifndef CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cells_set.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class SpaceDeformationTool
{
public:
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

	std::vector<Vertex> object_influence_area_; 

	SpaceDeformationTool()
	{
	}

	virtual ~SpaceDeformationTool()
	{
	}

	virtual void create_space_tool()
	{
	}

	virtual void update_deformation_object()
	{
	}

	void set_deformation_type(const std::string new_type){
		deformation_type_ = new_type; 
	}

protected:
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> position_weights_;

	

	std::string deformation_type_;

private:
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_
