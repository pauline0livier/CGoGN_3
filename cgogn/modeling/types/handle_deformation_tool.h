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

#ifndef CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_

#include <cgogn/modeling/types/space_deformation_tool.h>
#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH, typename GRAPH>
class HandleDeformationTool : public SpaceDeformationTool<MESH>
{
	template <typename T>
	using Attribute = typename mesh_traits<GRAPH>::template Attribute<T>;
	using Vertex = typename mesh_traits<GRAPH>::Vertex;
	using Edge = typename mesh_traits<GRAPH>::Face;

	 using Vec2 = geometry::Vec2;
    using Vec3 = geometry::Vec3;

public:
	HandleDeformationTool()
	{
		handle_vertex_position_ = nullptr; 
	}

	~HandleDeformationTool()
		{
		}

	void create_space_tool()
	{
		std::cout << "handle" << std::endl; 
	}

private :
	std::shared_ptr<Attribute<Vec3>> handle_vertex_position_;

}; 

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_
