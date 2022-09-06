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

#ifndef CGOGN_MODELING_AXIS_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_AXIS_DEFORMATION_TOOL_H_

#include <cgogn/modeling/types/space_deformation_tool.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class AxisDeformationTool : public SpaceDeformationTool<MESH>
{

	using Graph = cgogn::IncidenceGraph;

public:

	Graph* control_axis_; 

	AxisDeformationTool():SpaceDeformationTool<MESH>(), control_axis_vertex_position_(nullptr)
	{
		
	}

	~AxisDeformationTool()
	{

	}

	void create_space_tool(Graph* g, 
		Graph::Attribute<Vec3>* vertex_position, 
		Graph::Attribute<Scalar>* vertex_radius, 
		const std::vector<Vec3>& vertices_positions)
	{
		control_axis_ = g; 
		cgogn::modeling::create_axis(*g, vertex_position, vertex_radius, vertices_positions); 

		control_axis_vertex_position_ = cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		std::shared_ptr<Graph::Attribute<uint32>> position_indices = cgogn::add_attribute<uint32, Graph::Vertex>(*control_axis_, "position_indices");
		//cgogn::modeling::set_graph_attribute_position_indices(*control_axis_, position_indices.get()); 
	}

private : 

	std::shared_ptr<Graph::Attribute<Vec3>> control_axis_vertex_position_; 

}; 

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_H_
