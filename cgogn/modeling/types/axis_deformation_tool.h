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
#include <algorithm>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class AxisDeformationTool : public SpaceDeformationTool<MESH>
{

	using Graph = cgogn::IncidenceGraph;

public:
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using MeshVertex = typename mesh_traits<MESH>::Vertex;
	using MeshFace = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

	Graph* control_axis_;
	std::shared_ptr<Graph::Attribute<Vec3>> control_axis_vertex_position_;

	AxisDeformationTool() : SpaceDeformationTool<MESH>(), control_axis_vertex_position_(nullptr)
	{
	}

	~AxisDeformationTool()
	{
	}

	void create_space_tool(Graph* g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius,
						   const std::vector<Vec3>& vertex_coords)
	{
		control_axis_ = g;
		axis_skeleton_ = cgogn::modeling::create_axis(*g, vertex_position, vertex_radius, vertex_coords);

		control_axis_vertex_position_ = cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		std::shared_ptr<Graph::Attribute<uint32>> position_indices =
			cgogn::add_attribute<uint32, Graph::Vertex>(*control_axis_, "position_indices");
		// cgogn::modeling::set_graph_attribute_position_indices(*control_axis_, position_indices.get());
	}

	void set_influence_cage_axis(MESH* m, 
								CMap2::Attribute<Vec3>* 				   vertex_position,
								 CMap2::Attribute<Vec3>* local_vertex_position,
								 Graph::Attribute<uint32>* skeleton_vertex, const Vec3& bb_min, const Vec3& bb_max,
								 Vec3& main_direction, Vec3& normal, const double& width)
	{
		this->influence_cage_ = m;
		axis_normal_ = normal;

		Vec3 v = axis_normal_.cross(main_direction);
		v.normalize();

		local_frame_.row(0) = main_direction;
		local_frame_.row(1) = v;
		local_frame_.row(2) = axis_normal_;
		local_frame_inverse_ = local_frame_.inverse();

		const int array_size = 4 * axis_skeleton_.size();

		std::vector<Vec3> vertex_coords(array_size);

		std::vector<Vec3> local_vertex_coords(array_size);

		for (unsigned int v = 0; v < axis_skeleton_.size(); v++)
		{
			const Vec3 position = value<Vec3>(*control_axis_, control_axis_vertex_position_, axis_skeleton_[v]);

			const Vec3 local_bb_min = local_frame_ * (bb_min - position);
			const Vec3 local_bb_max = local_frame_ * (bb_max - position); 

			double n_min = std::min(std::abs(local_bb_min[2]), std::abs(local_bb_max[2])); 

			double n_max = std::max(std::abs(local_bb_min[2]), std::abs(local_bb_max[2])); 

			double front_n, back_n; 
			if (local_bb_min[2] > local_bb_max[2]){
				front_n = n_min; 
				back_n = -n_max; 
			} else {
				back_n = n_min; 
				front_n = -n_max;
			}

			//double front_y, back_y; 

			if (v == 1)
			{
				local_vertex_coords[4 * v] = {0.0, -width, front_n};
				local_vertex_coords[4 * v + 1] = {0.0, width, front_n};

				local_vertex_coords[4 * v + 2] = {0.0, -width, back_n};
				local_vertex_coords[4 * v + 3] = {0.0, width, back_n};
			}
			else
			{
				double min_x, max_x;
				if (local_bb_min[0] < local_bb_max[0])
				{
					min_x = local_bb_min[0];
					max_x = local_bb_max[0];
				}
				else
				{
					max_x = local_bb_min[0];
					min_x = local_bb_max[0];
				}

				if (v == 0)
				{
					local_vertex_coords[4 * v] = {max_x, -width, front_n};
					local_vertex_coords[4 * v + 1] = {max_x, width, front_n};

					local_vertex_coords[4 * v + 2] = {max_x, -width, back_n};
					local_vertex_coords[4 * v + 3] = {max_x, width, back_n};
				}
				else
				{
					local_vertex_coords[4 * v] = {min_x, -width, front_n};
					local_vertex_coords[4 * v + 1] = {min_x, width, front_n};

					local_vertex_coords[4 * v + 2] = {min_x, -width, back_n};
					local_vertex_coords[4 * v + 3] = {min_x, width, back_n};
				}
			}

			vertex_coords[4 * v] = (local_frame_inverse_ * local_vertex_coords[4 * v]) + position;
			vertex_coords[4 * v + 1] = (local_frame_inverse_ * local_vertex_coords[4 * v + 1]) + position;

			vertex_coords[4 * v + 2] = (local_frame_inverse_ * local_vertex_coords[4 * v + 2]) + position;
			vertex_coords[4 * v + 3] = (local_frame_inverse_ * local_vertex_coords[4 * v + 3]) + position;
		}

		cgogn::modeling::create_axis_box(*m, vertex_position, local_vertex_position, skeleton_vertex, vertex_coords,
										 local_vertex_coords);

		this->influence_cage_vertex_position_ =
			cgogn::get_attribute<Vec3, MeshVertex>(*m, "position");

		std::shared_ptr<Attribute<uint32>> position_indices =
			cgogn::add_attribute<uint32, MeshVertex>(*(this->influence_cage_), "position_indices");
		cgogn::modeling::set_attribute_position_indices(*(this->influence_cage_),
														position_indices.get());

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, MeshVertex>(*(this->influence_cage_), "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*(this->influence_cage_),
													   marked_vertices.get());
	}

private:
	std::vector<Graph::Vertex> axis_skeleton_;
	Vec3 axis_normal_;
	Eigen::Matrix3d local_frame_;
	Eigen::Matrix3d local_frame_inverse_;
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_H_
