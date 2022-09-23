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

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class HandleDeformationTool : public SpaceDeformationTool<MESH>
{
	using Graph = cgogn::IncidenceGraph;

public:
template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using MeshVertex = typename mesh_traits<MESH>::Vertex;
	using MeshFace = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

	Graph* control_handle_;

	std::shared_ptr<Graph::Attribute<Vec3>> control_handle_vertex_position_;
	std::shared_ptr<boost::synapse::connection> handle_attribute_update_connection_;

	HandleDeformationTool() : SpaceDeformationTool<MESH>(), control_handle_vertex_position_(nullptr)
	{
	}

	~HandleDeformationTool()
	{
	}

	void create_space_tool(Graph* g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius,
						   const Vec3& center1, const Vec3& center2)
	{
		control_handle_ = g;
		handle_vertex_ = cgogn::modeling::create_handle(*g, vertex_position, vertex_radius, center1, center2);
 
		control_handle_vertex_position_ = cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		/*std::shared_ptr<Graph::Attribute<uint32>> position_indices =
			cgogn::add_attribute<uint32, Graph::Vertex>(*control_handle_, "position_indices");*/
		//cgogn::modeling::set_graph_attribute_position_indices(*control_handle_, position_indices.get());
	}

	void set_influence_cage_handle(MESH* m, CMap2::Attribute<Vec3>* vertex_position,
								   CMap2::Attribute<Vec3>* local_vertex_position, const Vec3& bb_min,
								   const Vec3& bb_max, const Vec3& handle_position, Vec3& normal)
	{
		SpaceDeformationTool<MESH>::influence_cage_ = m;
		handle_normal_ = normal;
		handle_position_ = handle_position; 

		// plane ax+by+cz+d = 0, with (a,b,c) = handle_normal and (x,y,z) = bb_min

		// bbmin belongs to the plane
		float d_handle_min = -(handle_normal_.dot(bb_min));

		// find alpha such that center + alpha*handle_normal belongs to plane
		float alpha_min = -d_handle_min - (handle_normal_.dot(handle_position));

		Eigen::Vector3d center_min_plane = handle_position + alpha_min * handle_normal_;

		Vec3 u = bb_min - center_min_plane;
		u.normalize();

		Vec3 v = handle_normal_.cross(u);
		v.normalize();

		local_frame_.row(0) = u;
		local_frame_.row(1) = v;
		local_frame_.row(2) = handle_normal_;
		local_frame_inverse_ = local_frame_.inverse();

		const Vec3 local_bb_min = local_frame_ * (bb_min - handle_position_);
		const Vec3 local_bb_max = local_frame_ * (bb_max - handle_position_);

		const double radius = local_bb_min.norm(); 

		Vec3 min_depth, max_depth; 
		if (local_bb_min[2] > local_bb_max[2]){
			cgogn::modeling::create_handle_box(*m, vertex_position, local_vertex_position, handle_position, radius, local_frame_inverse_, local_bb_min[2], local_bb_max[2]);
		} else {
			cgogn::modeling::create_handle_box(*m, vertex_position, local_vertex_position, handle_position, radius, local_frame_inverse_, local_bb_max[2], local_bb_min[2]);
		}
		

		SpaceDeformationTool<MESH>::influence_cage_vertex_position_ = cgogn::get_attribute<Vec3, MeshVertex>(*m, "position");

		std::shared_ptr<Attribute<uint32>> position_indices =
			cgogn::add_attribute<uint32, MeshVertex>(*SpaceDeformationTool<MESH>::influence_cage_, "position_indices");
		cgogn::modeling::set_attribute_position_indices(*SpaceDeformationTool<MESH>::influence_cage_, position_indices.get());

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, MeshVertex>(*SpaceDeformationTool<MESH>::influence_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*SpaceDeformationTool<MESH>::influence_cage_, marked_vertices.get());
	}

	void update_influence_cage_position()
	{
		handle_position_ = value<Vec3>(*control_handle_, control_handle_vertex_position_, handle_vertex_);

		std::shared_ptr<Attribute<Vec3>> local_position =
			cgogn::get_attribute<Vec3, MeshVertex>(*SpaceDeformationTool<MESH>::influence_cage_, "local_position");

		foreach_cell(*SpaceDeformationTool<MESH>::influence_cage_, [&](MeshVertex v) -> bool {
			Vec3 local_point = value<Vec3>(*SpaceDeformationTool<MESH>::influence_cage_, local_position, v); 

			value<Vec3>(*SpaceDeformationTool<MESH>::influence_cage_, SpaceDeformationTool<MESH>::influence_cage_vertex_position_, v) = local_frame_inverse_*local_point + handle_position_;

			return true;
		});
	}

	void set_up_attenuation(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{

		uint32 nbv_object = nb_cells<MeshVertex>(object);

		SpaceDeformationTool<MESH>::attenuation_.resize(nbv_object);
		SpaceDeformationTool<MESH>::attenuation_.setZero();

		compute_attenuation_cage(object);

	}

private:
	Graph::Vertex handle_vertex_;
	Vec3 handle_position_; // also frame_origin 
	Vec3 handle_normal_;
	Eigen::Matrix3d local_frame_;
	Eigen::Matrix3d local_frame_inverse_;

	void compute_attenuation_cage(MESH& object)
	{
		std::shared_ptr<Attribute<uint32>> cage_face_indices =
			add_attribute<uint32, MeshFace>(*SpaceDeformationTool<MESH>::influence_cage_, "face_indices");

		cgogn::modeling::set_attribute_face_indices(*SpaceDeformationTool<MESH>::influence_cage_,
													cage_face_indices.get());

		std::shared_ptr<Attribute<uint32>> object_position_indices =
			get_attribute<uint32, MeshVertex>(object, "position_indices");

		std::shared_ptr<Attribute<uint32>> i_cage_position_indices =
			get_attribute<uint32, MeshVertex>(*SpaceDeformationTool<MESH>::influence_cage_, "position_indices");

		uint32 nbf_cage = 2 * nb_cells<MeshFace>(*SpaceDeformationTool<MESH>::influence_cage_);
		uint32 nbv_cage = nb_cells<MeshVertex>(*SpaceDeformationTool<MESH>::influence_cage_);

		float h = 0.0f;
		float max_dist = 0.0f; 
		std::vector<Vec2> attenuation_points;
		SpaceDeformationTool<MESH>::influence_area_->foreach_cell([&](MeshVertex v) {
			uint32 surface_point_idx = value<uint32>(object, object_position_indices, v);

			float i_dist = SpaceDeformationTool<MESH>::cage_influence_distance(surface_point_idx, nbf_cage, nbv_cage);

			SpaceDeformationTool<MESH>::attenuation_(surface_point_idx) = (float)sin(0.5*M_PI * (i_dist ));
		
		});

	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_
