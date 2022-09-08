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

#ifndef CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/modeling/types/space_deformation_tool.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class CageDeformationTool : public SpaceDeformationTool<MESH>
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

public:
	MESH* control_cage_;
	std::shared_ptr<Attribute<Vec3>> control_cage_vertex_position_;

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;

	CageDeformationTool() : SpaceDeformationTool<MESH>(), m_hFactor(-1.0f), control_cage_vertex_position_(nullptr)
	{
	}

	~CageDeformationTool()
	{
	}

	void create_space_tool(MESH* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max)
	{
		control_cage_ = m;
		cgogn::modeling::create_box(*m, vertex_position, bb_min, bb_max);

		control_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		std::shared_ptr<Attribute<uint32>> position_indices =
			cgogn::add_attribute<uint32, Vertex>(*control_cage_, "position_indices");
		cgogn::modeling::set_attribute_position_indices(*control_cage_, position_indices.get());

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, Vertex>(*control_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*control_cage_, marked_vertices.get());
	}

	void set_center_control_cage(Vec3& center)
	{
		center_control_cage_ = center;
	}

	void update_influence_cage_position()
	{
		foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, v);

			value<Vec3>(*SpaceDeformationTool<MESH>::influence_cage_,
						SpaceDeformationTool<MESH>::influence_cage_vertex_position_, v) =
				((cage_point - center_control_cage_) * 1.5) + center_control_cage_;

			return true;
		});
	}

	void set_up_attenuation(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
 
		std::shared_ptr<Attribute<uint32>> object_position_indices =
			get_attribute<uint32, Vertex>(object, "position_indices");

		uint32 nbv_object = nb_cells<Vertex>(object);

		control_area_validity_.resize(nbv_object);
		control_area_validity_.setZero();

		foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, vertex_position, v);

			bool inside_cage = local_mvc_pt_control_area(surface_point); 

			if (inside_cage)
			{
				uint32 surface_point_idx = value<uint32>(object, object_position_indices, v);
				control_area_validity_(surface_point_idx) = 1.0f;
			}

			return true;
		});

		SpaceDeformationTool<MESH>::attenuation_.resize(nbv_object);
		SpaceDeformationTool<MESH>::attenuation_.setZero();

		compute_attenuation_cage(object);

	}

private:
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> control_cage_coords_;

	Eigen::VectorXd mixing_factor_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_comp_;

	float m_hFactor;

	Vec3 center_control_cage_;

	Eigen::VectorXd control_area_validity_;

	void compute_attenuation_cage(MESH& object)
	{
		 
		std::shared_ptr<Attribute<uint32>> cage_face_indices =
			add_attribute<uint32, Face>(*SpaceDeformationTool<MESH>::influence_cage_, "face_indices");

		cgogn::modeling::set_attribute_face_indices(*SpaceDeformationTool<MESH>::influence_cage_,
													cage_face_indices.get());

		std::shared_ptr<Attribute<uint32>> object_position_indices =
			get_attribute<uint32, Vertex>(object, "position_indices");

		std::shared_ptr<Attribute<uint32>> cage_position_indices =
			get_attribute<uint32, Vertex>(*control_cage_, "position_indices");

		std::shared_ptr<Attribute<uint32>> i_cage_position_indices =
			get_attribute<uint32, Vertex>(*SpaceDeformationTool<MESH>::influence_cage_, "position_indices");

		uint32 nbf_cage = 2 * nb_cells<Face>(*SpaceDeformationTool<MESH>::influence_cage_);
		uint32 nbv_cage = nb_cells<Vertex>(*SpaceDeformationTool<MESH>::influence_cage_);

		// first loop to find h
		//
		float h = 0.0f;
		float max_dist = 0.0f; 
		std::vector<Vec2> attenuation_points;
		SpaceDeformationTool<MESH>::influence_area_->foreach_cell([&](Vertex v) {
			uint32 surface_point_idx = value<uint32>(object, object_position_indices, v);

			float i_dist = SpaceDeformationTool<MESH>::cage_influence_distance(surface_point_idx, nbf_cage, nbv_cage);

			SpaceDeformationTool<MESH>::attenuation_(surface_point_idx) = (float)sin(0.5*M_PI * (i_dist ));
			/*if (control_area_validity_(surface_point_idx) == 1.0f)
			{
		
				if (i_dist > h)
				{
					h = i_dist;
				}

				SpaceDeformationTool<MESH>::attenuation_(surface_point_idx) = 1.0f;
			}
			else
			{
				
				if (i_dist > max_dist){
					max_dist = i_dist; 
				}
				attenuation_points.push_back({surface_point_idx, i_dist});
			}*/
		});

		
		/*for (unsigned int i = 0; i < attenuation_points.size(); i++)
		{
			//SpaceDeformationTool<MESH>::attenuation_(attenuation_points[i][0]) = 0.5f * ((float)sin(M_PI * ((attenuation_points[i][1] / h) - 0.5f))) + 0.5f;
			SpaceDeformationTool<MESH>::attenuation_(attenuation_points[i][0]) = (float)sin(0.5*M_PI * (attenuation_points[i][1] / h));
		}*/
	}

	bool local_mvc_pt_control_area(Vec3 pt)
	{
		std::shared_ptr<Attribute<uint32>> i_position_indices =
			get_attribute<uint32, Vertex>(*control_cage_, "position_indices");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked =
			get_attribute<bool, Vertex>(*control_cage_, "marked_vertices");

		DartMarker dm(*control_cage_);

		bool checked = true;
		for (Dart d = control_cage_->begin(), end = control_cage_->end(); d != end; d = control_cage_->next(d))
		{
			Vertex cage_vertex = CMap2::Vertex(d);
			bool vc_marked = value<bool>(*control_cage_, cage_vertex_marked, cage_vertex);

			if (!dm.is_marked(d) && !vc_marked)
			{

				const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, cage_vertex);

				float mvc_value = cgogn::modeling::compute_mvc(pt, d, *control_cage_, cage_point, control_cage_vertex_position_.get());

				dm.mark(d);

				value<bool>(*control_cage_, cage_vertex_marked, cage_vertex) = true;

				if (mvc_value < 0)
				{
					checked = false;
					break;
				}
			}
		}

		parallel_foreach_cell(*control_cage_, [&](Vertex vc) -> bool {
			value<bool>(*control_cage_, cage_vertex_marked, vc) = false;
			return true;
		});

		return checked;
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
