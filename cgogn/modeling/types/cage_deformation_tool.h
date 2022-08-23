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

#include <cgogn/modeling/types/space_deformation_tool.h>
#include <cgogn/core/types/cmap/cmap2.h>


#include <cgogn/modeling/algos/deformation/deformation_cage.h>


namespace cgogn
{

namespace modeling
{

template <typename MESH>
class CageDeformationTool: public SpaceDeformationTool<MESH>
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
    using Vec3 = geometry::Vec3;

public:

	MESH* control_cage_; 

	CageDeformationTool(): SpaceDeformationTool<MESH>(), m_hFactor(-1.0f), control_cage_vertex_position_(nullptr)
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

		std::shared_ptr<Attribute<uint32>> position_indices = cgogn::add_attribute<uint32, Vertex>(*control_cage_, "position_indices");
		cgogn::modeling::set_attribute_position_indices(*control_cage_, position_indices.get()); 
    }

	void set_influence_cage(MESH* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max){
		SpaceDeformationTool<MESH>::influence_cage_ = m; 
		cgogn::modeling::create_box(*m, vertex_position, bb_min, bb_max); 

		SpaceDeformationTool<MESH>::influence_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position"); 

		std::shared_ptr<Attribute<uint32>> position_indices = cgogn::add_attribute<uint32, Vertex>(*SpaceDeformationTool<MESH>::influence_cage_, "position_indices");
		cgogn::modeling::set_attribute_position_indices(*SpaceDeformationTool<MESH>::influence_cage_, position_indices.get()); 

		std::shared_ptr<Attribute<bool>> marked_vertices = cgogn::add_attribute<bool, Vertex>(*SpaceDeformationTool<MESH>::influence_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*SpaceDeformationTool<MESH>::influence_cage_, marked_vertices.get()); 
	}
	

	void set_center_control_cage(Vec3& center){
		center_control_cage_ = center; 
	}

	void update_influence_area(const MESH& m, CMap2::Attribute<Vec3>* vertex_position){
		
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(m, vertex_position, v);

			bool inside_cage = local_mvc_pt_surface(surface_point);

			if (inside_cage)
			{
				SpaceDeformationTool<MESH>::influence_area_->select(v);
			}

			return true;
		});
	}


private:

	std::shared_ptr<Attribute<Vec3>> control_cage_vertex_position_; 

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> ctrl_cage_coords_;

	Eigen::VectorXd mixing_factor_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_comp_;

	float m_hFactor;

	Vec3 center_control_cage_; 

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;


	bool local_mvc_pt_surface(Vec3 pt)
	{
		std::shared_ptr<Attribute<uint32>> i_position_indices = get_attribute<uint32, Vertex>(*SpaceDeformationTool<MESH>::influence_cage_, "position_indices");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked = get_attribute<bool, Vertex>(*SpaceDeformationTool<MESH>::influence_cage_, "marked_vertices");

		DartMarker dm(*SpaceDeformationTool<MESH>::influence_cage_);

		bool checked = true;
		for (Dart d = SpaceDeformationTool<MESH>::influence_cage_->begin(), end = SpaceDeformationTool<MESH>::influence_cage_->end(); d != end; d = SpaceDeformationTool<MESH>::influence_cage_->next(d))
		{
			Vertex cage_vertex = CMap2::Vertex(d);
			bool vc_marked = value<bool>(*SpaceDeformationTool<MESH>::influence_cage_, cage_vertex_marked, cage_vertex);

			if (!dm.is_marked(d) && !vc_marked)
			{

				const Vec3& cage_point = value<Vec3>(*SpaceDeformationTool<MESH>::influence_cage_, SpaceDeformationTool<MESH>::influence_cage_vertex_position_, cage_vertex);

				float mvc_value = compute_mvc(pt, d, *SpaceDeformationTool<MESH>::influence_cage_, cage_point, SpaceDeformationTool<MESH>::influence_cage_vertex_position_.get());

				dm.mark(d);

				value<bool>(*SpaceDeformationTool<MESH>::influence_cage_, cage_vertex_marked, cage_vertex) = true;

				if (mvc_value < 0)
				{
					checked = false;
					break;
				}
			}
		}

		parallel_foreach_cell(*SpaceDeformationTool<MESH>::influence_cage_, [&](Vertex vc) -> bool {
			value<bool>(*SpaceDeformationTool<MESH>::influence_cage_, cage_vertex_marked, vc) = false;
			return true;
		});

		return checked;
	}
	
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
