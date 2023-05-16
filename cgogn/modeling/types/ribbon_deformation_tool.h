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

#ifndef CGOGN_MODELING_RIBBON_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_RIBBON_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cells_set.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

#include <cgogn/geometry/algos/normal.h>

#include <cgogn/modeling/types/handle_deformation_tool.h>
#include <cgogn/modeling/types/axis_deformation_tool.h>
#include <cgogn/modeling/types/cage_deformation_tool.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>

/**
 * @Class CageDeformationTool
 * Represents the local cage deformation tool
 */
class RibbonDeformationTool
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Graph = IncidenceGraph;

	template <typename T>
	using GraphAttribute = typename mesh_traits<IncidenceGraph>::
													template Attribute<T>;
	using GraphVertex = IncidenceGraph::Vertex;
	using GraphEdge = IncidenceGraph::Edge;
	using GraphFace = IncidenceGraph::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;


public:
	MESH* control_ribbon_;
	std::shared_ptr<Attribute<Vec3>> control_ribbon_vertex_position_;

	std::shared_ptr<boost::synapse::connection> 
								ribbon_attribute_update_connection_;

	MatrixWeights object_weights_;

	Eigen::Matrix3d local_frame_;  

	Vec3 control_ribbon_local_bb_min_; 
	Vec3 control_ribbon_local_bb_max_; 

	std::shared_ptr<Attribute<uint32>> control_ribbon_vertex_index_;

	RibbonDeformationTool() : control_ribbon_vertex_position_(nullptr)
	{
		
	}

	~RibbonDeformationTool()
	{
	}

	/// @brief create local cage mesh 
	/// initialize the triangles set of this cage 
	/// initialize the local_direction_control_planes of this local cage
	/// generate the virtual cubes
	/// TODO: see to generate the virtual cubes only if needed 
	/// @param m default mesh 
	/// @param m_vertex_position position of the vertices of the default mesh 
	/// @param bb_min bounding box minimum position
	/// @param bb_max bounding box maximum position
	/// @param normal bounding box normal 
	void create_space_tool(MESH* m, 
					CMap2::Attribute<Vec3>* m_vertex_position, 
					const Vec3& bb_min, const Vec3& bb_max,
					const std::tuple<Vec3, Vec3, Vec3>& main_directions)
	{
		control_ribbon_ = m;

		local_frame_.row(0) = std::get<0>(main_directions);
		local_frame_.row(1) = std::get<1>(main_directions);
		local_frame_.row(2) = std::get<2>(main_directions);

		create_cylinder_tool(*m, m_vertex_position, bb_min, bb_max); 

		control_ribbon_vertex_position_ = 
						cgogn::get_attribute<Vec3, Vertex>(*control_ribbon_, "position");

		control_ribbon_vertex_index_ = 
					cgogn::add_attribute<uint32, Vertex>(*control_ribbon_, 
														"vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*control_ribbon_, 
										control_ribbon_vertex_index_.get());

		control_ribbon_bb_min_ = bb_min;
		control_ribbon_bb_max_ = bb_max;

		control_ribbon_local_bb_min_ = local_frame_ * bb_min; 
		control_ribbon_local_bb_max_ = local_frame_ * bb_max; 

		set_control_ribbon_init_positions(); 
	}

	void reset_deformation()
	{
		foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
			uint32 ribbon_point_index = value<uint32>(*control_ribbon_, 
											control_ribbon_vertex_index_, cv);

			value<Vec3>(*control_ribbon_, control_ribbon_vertex_position_, cv) = 
											init_positions_[ribbon_point_index]; 

			return true;
		});
	}

	/// @brief set the center of the control cage 
	/// @param center 
	void set_center_control_cage(Vec3& center)
	{
		control_ribbon_center_ = center;
	}

	/// @brief initialize the binding of the model to deform 
	/// so far only MVC deformation 
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model 
	/// @param object_vertex_index index of the vertices of the model 
	void init_bind_object(MESH& object, 
				CMap2::Attribute<Vec3>* object_vertex_position, 
				CMap2::Attribute<uint32>* object_vertex_index)
	{
		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		object_weights_.position_.resize(nbv_object, nbv_cage);
		object_weights_.position_.setZero();
	}


	/// @brief update the binding of the model 
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param object_vertex_index index of the vertices of the model
	void update_bind_object(MESH& object, 
				CMap2::Attribute<Vec3>* object_vertex_position, 
				CMap2::Attribute<uint32>* object_vertex_index)
	{

	}


	/// @brief deform the model following the chosen deformation type 
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param object_vertex_index index of the vertices of the model
	void deform_object(MESH& object,
					CMap2::Attribute<Vec3>* object_vertex_position,
					CMap2::Attribute<uint32>* object_vertex_index, 
					std::unordered_map<std::string, 
					std::shared_ptr<modeling::HandleDeformationTool<MESH>>>& handle_container, 
					std::unordered_map<std::string, 
					std::shared_ptr<modeling::AxisDeformationTool<MESH>>>& axis_container, 
					std::unordered_map<std::string, 
					std::shared_ptr<modeling::CageDeformationTool<MESH>>>& cage_container)
	{
		const Vec3 displacement = Vec3();  
		foreach_cell(object, [&](Vertex v) -> bool {
			uint32 object_point_index = 
							value<uint32>(object, object_vertex_index, v);

			value<Vec3>(object, object_vertex_position, v) += 
								object_weights[object_point_index]*displacement;
			
			return true; 
			}); 

	}


private:
	std::vector<CMap2::Vertex> control_ribbon_vertices_; 

	Vec3 control_ribbon_center_;
	Vec3 control_ribbon_bb_min_;
	Vec3 control_ribbon_bb_max_;

	std::vector<Vec3> init_positions_;


	/// @brief set init positions to reset the deformation
	///
	void set_control_ribbon_init_positions()
	{
		uint32 nbv_cage = nb_cells<Vertex>(*control_ribbon_);

		init_positions_.resize(nbv_cage);

		foreach_cell(*control_ribbon_, [&](Vertex cv) -> bool {
			const Vec3& ribbon_point = value<Vec3>(*control_ribbon_, 
										control_ribbon_vertex_position_, cv);

			uint32 ribbon_point_index = value<uint32>(*control_ribbon_, 
											control_ribbon_vertex_index_, cv);

			init_positions_[ribbon_point_index] = ribbon_point; 

			return true;
		});
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_RIBBON_DEFORMATION_TOOL_H_
