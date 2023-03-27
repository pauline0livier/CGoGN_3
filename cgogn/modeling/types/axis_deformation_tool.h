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
 * Inc., 51 Franklin Street, Fifth Floor, Bosston, MA  02110-1301 USA.
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *                            
 *******************************************************************************/

#ifndef CGOGN_MODELING_AXIS_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_AXIS_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cells_set.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>
#include <cgogn/rendering/types.h>

namespace cgogn
{

namespace modeling
{

/**
 * @Class Axis Deformation Tool 
 * Axis represented as an Incidence graph 
 * Create axis, deform and propagate deformation to object 
 * and associate tools 
*/
template <typename MESH>
class AxisDeformationTool
{

	using Graph = cgogn::IncidenceGraph;

public:
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using MeshVertex = typename mesh_traits<MESH>::Vertex;
	using MeshFace = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	Graph* control_axis_;

	Eigen::SparseMatrix<double, Eigen::RowMajor> axis_weights_; 

	std::shared_ptr<Graph::Attribute<Vec3>> control_axis_vertex_position_;
	std::shared_ptr<boost::synapse::connection> 
										axis_attribute_update_connection_;

	Eigen::MatrixXf global_cage_weights_;
	Eigen::MatrixXf global_cage_normal_weights_;

	std::vector<MeshVertex> object_influence_area_; 

	AxisDeformationTool(): control_axis_vertex_position_(nullptr)
	{
	}

	~AxisDeformationTool()
	{
	}

	/**
	* create axis space tool 
	*	set axis from user's selection
	*  	create clone axis inside the object 
	*/
	void create_space_tool(Graph* g, 
			Graph::Attribute<Vec3>* vertex_position, 
			Graph::Attribute<Scalar>* vertex_radius,
			const std::vector<Vec3>& vertex_coords, 
			const std::vector<Vec3>& vertex_normals, 
			const std::vector<Vec3>& inside_axis_position)
	{
		control_axis_ = g;
		axis_skeleton_ = 
			cgogn::modeling::create_axis(*g, vertex_position, 
										vertex_radius, vertex_coords);

		control_axis_vertex_position_ = 
				cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		inside_axis_position_ = inside_axis_position; 

		std::shared_ptr<Graph::Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, Graph::Vertex>(*control_axis_,
															 "vertex_index");

		cgogn::modeling::set_attribute_vertex_index_graph(*control_axis_, 
														vertex_index.get());
	}

	/**
	 * set deformation type 
	 * so far between rigid or loose deformation 
	*/
	void set_deformation_type(const std::string new_type){
		deformation_type_ = new_type; 
	}

	/**
	 * bind axis to object 
	 * rigid deformation 
	 * weights inspired from linear blend skinning weights 
	*/
	void set_binding_rigid(MESH& object, 
				const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		uint32 nb_bones = axis_skeleton_.size() +1;
		axis_weights_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(nbv_object, nb_bones); 

		compute_weights(object, vertex_position);
	}

	/**
	 * bind axis to object 
	 * loose deformation  
	*/
	void set_binding_loose(MESH& object, 
				const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		uint32 nb_bones = axis_skeleton_.size() +1;

		axis_weights_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(nbv_object, nb_bones); 

		compute_weights(object, vertex_position);
	}

	/**
	 * set axis transformation from user's input 
	 * @param {std::vector<rendering::Transfo3d} axis_transformation
	 * 	transformation of each axis's bone
	*/
	void set_axis_transformation(
			const std::vector<rendering::Transfo3d>& axis_transformation){

		axis_transformation_ = axis_transformation; 
	}

	/**
	 * deform the influence area on the object 
	 * use the weights and Linear Blend Skinning formalism 
	 * check if weights on virtual bones with axis_fixed_point value
	*/
	void deform_object(MESH& object, 
					CMap2::Attribute<Vec3>* object_vertex_position, 
					CMap2::Attribute<uint32>* object_vertex_index)
	{
		 
		const std::size_t influence_area_length = 
									this->object_influence_area_.size(); 

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			uint32 v_index = value<uint32>(object, object_vertex_index, v);
			Vec3 v_position = value<Vec3>(object, object_vertex_position, v);

			for (size_t t = 0; t < axis_transformation_.size(); t++){
				Vec3 local_transformation = axis_transformation_[t]*v_position; 

				value<Vec3>(object, object_vertex_position, v) += 
									axis_weights_.coeff(v_index, t)*local_transformation; 

			}

		}

	}

	/**
	* get axis skeleton
	* @returns {std::vector<Graph::Vertex>} axis vertices 
	*/
	std::vector<Graph::Vertex> get_axis_skeleton(){
		return axis_skeleton_; 
	}

	
private:
	std::vector<Graph::Vertex> axis_skeleton_; 
	Vec3 axis_normal_;
	Eigen::Matrix3d local_frame_;
	Eigen::Matrix3d local_frame_inverse_;

	std::string deformation_type_; 

	std::vector<rendering::Transfo3d> axis_transformation_; 

	std::vector<Vec3> inside_axis_position_; 

	/**
	 * compute weights of the influence area on the object 
	*/
	void compute_weights(MESH& object, 
		const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		const std::size_t influence_area_length = 
								this->object_influence_area_.size(); 

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			Vec3 surface_point = value<Vec3>(object, vertex_position, v);

			uint32 surface_point_index = 
							value<uint32>(object, object_vertex_index, v);

			axis_weights_.row(surface_point_index) = 
				cgogn::modeling::weight_partial_skeleton(inside_axis_position_, 
												surface_point);
		
		}

	}


};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_AXIS_DEFORMATION_H_
