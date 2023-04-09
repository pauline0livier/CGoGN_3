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

#include <cgogn/modeling/types/dual_quaternion.h>

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

	Eigen::SparseMatrix<double, Eigen::RowMajor> object_weights_; 

	std::shared_ptr<Graph::Attribute<Vec3>> control_axis_vertex_position_;
	std::shared_ptr<boost::synapse::connection> 
										axis_attribute_update_connection_;

	Eigen::MatrixXf global_cage_weights_;
	Eigen::MatrixXf global_cage_normal_weights_;

	std::vector<MeshVertex> object_influence_area_; 

	std::string deformation_type_; 

	AxisDeformationTool(): control_axis_vertex_position_(nullptr)
	{
	}

	~AxisDeformationTool()
	{
	}

	/// @brief create axis space tool
	/// set axis from user's selection
	/// create clone axis inside the model to deform 
	/// set the tranformations to identity 
	/// @param g default graph 
	/// @param g_vertex_position position of the vertices of the default graph
	/// @param vertex_radius radius of the vertices of the default graph 
	/// @param vertex_coordinates positions to set on the vertices
	/// @param vertex_normals normals to set on the vertices 
	/// @param inside_axis_position positions inside the mesh 
	void create_space_tool(Graph* g, 
			Graph::Attribute<Vec3>* g_vertex_position, 
			Graph::Attribute<Scalar>* g_vertex_radius,
			const std::vector<Vec3>& vertex_coordinates, 
			const std::vector<Vec3>& vertex_normals, 
			const std::vector<Vec3>& inside_axis_position)
	{
		control_axis_ = g;

		axis_skeleton_ = 
			cgogn::modeling::create_axis(*g, g_vertex_position, 
										g_vertex_radius, vertex_coordinates, Scalar(3));

		control_axis_vertex_position_ = 
				cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		inside_axis_position_ = inside_axis_position; 

		std::shared_ptr<Graph::Attribute<uint32>> a_vertex_index =
			cgogn::add_attribute<uint32, Graph::Vertex>(*control_axis_,
															 "vertex_index");

		cgogn::modeling::set_attribute_vertex_index_graph(*control_axis_, 
														a_vertex_index.get());

		size_t number_of_transformations = axis_skeleton_.size() +1; 
		axis_transformation_.resize(number_of_transformations); 
		dual_quaternion_transformations_.resize(number_of_transformations); 

		for (size_t t = 0; t < axis_transformation_.size(); t++){
			axis_transformation_[t].setIdentity();

			dual_quaternion_transformations_[t] = modeling::DualQuaternion(); 
		}

		starting_positions_ = vertex_coordinates; 
		normal_ = {0.0, 0.0, 1.0}; 
	}

	/// @brief set axis transformation from user's input 
	/// @param axis_transformation set of tranformations (one per bone)
	void set_axis_transformation(
			const std::vector<rendering::Transfo3d>& axis_transformation){

		axis_transformation_ = axis_transformation; 
	}

	void set_dual_quaternions_transformation(const std::vector<DualQuaternion>& dual_quat_transformation)
	{
		dual_quaternion_transformations_ = dual_quat_transformation; 
	}

	/// @brief set deformation type 
	/// so far between LBS or DQS deformation 
	/// @param new_type 
	void set_deformation_type(const std::string new_type){
		deformation_type_ = new_type; 
	}

	/// @brief reset deformation
	/// useful to change deformation type
	/// TODO: find transformation from current to init positions
	void reset_deformation()
	{

		/*for (size_t v = 0; v < axis_skeleton_.size(); v++)
		{

			Vec3 init_position = starting_positions_[v].normalized(); 
			Vec3 current_position = value<Vec3>(*control_axis_, control_axis_vertex_position_, axis_skeleton_[v]).normalized(); 

			if (init_position == current_position)
			{
				std::cout << "identity" << std::endl; 
			} else {
				Vec3 cross_product = current_position.cross(init_position); 

				double sine = cross_product.norm(); 

				double delta = asin(sine); 

				std::cout << delta << std::endl; 


			}

			

		}  */

	}

	/// @brief initialize binding of object 
	/// so far only rigid deformation type 
	/// @param object 
	/// @param object_vertex_position 
	void init_bind_object(MESH& object, 
				const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		uint32 nb_bones = axis_skeleton_.size() +1;
		object_weights_ = Eigen::SparseMatrix<double, Eigen::RowMajor>(nbv_object, nb_bones); 

		
		bind_object_rigid(object, object_vertex_position); 
		
	}

	/// @brief update binding of object 
	/// so far only rigid deformation type 
	/// @param object 
	/// @param object_vertex_position 
	void bind_object(MESH& object, 
				const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{
		object_weights_.setZero(); 

		bind_object_rigid(object, object_vertex_position);
		
	}

	void deform_object(MESH& object, 
					CMap2::Attribute<Vec3>* object_vertex_position, 
					CMap2::Attribute<uint32>* object_vertex_index)
	{
		if (deformation_type_ == "LBS")
		{
			deform_object_LBS(object, object_vertex_position, 
								object_vertex_index); 
		}

		if (deformation_type_ == "DQS")
		{
			std::shared_ptr<Attribute<Vec3>> object_vertex_normal =
			get_attribute<Vec3, MeshVertex>(object, "normal");

			deform_object_DQS(object, object_vertex_position, 
								object_vertex_index, object_vertex_normal.get());
		}
	}

	

	

	/// @brief 
	/// @return vector composed of the vertices of the axis
	std::vector<Graph::Vertex> get_axis_skeleton(){
		return axis_skeleton_; 
	}

	
private:
	std::vector<Graph::Vertex> axis_skeleton_; 
	Vec3 axis_normal_;

	std::vector<rendering::Transfo3d> axis_transformation_; 

	std::vector<modeling::DualQuaternion> dual_quaternion_transformations_; 

	std::vector<Vec3> inside_axis_position_; 

	std::vector<Vec3> starting_positions_; 

	Vec3 normal_; 

	/// @brief bind the zone of influence of the object  
	/// compute weights using the projection of the points on the axis
	/// @param object 
	/// @param vertex_position 
	void bind_object_rigid(MESH& object, 
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

			object_weights_.row(surface_point_index) = 
				cgogn::modeling::weight_partial_skeleton(inside_axis_position_, 
												surface_point);
		
		}

	}

	/// @brief deform the assigned zone of influence of the object
	/// follow Linear Blend Skinning (LBS) formalism  
	/// @param object model to deform
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	void deform_object_LBS(MESH& object, 
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

			Vec3 new_position = {0.0, 0.0, 0.0}; 

			for (size_t t = 0; t < axis_transformation_.size(); t++){
 
				Vec3 local_transformation = axis_transformation_[t]*v_position; 

				new_position += object_weights_.coeff(v_index, t)*
											local_transformation;

			}

			value<Vec3>(object, object_vertex_position, v) = new_position; 

		}

	}

	/// @brief deform the assigned zone of influence of the object
	/// follow Dual Quaternion Skinning (DQS) formalism  
	/// @param object model to deform
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	void deform_object_DQS(MESH& object, 
					CMap2::Attribute<Vec3>* object_vertex_position, 
					CMap2::Attribute<uint32>* object_vertex_index, 
					CMap2::Attribute<Vec3>* object_vertex_normal)
	{

		const std::size_t influence_area_length = 
									this->object_influence_area_.size(); 

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			uint32 v_index = value<uint32>(object, object_vertex_index, v);
			Vec3 v_position = value<Vec3>(object, object_vertex_position, v); 
			Vec3 v_normal = value<Vec3>(object, object_vertex_normal, v);

			int  k0 = -1;
			double w0 = 0.f;
			DualQuaternion dq_blend;
			Eigen::Quaterniond q0;

			/*if(axis_skeleton_.size() != 0)
			{
				k0 = axis_skeleton_[0].index_;
				w0 = object_weights_.coeff(v_index, 0);
			}else
				dq_blend = DualQuaternion::identity();

			if(k0 != -1) dq_blend = dual_quaternion_transformations_[k0] * w0;

			int pivot = k0;*/

			//q0 = dual_quaternion_transformations_[pivot].rotation();

			for (size_t t = 0; t < dual_quaternion_transformations_.size(); t++){
 
				//const int k = axis_skeleton_[j].index_;
				double w = object_weights_.coeff(v_index, t);
				//const DualQuaternion& dq = (k == -1) ? DualQuaternion::identity() : dual_quaternion_transformations_[k];
				const DualQuaternion& dq = dual_quaternion_transformations_[t]; 

				/*if( dq.rotation().dot( q0 ) < 0.f )
					w *= -1.0;*/

				dq_blend = dq_blend + dq * w;

			}

			Vec3 vi = dq_blend.transform( v_position );
			
			value<Vec3>(object, object_vertex_position, v) = vi;
			
			dq_blend.normalize(); 
			value<Vec3>(object, object_vertex_normal, v) = dq_blend.rotate_quaternion( dq_blend.get_non_dual_part(),v_normal );

		}

	}


};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_AXIS_DEFORMATION_H_
