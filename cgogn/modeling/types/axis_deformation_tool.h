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

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> axis_weights_; 
	Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> axis_fixed_point_; 

	std::shared_ptr<Graph::Attribute<Vec3>> control_axis_vertex_position_;
	std::shared_ptr<boost::synapse::connection> axis_attribute_update_connection_;

	Eigen::MatrixXf global_cage_weights_;
	Eigen::MatrixXf global_cage_normal_weights_;

	std::vector<MeshVertex> object_influence_area_; 

	AxisDeformationTool(): control_axis_vertex_position_(nullptr)
	{
	}

	~AxisDeformationTool()
	{
	}

	void create_space_tool(Graph* g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius,
						   const std::vector<Vec3>& vertex_coords, const std::vector<Vec3>& vertex_normals, const std::vector<Vec3>& inside_axis_position)
	{
		control_axis_ = g;
		axis_skeleton_ = cgogn::modeling::create_axis(*g, vertex_position, vertex_radius, vertex_coords);

		control_axis_vertex_position_ = cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		inside_axis_position_ = inside_axis_position; 

		std::shared_ptr<Graph::Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, Graph::Vertex>(*control_axis_, "vertex_index");

		cgogn::modeling::set_attribute_vertex_index_graph(*control_axis_, vertex_index.get());
	}

	void set_deformation_type(const std::string new_type){
		deformation_type_ = new_type; 
	}

	void set_binding_rigid(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		uint32 nb_bones = axis_skeleton_.size() -1; 
		axis_weights_.resize(nbv_object, nb_bones);
		axis_weights_.setZero();

		axis_fixed_point_.resize(nbv_object, nb_bones);
		axis_fixed_point_.setZero();

		compute_weights(object, vertex_position);
	}

	void set_binding_loose(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		uint32 nb_bones = axis_skeleton_.size() -1;

		axis_weights_.resize(nbv_object, nb_bones);
		axis_weights_.setZero();

		axis_fixed_point_.resize(nbv_object, nb_bones);
		axis_fixed_point_.setZero();

		compute_weights(object, vertex_position);
	}

	void set_axis_transformation(const std::vector<rendering::Transfo3d>& axis_transformation){
		axis_transformation_ = axis_transformation; 
	}

	void deform_object(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position, CMap2::Attribute<uint32>* object_vertex_index)
	{
		 
		const std::size_t influence_area_length = this->object_influence_area_.size(); 

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			uint32 v_index = value<uint32>(object, object_vertex_index, v);
			Vec3 v_position = value<Vec3>(object, object_vertex_position, v);

			Vec3 left_influence = axis_transformation_[0]*v_position; 
			Vec3 right_influence = axis_transformation_[1]*v_position;

			if (axis_fixed_point_(v_index,0)){
				left_influence = v_position; 
			}

			if (axis_fixed_point_(v_index,1)){
				right_influence = v_position; 
			}

			double weight0 = axis_weights_(v_index, 0);

			double weight1 = axis_weights_(v_index, 1);

			value<Vec3>(object, object_vertex_position, v) = weight0 * left_influence + weight1 * right_influence;


			/*value<Vec3>(object, object_vertex_position, v) += attenuation_[vidx] * new_deformation;*/
		}

	}

	std::vector<Graph::Vertex> get_axis_skeleton(){
		return axis_skeleton_; 
	}

	
private:
	std::vector<Graph::Vertex> axis_skeleton_;
	//Graph::Vertex axis_center_; 
	Vec3 axis_normal_;
	Eigen::Matrix3d local_frame_;
	Eigen::Matrix3d local_frame_inverse_;

	std::string deformation_type_; 

	std::vector<rendering::Transfo3d> axis_transformation_; 

	std::vector<Vec3> inside_axis_position_; 

	void compute_weights(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		const std::size_t influence_area_length = this->object_influence_area_.size(); 

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			Vec3 surface_point = value<Vec3>(object, vertex_position, v);

			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			std::pair<Eigen::Vector2d, std::vector<bool>> result_weights = cgogn::modeling::weight_two_bones(inside_axis_position_[0], inside_axis_position_[1], inside_axis_position_[2], surface_point); 

			//axis_weights_.row(surface_point_idx) = local_weight; 
			axis_weights_(surface_point_idx, 0) = result_weights.first[0]; 
			axis_weights_(surface_point_idx, 1) = result_weights.first[1];

			axis_fixed_point_(surface_point_idx, 0) = result_weights.second[0]; 
			axis_fixed_point_(surface_point_idx, 1) = result_weights.second[1];
		}

	}


};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_H_


// DeBUG

/*void init_space_tool(Graph* g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius,
						   const Vec3& vertex_coord, const Vec3& vertex_normal)
	{
		control_axis_ = g;
		axis_center_ = cgogn::modeling::create_handle(*g, vertex_position, vertex_radius, vertex_coord);

		control_axis_vertex_position_ = cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");
	}*/

	/*void set_influence_cage_axis(MESH* m, 
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

		std::shared_ptr<Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, MeshVertex>(*(this->influence_cage_), "vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*(this->influence_cage_),
														vertex_index.get());

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, MeshVertex>(*(this->influence_cage_), "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*(this->influence_cage_),
													   marked_vertices.get());
	}*/
