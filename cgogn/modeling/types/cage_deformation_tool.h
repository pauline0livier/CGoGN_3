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

#include <cgogn/core/types/cells_set.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

#include <cgogn/geometry/algos/normal.h>

#include <cgogn/modeling/types/handle_deformation_tool.h>
#include <cgogn/modeling/types/axis_deformation_tool.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>

/**
 * @Class CageDeformationTool
 * Represents the local cage deformation tool
 */
class CageDeformationTool
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

	/**
	 * point structure
	 * point can be on the local cage
	 * if it is the case: 
	 * 		- inside_control_cage = true,
	 * 		- control_cage_index and vertex defined, 
	 * otherwise:
	 * 		- inside_control_cage = false, 
	 * 		- vertex of closest control cage point 
	 * 		- shift_vector defined (a,b,c)
	 * 		a,b,c in {-1, 0, 1}, encode the direction and if after or before  
	 * 		TODO: check if update is activated if the point is common to 
	 * 				at least two cages 
	*/
	struct Point
	{
		Vec3 position;
		uint32_t control_cage_index;
		bool inside_control_cage;
		Vertex vertex; 
		Vec3 shift_vector; 
		bool update; 
	};

	/**
	 * triangle structure
	 * virtual_cage_indices for virtual cubes 
	*/
	struct Triangle
	{
		std::vector<Point> points;
		std::vector<uint32_t> virtual_cage_indices;
		Vec3 normal;
		std::pair<Vec3, Vec3> edges;
	};

	/**
	 * virtual cube structure 
	 * can contain point belonging 
	 * to the control cage 
	 * Constraint: not all points belong to the local cage
	 * 
	*/
	struct Virtual_cube
	{
		std::vector<Triangle> triangles; 
		std::vector<Point> points;
	};

	/**
	 * delimit the local cage in terms of plane intersections
	 * Each local_direction_control_planes is specific to a direction
	 * line equation a*x + b*y + c*z = d
	 * d_min and d_max for the extrema heights 
	 * triangles_d_min and triangles_d_max represent the square faces 
	 * shift_after_d_max, translation  
	 * shift_before_d_min, translation 
	 * 	that is the position of the next virtual cube planes along this direction
	*/
	struct Local_direction_control_planes
	{
		double d_min;
		double d_max;
		double d_gap;

		std::pair<Triangle, Triangle> triangles_d_min;
		std::pair<Triangle, Triangle> triangles_d_max;

		Vec3 direction;
		Vec3 shift_after_d_max;
		Vec3 shift_before_d_min;
	};

	struct Handle_data
	{
		VectorWeights weights_; 
		int cage_index_; 
	}; 

	struct Tool_data
	{
		MatrixWeights weights_; 
		std::unordered_map<uint32, int> cage_index_; 
	}; 

public:
	MESH* control_cage_;
	std::shared_ptr<Attribute<Vec3>> control_cage_vertex_position_;

	std::shared_ptr<Attribute<Vec3>> influence_cage_vertex_position_;

	std::shared_ptr<boost::synapse::connection> 
								cage_attribute_update_connection_;

	Eigen::VectorXd attenuation_;

	MatrixWeights object_weights_;
	Eigen::VectorXi object_activation_cage_; 

	Eigen::Matrix3d local_frame_;  

	Vec3 cage_local_bb_min_; 
	Vec3 cage_local_bb_max_; 

	std::unordered_map<std::string, Handle_data> local_handle_data_;
	std::unordered_map<std::string, Tool_data> local_cage_data_; 

	
	std::string deformation_type_;

	std::string split_deformation_type_; 

	bool need_full_bind_; 

	std::shared_ptr<Attribute<uint32>> control_cage_vertex_index_;

	CageDeformationTool() : control_cage_vertex_position_(nullptr)
	{
		virtual_cubes_.resize(26); 
	}

	~CageDeformationTool()
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
		control_cage_ = m;

		Vec3 temp_min = bb_min; 
		Vec3 temp_max = bb_max; 
		if (bb_min[2] == bb_max[2] ){
			temp_min[2] = -0.1; 
			temp_max[2] = 0.1; 
		}

		local_frame_.row(0) = std::get<0>(main_directions);
		local_frame_.row(1) = std::get<1>(main_directions);
		local_frame_.row(2) = std::get<2>(main_directions);

		//create_cage_box(*m, m_vertex_position, bb_min, bb_max, local_frame_); 
		control_cage_vertices_ = create_cage_box(*m, m_vertex_position, temp_min, temp_max, local_frame_); 

		control_cage_vertex_position_ = 
						cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		control_cage_bb_min_ = bb_min;
		control_cage_bb_max_ = bb_max;

		cage_local_bb_min_ = local_frame_ * bb_min; 
		cage_local_bb_max_ = local_frame_ * bb_max; 

		control_cage_vertex_index_ = 
					cgogn::add_attribute<uint32, Vertex>(*control_cage_, 
														"vertex_index");
		
		cgogn::modeling::set_attribute_vertex_index(*control_cage_, 
										control_cage_vertex_index_.get());


		init_triangles();

		set_start_positions(); 

		init_control_cage_plane(main_directions, 5); //3, 25

		init_virtual_cubes();
	}

	/// @brief set the deformation type 
	/// @param new_type so far MVC or Green
	void set_deformation_type(const std::string& new_type)
	{
		deformation_type_ = new_type;
	}

	/// @brief set the split deformation type 
	/// @param new_type by strip or attenuation
	void set_split_deformation_type(const std::string& new_type)
	{
		split_deformation_type_ = new_type; 
	}

	/// @brief require full re-binding of the local cage
	/// happen when the local cage position are updated
	/// outside of the local cage deformation
	void require_full_binding()
	{
		need_full_bind_ = true; 
	}


	void reset_deformation()
	{
		foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
			uint32 cage_point_index = value<uint32>(*control_cage_, 
											control_cage_vertex_index_, cv);

			value<Vec3>(*control_cage_, control_cage_vertex_position_, cv) = start_positions_[cage_point_index]; 

			return true;
		});
	}

	/// @brief set the center of the control cage 
	/// @param center 
	void set_center_control_cage(Vec3& center)
	{
		control_cage_center_ = center;
	}

	/// @brief update virtual cages after local cage deformation 
	void update_virtual_cages(){
		if (split_deformation_type_ == "Strip")
		{
			update_virtual_cages_strip(); 
		} 

		if (split_deformation_type_ == "Attenuation")
		{
			update_virtual_cages_attenuation(); 
		}
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

		object_activation_cage_.resize(nbv_object); 
		object_activation_cage_.setZero(); 

		if (deformation_type_ == "MVC")
		{
			bind_object_mvc(object, object_vertex_position, 
							object_vertex_index);
		}

		if (deformation_type_ == "Green"){
			uint32 nbt_cage = cage_triangles_.size();
			object_weights_.normal_.resize(nbv_object, nbt_cage);
			object_weights_.normal_.setZero();

			bind_object_green(object, object_vertex_position, 
							object_vertex_index);
		}
	}


	/// @brief update the binding of the model 
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model
	/// @param object_vertex_index index of the vertices of the model
	void update_bind_object(MESH& object, 
				CMap2::Attribute<Vec3>* object_vertex_position, 
				CMap2::Attribute<uint32>* object_vertex_index)
	{

		if (deformation_type_ == "MVC")
		{
			update_bind_object_mvc(object, object_vertex_position, 
							object_vertex_index);
		}
	}

	void bind_connecting_tools(std::unordered_map<std::string, 
					std::shared_ptr<modeling::HandleDeformationTool<MESH>>>& handle_container, 
					std::unordered_map<std::string, 
					std::shared_ptr<modeling::AxisDeformationTool<MESH>>>& axis_container, 
					std::unordered_map<std::string, 
					std::shared_ptr<modeling::CageDeformationTool<MESH>>>& cage_container )
	{
		if (local_handle_data_.size() > 0)
		{
			for (auto& [name_handle, data] : local_handle_data_)
			{
				std::shared_ptr<modeling::HandleDeformationTool<MESH>> local_hdt = handle_container[name_handle];

				Vec3 handle_position = local_hdt->get_handle_position(); 

				update_bind_handle(name_handle, handle_position);
			}
		}

		if (local_cage_data_.size() > 0)
		{
			for (auto& [name_cage, cData] : local_cage_data_)
			{
				std::shared_ptr<modeling::CageDeformationTool<MESH>> local_cdt = cage_container[name_cage];

				update_bind_local_cage(*(local_cdt->control_cage_), name_cage, 
									local_cdt->control_cage_vertex_position_, local_cdt->control_cage_vertex_index_);
			}
				
		}
	}

	void update_local_bounding_box()
	{
		Vec3 local_min = {1000.0, 1000.0, 1000.0};
		Vec3 local_max = {-1000.0, -1000.0, -1000.0};

		foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = 
					value<Vec3>(*control_cage_, 
								control_cage_vertex_position_, cv);
				
				const Vec3& local_point = local_frame_*cage_point; 
				for (size_t j = 0; j < 3; j++)
				{
					if (local_point[j] < local_min[j])
					{
						local_min[j] = local_point[j];
					}

					if (local_point[j] > local_max[j])
					{
						local_max[j] = local_point[j];
					}
				}
				return true;
			});

			
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
		update_local_bounding_box(); 

		if (need_full_bind_)
		{
			update_virtual_cages(); 

			update_bind_object(object, object_vertex_position, object_vertex_index); 

			bind_connecting_tools(handle_container, axis_container, cage_container); 
			need_full_bind_ = false;
		}

		if (deformation_type_ == "MVC")
		{
			deform_object_mvc(object, object_vertex_position, 
								object_vertex_index);
		}

		if (deformation_type_ == "Green")
		{
			deform_object_green(object, object_vertex_position, 
								object_vertex_index);
		}
	}

	// @brief initialize binding for handle tool 
	/// @param graph_name name of the handle 
	/// @param handle_position position of the handle 
	void init_bind_handle(const std::string& graph_name, 
						const Vec3& handle_position)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		Handle_data handle_data; 

		VectorWeights handle_weights; 
		 
		handle_weights.position_.resize(nbv_cage);
		handle_weights.position_.setZero();

		handle_data.weights_ = handle_weights; 

		local_handle_data_[graph_name] = handle_data;

		if (deformation_type_ == "MVC"){
			bind_handle_mvc(graph_name, handle_position); 
		}

	}

	/// @brief update binding for handle tool 
	/// used when deformation type is changed
	/// @param graph_name handle name 
	/// @param handle_position handle position 
	void update_bind_handle(const std::string& graph_name, 
					const Vec3& handle_position)
	{
		local_handle_data_[graph_name].weights_.position_.setZero();

		if (deformation_type_ == "MVC"){
			update_bind_handle_mvc(graph_name, handle_position); 
		}
	}

	/// @brief deform the handle 
	/// @param g handle 
	/// @param graph_name name of the handle  
	/// @param graph_vertex_position position of the handle's vertex
	/// @param handle_vertex handle vertex 
	void deform_handle(Graph& g, const std::string& graph_name, 
			std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position, 
					const GraphVertex& handle_vertex)
	{
		if (deformation_type_ == "MVC")
		{
			deform_handle_mvc(g, graph_name, graph_vertex_position, 
								handle_vertex); 
		}

	}


	void init_bind_local_cage(MESH& local_cage, 
							const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		uint32 nbv_target_cage = nb_cells<Vertex>(local_cage);
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		Tool_data cage_data; 

		MatrixWeights local_cage_weights; 
		local_cage_weights.position_.resize(nbv_target_cage, nbv_cage);
		local_cage_weights.position_.setZero();

		cage_data.weights_ = local_cage_weights; 

		local_cage_data_[cage_name] = cage_data;

		bind_local_cage_mvc(local_cage, cage_name, local_cage_vertex_position, local_cage_vertex_index);


	}

	/// @brief update binding for local cage tool 
	/// used when deformation type is changed
	void update_bind_local_cage(MESH& local_cage, const std::string& cage_name, 
	std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
	std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		local_cage_data_[cage_name].weights_.position_.setZero();

		update_bind_local_cage_mvc(local_cage, cage_name, local_cage_vertex_position, local_cage_vertex_index); 
	}


	void deform_local_cage(MESH& local_cage, const std::string& cage_name, 
	std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
	std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		
		deform_local_cage_mvc(local_cage, cage_name, 
						local_cage_vertex_position, local_cage_vertex_index); 

	}

	/// @brief update the control cage dimensions 
	/// call when the content spatial tool deformation goes outside the cage
	/// @param bb_min 
	/// @param bb_max 
	void update_control_cage()
	{
		Eigen::Matrix3d inverse_frame = local_frame_.inverse(); 

		const Vec3 cage_bb_min = inverse_frame*cage_local_bb_min_; 
		const Vec3 cage_bb_max = inverse_frame*cage_local_bb_max_;

		modeling::update_bounding_box(*control_cage_, 
									control_cage_vertex_position_.get(), 
									control_cage_vertices_, 
									cage_bb_min, cage_bb_max);
	}


private:
	std::vector<CMap2::Vertex> control_cage_vertices_; 

	std::vector<Triangle> cage_triangles_;

	Vec3 control_cage_center_;
	Vec3 control_cage_bb_min_;
	Vec3 control_cage_bb_max_;

	Vec3 influence_cage_bb_min_;
	Vec3 influence_cage_bb_max_;

	Local_direction_control_planes local_x_direction_control_planes_;
	Local_direction_control_planes local_y_direction_control_planes_;
	Local_direction_control_planes local_z_direction_control_planes_;

	Eigen::VectorXd control_area_validity_;

	std::vector<Virtual_cube> virtual_cubes_; 	

	Eigen::Matrix<Vec2, Eigen::Dynamic, Eigen::Dynamic> normal_weights_;

	std::vector<Vec3> start_positions_;

	bool need_update_weights_;


	/// @brief update virtual cages after local cage deformation
	void update_virtual_cages_strip()
	{
		for (size_t c = 0; c < virtual_cubes_.size(); c++){
			for (size_t p = 0; p < virtual_cubes_[c].points.size(); p++){
				Point local_point = virtual_cubes_[c].points[p];    

				virtual_cubes_[c].points[p].position = value<Vec3>(*control_cage_, 
						control_cage_vertex_position_, local_point.vertex);

				if (!local_point.inside_control_cage){  
					
					virtual_cubes_[c].points[p].position += 
					(local_point.shift_vector[0]*
					local_x_direction_control_planes_.shift_after_d_max); 

					virtual_cubes_[c].points[p].position += 
					(local_point.shift_vector[1]*
					local_y_direction_control_planes_.shift_after_d_max);

					virtual_cubes_[c].points[p].position += 
					(local_point.shift_vector[2]*
					local_z_direction_control_planes_.shift_after_d_max);

				} 
			}
		}
	}

	/// @brief update virtual cages after local cage deformation
	/// update only local cage points
	void update_virtual_cages_attenuation()
	{
		for (size_t c = 0; c < virtual_cubes_.size(); c++){
			for (size_t p = 0; p < virtual_cubes_[c].points.size(); p++){
				Point local_point = virtual_cubes_[c].points[p]; 

				if (local_point.inside_control_cage)
				{
					virtual_cubes_[c].points[p].position = value<Vec3>(*control_cage_, 
						control_cage_vertex_position_, local_point.vertex);
				}   
			}
		}
	}


	/// @brief set start positions to reset the deformation
	///
	void set_start_positions()
	{
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		start_positions_.resize(nbv_cage);

		foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
			const Vec3& cage_point = value<Vec3>(*control_cage_, 
										control_cage_vertex_position_, cv);

			uint32 cage_point_index = value<uint32>(*control_cage_, 
											control_cage_vertex_index_, cv);

			start_positions_[cage_point_index] = cage_point; 

			return true;
		});
	}

	/// @brief initialize the set of triangles
	/// from the square faces of the control cage 
	void init_triangles()
	{
		foreach_cell(*control_cage_, [&](Face fc) -> bool {
			std::vector<CMap2::Vertex> face_vertices_ = 
								incident_vertices(*control_cage_, fc);

			std::vector<Point> points(4);
			for (std::size_t p = 0; p < 4; p++)
			{
				Point point_i;
				point_i.vertex = face_vertices_[p];
				point_i.position = value<Vec3>(*control_cage_, 
								control_cage_vertex_position_, point_i.vertex);
				point_i.control_cage_index = 
					value<uint32_t>(*control_cage_, 
									control_cage_vertex_index_, point_i.vertex);
				point_i.inside_control_cage = true;
				point_i.shift_vector = {0.0, 0.0, 0.0}; 
				points[p] = point_i; 
			}

			Triangle triangle1;
			triangle1.points = {points[1], points[3], points[0]};
			triangle1.normal = get_normal_from_triangle(triangle1);  
			triangle1.edges = get_edges_from_triangle(triangle1); 
			
			std::make_pair(triangle1.points[1].position 
												- triangle1.points[0].position, 
											triangle1.points[2].position 
												- triangle1.points[1].position);
			cage_triangles_.push_back(triangle1);

			Triangle triangle2;
			triangle2.points = {points[1], points[2], points[3]};
			triangle2.normal = get_normal_from_triangle(triangle2); 
			triangle2.edges = get_edges_from_triangle(triangle2); 
			cage_triangles_.push_back(triangle2);

			return true;
		});
	}

	/// @brief initialize the control cage planes
	/// delimit the control cage area in terms of planes
	/// Plane of equation ax + by + cz = d
	/// Create the three local direction control planes 
	void init_control_cage_plane(const std::tuple<Vec3,Vec3,Vec3>& main_directions, 
									const double& scale)
	{
		const Vec3 x_dir = std::get<0>(main_directions), 
		y_dir = std::get<1>(main_directions), 
		z_dir = std::get<2>(main_directions);

		double d_x_min = 1000.0, d_x_max = -1000.0, 
		d_y_min = 1000.0, d_y_max = -1000.0, 
		d_z_min = 1000.0, d_z_max = -1000.0;

		std::pair<Triangle, Triangle> face_x_min, face_x_max, 
		face_y_min, face_y_max, face_z_min, face_z_max;

		for (std::size_t i = 0; i < 6; i++)
		{ 
			const Triangle triangle1 = cage_triangles_[2 * i];
			const Triangle triangle2 = cage_triangles_[2 * i + 1];

			const Vec3 normal = triangle1.normal;

			const Vec3 triangle_point = triangle1.points[0].position; 

			const double x_projection = normal.dot(x_dir), 
			y_projection = normal.dot(y_dir), 
			z_projection = normal.dot(z_dir); 

			if (abs(x_projection) > 0.5)
			{ 
				const double d = triangle_point.dot(x_dir);
				if (d < d_x_min)
				{
					d_x_min = d;
					face_x_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_x_max)
				{
					d_x_max = d;
					face_x_max = std::make_pair(triangle1, triangle2);
				}
			}
			else if (abs(y_projection) > 0.5)
			{ 
				const double d = triangle_point.dot(y_dir);
				if (d < d_y_min)
				{
					d_y_min = d;
					face_y_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_y_max)
				{
					d_y_max = d;
					face_y_max = std::make_pair(triangle1, triangle2);
				}
			}
			else if (abs(z_projection) > 0.5)
			{
				const double d = triangle_point.dot(z_dir);
				if (d < d_z_min)
				{
					d_z_min = d;
					face_z_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_z_max)
				{
					d_z_max = d;
					face_z_max = std::make_pair(triangle1, triangle2);
				}
			}	
		}

		const double gap_x = scale*(d_x_max - d_x_min); 
		local_x_direction_control_planes_.d_min = d_x_min;
		local_x_direction_control_planes_.d_max = d_x_max;
		local_x_direction_control_planes_.d_gap = gap_x;
		local_x_direction_control_planes_.direction = x_dir;
		local_x_direction_control_planes_.triangles_d_min = face_x_min;
		local_x_direction_control_planes_.triangles_d_max = face_x_max;
		local_x_direction_control_planes_.shift_after_d_max = (gap_x)*x_dir;
		local_x_direction_control_planes_.shift_before_d_min = (-gap_x)*x_dir; 

		const double gap_y = scale*(d_y_max - d_y_min);
		local_y_direction_control_planes_.d_min = d_y_min;
		local_y_direction_control_planes_.d_max = d_y_max;
		local_y_direction_control_planes_.d_gap = gap_y;
		local_y_direction_control_planes_.direction = y_dir;
		local_y_direction_control_planes_.triangles_d_min = face_y_min;
		local_y_direction_control_planes_.triangles_d_max = face_y_max;
		local_y_direction_control_planes_.shift_after_d_max = (gap_y)*y_dir;
		local_y_direction_control_planes_.shift_before_d_min = (-gap_y)*y_dir;

		const double gap_z = scale*(d_z_max - d_z_min); 
		local_z_direction_control_planes_.d_min = d_z_min;
		local_z_direction_control_planes_.d_max = d_z_max;
		local_z_direction_control_planes_.d_gap = gap_z;
		local_z_direction_control_planes_.direction = z_dir;
		local_z_direction_control_planes_.triangles_d_min = face_z_min;
		local_z_direction_control_planes_.triangles_d_max = face_z_max;
		local_z_direction_control_planes_.shift_after_d_max = (gap_z)*z_dir;
		local_z_direction_control_planes_.shift_before_d_min = (-gap_z)*z_dir;

	}

	/// @brief initialize the virtual cubes 
	void init_virtual_cubes()
	{
		init_face_adjacent_virtual_cubes();

		init_edge_adjacent_virtual_cubes();

		init_vertex_adjacent_virtual_cubes();
	}

	/// @brief initialize the virtual cubes 
	/// that share a face with the control cage 
	void init_face_adjacent_virtual_cubes()
	{
		Virtual_cube face_adjacent0 = get_virtual_cube_triangles(
				local_x_direction_control_planes_.triangles_d_min,
				local_x_direction_control_planes_.shift_before_d_min,
				{-1.0, 0.0, 0.0});

		virtual_cubes_[0] = face_adjacent0;	 

		Virtual_cube face_adjacent1 = get_virtual_cube_triangles(
				local_x_direction_control_planes_.triangles_d_max,
				local_x_direction_control_planes_.shift_after_d_max,
				{1.0, 0.0, 0.0});

		virtual_cubes_[1] = face_adjacent1;

		Virtual_cube face_adjacent2 = get_virtual_cube_triangles(
				local_y_direction_control_planes_.triangles_d_min,
				local_y_direction_control_planes_.shift_before_d_min,
				{0.0, -1.0, 0.0});

		virtual_cubes_[2] = face_adjacent2;

		Virtual_cube face_adjacent3 = get_virtual_cube_triangles(
				local_y_direction_control_planes_.triangles_d_max,
				local_y_direction_control_planes_.shift_after_d_max,
				{0.0, 1.0, 0.0});

		virtual_cubes_[3] = face_adjacent3;

		Virtual_cube face_adjacent4 = get_virtual_cube_triangles(
				local_z_direction_control_planes_.triangles_d_min,
				local_z_direction_control_planes_.shift_before_d_min, 
				{0.0, 0.0, -1.0});

		virtual_cubes_[4] = face_adjacent4;

		Virtual_cube face_adjacent5 = get_virtual_cube_triangles(
				local_z_direction_control_planes_.triangles_d_max,
				local_z_direction_control_planes_.shift_after_d_max, 
				{0.0, 0.0, 1.0});

		virtual_cubes_[5] = face_adjacent5;
	}

	/// @brief initialize the virtual cubes 
	/// that share an edge with the control cage 
	void init_edge_adjacent_virtual_cubes()
	{
		const Vec3 shift_x_min = 
				local_x_direction_control_planes_.shift_before_d_min;
		const Vec3 shift_x_max = 
				local_x_direction_control_planes_.shift_after_d_max;
		
		const Vec3 unit_x_min = {-1.0, 0.0, 0.0}, 
		unit_x_max = {1.0, 0.0, 0.0}; 

		const Vec3 shift_y_min = 
				local_y_direction_control_planes_.shift_before_d_min;
		const Vec3 shift_y_max = 
				local_y_direction_control_planes_.shift_after_d_max;

		const Vec3 unit_y_min = {0.0, -1.0, 0.0}, 
		unit_y_max = {0.0, 1.0, 0.0};

		const Vec3 shift_z_min = 
				local_z_direction_control_planes_.shift_before_d_min;
		const Vec3 shift_z_max = 
				local_z_direction_control_planes_.shift_after_d_max;

		const Vec3 unit_z_min = {0.0, 0.0, -1.0}, 
		unit_z_max = {0.0, 0.0, 1.0};

		std::vector<Point> intersect_points0 = 
			find_intersection_points_face(
					local_x_direction_control_planes_.triangles_d_min, 
					local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face0 = 
			get_face_from_intersecting_edge(intersect_points0, shift_x_min,
											unit_x_min);

		Virtual_cube edge_adjacent0 = get_virtual_cube_triangles(new_face0, 
													shift_z_min, unit_z_min);

		virtual_cubes_[6] = edge_adjacent0;

		std::vector<Point> intersect_points1 = 
			find_intersection_points_face(
				local_x_direction_control_planes_.triangles_d_max, 
				local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face1 = 
			get_face_from_intersecting_edge(intersect_points1, shift_x_max, 
											unit_x_max);

		Virtual_cube edge_adjacent1 = get_virtual_cube_triangles(new_face1, 
													shift_z_min, unit_z_min);

		virtual_cubes_[7] = edge_adjacent1;

		std::vector<Point> intersect_points2 = 
			find_intersection_points_face(
				local_x_direction_control_planes_.triangles_d_min, 
				local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face2 = 
			get_face_from_intersecting_edge(intersect_points2, shift_x_min, 
											unit_x_min);

		Virtual_cube edge_adjacent2 = get_virtual_cube_triangles(new_face2, 
													shift_z_max, unit_z_max);
		virtual_cubes_[8] = edge_adjacent2;

		std::vector<Point> intersect_points3 = 
			find_intersection_points_face(
				local_x_direction_control_planes_.triangles_d_max, 
				local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face3 = 
			get_face_from_intersecting_edge(intersect_points3, shift_x_max,
											unit_x_max);

		Virtual_cube edge_adjacent3 = get_virtual_cube_triangles(new_face3, 
													shift_z_max, unit_z_max);

		virtual_cubes_[9] = edge_adjacent3;

		std::vector<Point> intersect_points4 = 
			find_intersection_points_face(
				local_y_direction_control_planes_.triangles_d_min, 
				local_x_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face4 = 
			get_face_from_intersecting_edge(intersect_points4, shift_y_min,
											unit_y_min);

		Virtual_cube edge_adjacent4 = get_virtual_cube_triangles(new_face4, 
													shift_x_min, unit_x_min);

		virtual_cubes_[10] = edge_adjacent4;

		std::vector<Point> intersect_points5 = 
			find_intersection_points_face(
				local_y_direction_control_planes_.triangles_d_max, 
				local_x_direction_control_planes_.triangles_d_min);
		std::pair<Triangle, Triangle> new_face5 = 
			get_face_from_intersecting_edge(intersect_points5, shift_y_max,
											unit_y_max);

		Virtual_cube edge_adjacent5 = get_virtual_cube_triangles(new_face5, 
													shift_x_min, unit_x_min);

		virtual_cubes_[11] = edge_adjacent5;

		std::vector<Point> intersect_points6 = 
			find_intersection_points_face(
				local_y_direction_control_planes_.triangles_d_min, 
				local_x_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face6 = 
			get_face_from_intersecting_edge(intersect_points6, shift_y_min,
											unit_y_min);

		Virtual_cube edge_adjacent6 = get_virtual_cube_triangles(new_face6, 
													shift_x_max, unit_x_max);
		virtual_cubes_[12] = edge_adjacent6;

		std::vector<Point> intersect_points7 = 
			find_intersection_points_face(
				local_y_direction_control_planes_.triangles_d_max, 
				local_x_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face7 = 
			get_face_from_intersecting_edge(intersect_points7, shift_y_max,
											unit_y_max);

		Virtual_cube edge_adjacent7 = get_virtual_cube_triangles(new_face7, 
													shift_x_max, unit_x_max);

		virtual_cubes_[13] = edge_adjacent7;

		std::vector<Point> intersect_points8 = 
			find_intersection_points_face(
				local_y_direction_control_planes_.triangles_d_min, 
				local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face8 = 
			get_face_from_intersecting_edge(intersect_points8, shift_y_min,
											unit_y_min);

		Virtual_cube edge_adjacent8 = get_virtual_cube_triangles(new_face8, 
													shift_z_min, unit_z_min);
		virtual_cubes_[14] = edge_adjacent8;

		std::vector<Point> intersect_points9 = 
			find_intersection_points_face(
				local_y_direction_control_planes_.triangles_d_max, 
				local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face9 = 
			get_face_from_intersecting_edge(intersect_points9, shift_y_max,
											unit_y_max);

		Virtual_cube edge_adjacent9 = get_virtual_cube_triangles(new_face9, 
													shift_z_min, unit_z_min);
		virtual_cubes_[15] = edge_adjacent9;

		std::vector<Point> intersect_points10 = 
			find_intersection_points_face(
				local_y_direction_control_planes_.triangles_d_min, 
				local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face10 = 
			get_face_from_intersecting_edge(intersect_points10, shift_y_min,
											unit_y_min);

		Virtual_cube edge_adjacent10 = get_virtual_cube_triangles(new_face10, 
													shift_z_max, unit_z_max);
		virtual_cubes_[16] = edge_adjacent10;

		std::vector<Point> intersect_points11 = 
			find_intersection_points_face(
				local_y_direction_control_planes_.triangles_d_max, 
				local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face11 = 
			get_face_from_intersecting_edge(intersect_points11, shift_y_max,
											unit_y_max);

		Virtual_cube edge_adjacent11 = get_virtual_cube_triangles(new_face11, 
													shift_z_max, unit_z_max);
		virtual_cubes_[17] = edge_adjacent11;
	}


	/// @brief initialize the virtual cubes 
	/// that share a vertex with the control cage 
	void init_vertex_adjacent_virtual_cubes()
	{
		const Vec3 shift_x_min = 
			local_x_direction_control_planes_.shift_before_d_min;
		const Vec3 shift_x_max = 
			local_x_direction_control_planes_.shift_after_d_max;

		const Vec3 unit_x_min = {-1.0, 0.0, 0.0}, 
		unit_x_max = {1.0, 0.0, 0.0};

		const Vec3 shift_y_min = 
			local_y_direction_control_planes_.shift_before_d_min;
		const Vec3 shift_y_max = 
			local_y_direction_control_planes_.shift_after_d_max;

		const Vec3 unit_y_min = {0.0, -1.0, 0.0}, 
		unit_y_max = {0.0, 1.0, 0.0};

		const Vec3 shift_z_min = 
			local_z_direction_control_planes_.shift_before_d_min;
		const Vec3 shift_z_max = 
			local_z_direction_control_planes_.shift_after_d_max;

		const Vec3 unit_z_min = {0.0, 0.0, -1.0}, 
		unit_z_max = {0.0, 0.0, 1.0};

		Point intersection_point0 = 
			find_intersection_point(
					local_x_direction_control_planes_.triangles_d_min,
					local_y_direction_control_planes_.triangles_d_min,
					local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face0 =
			get_face_from_intersecting_vertex(intersection_point0, 
											shift_x_min, shift_y_min,
											unit_x_min, unit_y_min);

		Virtual_cube vertex_adjacent0 = get_virtual_cube_triangles(new_face0, 
													shift_z_min, unit_z_min);
		virtual_cubes_[18] = vertex_adjacent0;

		Point intersection_point1 = 
			find_intersection_point(
				local_x_direction_control_planes_.triangles_d_max,
				local_y_direction_control_planes_.triangles_d_min,
				local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face1 =
			get_face_from_intersecting_vertex(intersection_point1, 
											shift_x_max, shift_y_min,
											unit_x_max, unit_y_min);

		Virtual_cube vertex_adjacent1 = get_virtual_cube_triangles(new_face1, 
													shift_z_min, unit_z_min);
		virtual_cubes_[19] = vertex_adjacent1;

		Point intersection_point2 = 
			find_intersection_point(
				local_x_direction_control_planes_.triangles_d_min,
				local_y_direction_control_planes_.triangles_d_min,
				local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face2 =
			get_face_from_intersecting_vertex(intersection_point2, 
											shift_x_min, shift_y_min,
											unit_x_min, unit_y_min);

		Virtual_cube vertex_adjacent2 = get_virtual_cube_triangles(new_face2, 
													shift_z_max, unit_z_max);
		virtual_cubes_[20] = vertex_adjacent2;

		Point intersection_point3 = 
			find_intersection_point(
				local_x_direction_control_planes_.triangles_d_max,
				local_y_direction_control_planes_.triangles_d_min,
				local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face3 =
			get_face_from_intersecting_vertex(intersection_point3, 
											shift_x_max, shift_y_min,
											unit_x_max, unit_y_min);
		
		Virtual_cube vertex_adjacent3 = get_virtual_cube_triangles(new_face3, 
													shift_z_max, unit_z_max);
		virtual_cubes_[21] = vertex_adjacent3;

		Point intersection_point4 = 
			find_intersection_point(
				local_x_direction_control_planes_.triangles_d_min,
				local_y_direction_control_planes_.triangles_d_max,
				local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face4 =
			get_face_from_intersecting_vertex(intersection_point4, 
											shift_x_min, shift_y_max,
											unit_x_min, unit_y_max);
		
		Virtual_cube vertex_adjacent4 = get_virtual_cube_triangles(new_face4, 
													shift_z_min, unit_z_min);
		virtual_cubes_[22] = vertex_adjacent4;

		Point intersection_point5 = 
			find_intersection_point(
				local_x_direction_control_planes_.triangles_d_max,
				local_y_direction_control_planes_.triangles_d_max,
				local_z_direction_control_planes_.triangles_d_min);

		std::pair<Triangle, Triangle> new_face5 =
			get_face_from_intersecting_vertex(intersection_point5, 
											shift_x_max, shift_y_max,
											unit_x_max, unit_y_max);
		
		Virtual_cube vertex_adjacent5 = get_virtual_cube_triangles(new_face5, 
													shift_z_min, unit_z_min);
		virtual_cubes_[23] = vertex_adjacent5;

		Point intersection_point6 = 
			find_intersection_point(
				local_x_direction_control_planes_.triangles_d_min,
				local_y_direction_control_planes_.triangles_d_max,
				local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face6 =
			get_face_from_intersecting_vertex(intersection_point6, 
											shift_x_min, shift_y_max,
											unit_x_min, unit_y_max);
		
		Virtual_cube vertex_adjacent6 = get_virtual_cube_triangles(new_face6, 
													shift_z_max, unit_z_max);

		virtual_cubes_[24] = vertex_adjacent6;

		Point intersection_point7 = 
			find_intersection_point(
				local_x_direction_control_planes_.triangles_d_max,
				local_y_direction_control_planes_.triangles_d_max,
				local_z_direction_control_planes_.triangles_d_max);

		std::pair<Triangle, Triangle> new_face7 =
			get_face_from_intersecting_vertex(intersection_point7, 
											shift_x_max, shift_y_max,
											unit_x_max, unit_y_max);
		
		Virtual_cube vertex_adjacent7 = get_virtual_cube_triangles(new_face7, 
													shift_z_max, unit_z_max);
		virtual_cubes_[25] = vertex_adjacent7;
	}

	/// @brief compute a virtual cube from a face and a shift vector
	/// create the opposite face by shifting the points of the face by the vector 
	/// link the set of points to create the virtual cube from it 
	/// @param face pair of triangles 
	/// @param shift_vector 
	/// @return virtual cube 
	Virtual_cube get_virtual_cube_triangles(
								const std::pair<Triangle, Triangle> face, 
								const Vec3& shift_vector, const Vec3& unit_dir)
	{
		Virtual_cube new_virtual_cube;

		std::vector<Triangle> virtual_cube_triangles;

		Triangle triangle1 = face.first;
		triangle1.virtual_cage_indices = {3, 1, 0};
		Triangle triangle2 = face.second;
		triangle2.virtual_cage_indices = {3, 2, 1};
		virtual_cube_triangles.push_back(triangle1);
		virtual_cube_triangles.push_back(triangle2); 

		Point point4;
		point4.position = triangle1.points[2].position + shift_vector;
		point4.inside_control_cage = false;
		point4.shift_vector = {0.0, 0.0, 0.0}; 
		point4.shift_vector += triangle1.points[2].shift_vector; 
		point4.shift_vector += unit_dir; 
		point4.vertex = triangle1.points[2].vertex; 

		Point point5;
		point5.position = triangle1.points[1].position + shift_vector;
		point5.inside_control_cage = false;
		point5.shift_vector = {0.0, 0.0, 0.0}; 
		point5.shift_vector += triangle1.points[1].shift_vector; 
		point5.shift_vector += unit_dir; 
		point5.vertex = triangle1.points[1].vertex;

		Point point6;
		point6.position = triangle1.points[0].position + shift_vector;
		point6.inside_control_cage = false;
		point6.shift_vector = {0.0, 0.0, 0.0}; 
		point6.shift_vector += triangle1.points[0].shift_vector; 
		point6.shift_vector += unit_dir; 
		point6.vertex = triangle1.points[0].vertex; 

		Point point7;
		point7.position = triangle2.points[1].position + shift_vector;
		point7.inside_control_cage = false;
		point7.shift_vector = {0.0, 0.0, 0.0}; 
		point7.shift_vector += triangle2.points[1].shift_vector; 
		point7.shift_vector += unit_dir; 
		point7.vertex = triangle2.points[1].vertex;  

		Triangle triangle3, triangle4;
		triangle3.points = {point4, point5, point7};
		triangle3.virtual_cage_indices = {4, 5, 7};
		triangle3.normal = get_normal_from_triangle(triangle3);
		triangle3.edges = get_edges_from_triangle(triangle3);

		triangle4.points = {point4, point7, point6};
		triangle4.virtual_cage_indices = {4, 7, 6};
		triangle4.normal = get_normal_from_triangle(triangle4);
		triangle4.edges = get_edges_from_triangle(triangle4);

		std::vector<Point> virtual_cube_points(8);
		virtual_cube_points[0] = triangle1.points[2];
		virtual_cube_points[1] = triangle1.points[1];
		virtual_cube_points[2] = triangle2.points[1];
		virtual_cube_points[3] = triangle1.points[0];

		virtual_cube_points[4] = triangle3.points[0];
		virtual_cube_points[5] = triangle3.points[1];
		virtual_cube_points[6] = triangle4.points[2];
		virtual_cube_points[7] = triangle3.points[2];

		new_virtual_cube.points = virtual_cube_points;

		Triangle triangle5, triangle6, triangle7, triangle8, 
		triangle9, triangle10, triangle11, triangle12;

		triangle5.points = 
			{triangle1.points[2], triangle3.points[1], triangle3.points[0]};
		triangle5.virtual_cage_indices = {0, 5, 4};
		triangle5.normal = get_normal_from_triangle(triangle5);
		triangle5.edges = get_edges_from_triangle(triangle5);

		triangle6.points = 
			{triangle1.points[2], triangle1.points[1], triangle3.points[1]};
		triangle6.virtual_cage_indices = {0, 1, 5};
		triangle6.normal = get_normal_from_triangle(triangle6);
		triangle6.edges = get_edges_from_triangle(triangle6);

		triangle7.points = 
			{triangle4.points[2], triangle2.points[1], triangle1.points[0]};
		triangle7.virtual_cage_indices = {6, 2, 3};
		triangle7.normal = get_normal_from_triangle(triangle7);
		triangle7.edges = get_edges_from_triangle(triangle7);

		triangle8.points = 
			{triangle4.points[2], triangle3.points[2], triangle2.points[1]};
		triangle8.virtual_cage_indices = {6, 7, 2};
		triangle8.normal = get_normal_from_triangle(triangle8);
		triangle8.edges = get_edges_from_triangle(triangle8);

		triangle9.points = 
			{triangle1.points[2], triangle4.points[2], triangle1.points[0]};
		triangle9.virtual_cage_indices = {0, 6, 3};
		triangle9.normal = get_normal_from_triangle(triangle9);
		triangle9.edges = get_edges_from_triangle(triangle9);

		triangle10.points = 
			{triangle1.points[2], triangle3.points[0], triangle4.points[2]};
		triangle10.virtual_cage_indices = {0, 4, 6};
		triangle10.normal = get_normal_from_triangle(triangle10);
		triangle10.edges = get_edges_from_triangle(triangle10);

		triangle11.points = 
			{triangle2.points[1], triangle3.points[1], triangle1.points[1]};
		triangle11.virtual_cage_indices = {2, 5, 1};
		triangle11.normal = get_normal_from_triangle(triangle11);
		triangle11.edges = get_edges_from_triangle(triangle11);

		triangle12.points = 
			{triangle2.points[1], triangle3.points[2], triangle3.points[1]};
		triangle12.virtual_cage_indices = {2, 7, 5};
		triangle12.normal = get_normal_from_triangle(triangle12);
		triangle12.edges = get_edges_from_triangle(triangle12);

		virtual_cube_triangles.push_back(triangle3);
		virtual_cube_triangles.push_back(triangle4);
		virtual_cube_triangles.push_back(triangle5);
		virtual_cube_triangles.push_back(triangle6);
		virtual_cube_triangles.push_back(triangle7);
		virtual_cube_triangles.push_back(triangle8);
		virtual_cube_triangles.push_back(triangle9);
		virtual_cube_triangles.push_back(triangle10);
		virtual_cube_triangles.push_back(triangle11);
		virtual_cube_triangles.push_back(triangle12);

		new_virtual_cube.triangles = virtual_cube_triangles;

		return new_virtual_cube;
	}


	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model 
	/// @param object_vertex_index index of the vertices of the model 
	void update_bind_object_mvc(MESH& object, 
			CMap2::Attribute<Vec3>* object_vertex_position, 
			CMap2::Attribute<uint32>* object_vertex_index)
	{

		foreach_cell(object, [&](Vertex v) -> bool {

			const Vec3& surface_point_position = 
							value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_index = 
							value<uint32>(object, object_vertex_index, v);

			if (object_activation_cage_[surface_point_index] == 27){
				Eigen::VectorXd local_row_weights; 
				compute_mvc_on_point_inside_cage(surface_point_position, 
												local_row_weights);

				object_weights_.position_.row(surface_point_index) = local_row_weights; 
			} else {
				Virtual_cube virtual_cube_target = virtual_cubes_[object_activation_cage_[surface_point_index]]; 

				Eigen::VectorXd local_row_weights; 
				compute_mvc_on_point_outside_cage(surface_point_position, 
												local_row_weights, 
												virtual_cube_target);

				object_weights_.position_.row(surface_point_index) = local_row_weights; 
			}

			return true; 
		}); 
	}



	/// @brief bind object with MVC deformation type
	/// loop on each point of the influence area 
	/// check if point inside control cage => classical MVC binding 
	/// otherwise: find corresponding virtual cube, compute MVC inside it
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model 
	/// @param object_vertex_index index of the vertices of the model 
	void bind_object_mvc(MESH& object, 
			CMap2::Attribute<Vec3>* object_vertex_position, 
			CMap2::Attribute<uint32>* object_vertex_index)
	{

		foreach_cell(object, [&](Vertex v) -> bool {

			const Vec3& surface_point_position = 
							value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_index = 
							value<uint32>(object, object_vertex_index, v);

			const double d_x = get_projection_on_direction(surface_point_position, local_x_direction_control_planes_.direction), 
			
			d_y = get_projection_on_direction(surface_point_position, local_y_direction_control_planes_.direction),
		
			d_z = get_projection_on_direction(surface_point_position, local_z_direction_control_planes_.direction); 
					

			const bool valid_x_dir = check_projection_in_area(d_x, local_x_direction_control_planes_.d_min, local_x_direction_control_planes_.d_max), 
					
			valid_y_dir = check_projection_in_area(d_y, local_y_direction_control_planes_.d_min, local_y_direction_control_planes_.d_max), 
					   
			valid_z_dir = check_projection_in_area(d_z, local_z_direction_control_planes_.d_min, local_z_direction_control_planes_.d_max); 

			if (valid_x_dir && valid_y_dir && valid_z_dir)
			{
				object_activation_cage_[surface_point_index] = 27; 
				Eigen::VectorXd local_row_weights; 
				compute_mvc_on_point_inside_cage(surface_point_position, 
												local_row_weights);

				object_weights_.position_.row(surface_point_index) = local_row_weights; 
			}
			else
			{
				Vec3 projection_values = {d_x, d_y, d_z};  
				
				std::vector<bool> valid_values(3); 
				valid_values[0] = valid_x_dir, valid_values[1] = valid_y_dir, 
				valid_values[2] = valid_z_dir; 

				Vec3 max_area = {local_x_direction_control_planes_.d_max, local_y_direction_control_planes_.d_max, local_z_direction_control_planes_.d_max}; 

				std::size_t virtual_cube_index = get_index_virtual_cube(projection_values, valid_values, max_area); 

				object_activation_cage_[surface_point_index] = virtual_cube_index; 

				Virtual_cube virtual_cube_target = virtual_cubes_[virtual_cube_index]; 

				Eigen::VectorXd local_row_weights; 
				compute_mvc_on_point_outside_cage(surface_point_position, 
												local_row_weights, 
												virtual_cube_target);

				object_weights_.position_.row(surface_point_index) = local_row_weights; 
			}
			return true; 
		}); 
	}


	/// @brief bind object with MVC deformation type
	/// loop on each point of the influence area 
	/// check if point inside control cage => classical MVC binding 
	/// otherwise: find corresponding virtual cube, compute MVC inside it
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model 
	/// @param object_vertex_index index of the vertices of the model 
	void bind_object_mvc_parallel(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position, 
						CMap2::Attribute<uint32>* object_vertex_index)
	{
		parallel_foreach_cell(object, [&](Vertex v) -> bool {

			const Vec3& surface_point_position = 
							value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_index = 
							value<uint32>(object, object_vertex_index, v);

			const double d_x = get_projection_on_direction(surface_point_position, local_x_direction_control_planes_.direction), 
			
			d_y = get_projection_on_direction(surface_point_position, local_y_direction_control_planes_.direction),
		
			d_z = get_projection_on_direction(surface_point_position, local_z_direction_control_planes_.direction); 
					

			const bool valid_x_dir = check_projection_in_area(d_x, local_x_direction_control_planes_.d_min, local_x_direction_control_planes_.d_max), 
					
			valid_y_dir = check_projection_in_area(d_y, local_y_direction_control_planes_.d_min, local_y_direction_control_planes_.d_max), 
					   
			valid_z_dir = check_projection_in_area(d_z, local_z_direction_control_planes_.d_min, local_z_direction_control_planes_.d_max); 

			if (valid_x_dir && valid_y_dir && valid_z_dir)
			{
				object_activation_cage_[surface_point_index] = 27; 
				Eigen::VectorXd local_row_weights; 
				compute_mvc_on_point_inside_cage(surface_point_position, 
												local_row_weights);

				object_weights_.position_.row(surface_point_index) = local_row_weights; 
			}
			else
			{
				Vec3 projection_values = {d_x, d_y, d_z};  
				
				std::vector<bool> valid_values(3); 
				valid_values[0] = valid_x_dir, valid_values[1] = valid_y_dir, 
				valid_values[2] = valid_z_dir; 

				Vec3 max_area = {local_x_direction_control_planes_.d_max, local_y_direction_control_planes_.d_max, local_z_direction_control_planes_.d_max}; 

				std::size_t virtual_cube_index = get_index_virtual_cube(projection_values, valid_values, max_area); 

				object_activation_cage_[surface_point_index] = virtual_cube_index; 

				Virtual_cube virtual_cube_target = virtual_cubes_[virtual_cube_index]; 

				Eigen::VectorXd local_row_weights; 
				compute_mvc_on_point_outside_cage(surface_point_position, 
												local_row_weights, 
												virtual_cube_target);

				object_weights_.position_.row(surface_point_index) = local_row_weights; 
			}
			return true; 
		}); 
	}

	/// @brief deform the influence area of the object with MVC deformation type
	/// rely object_fixed_data to handle case of point inside virtual cube
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model 
	/// @param object_vertex_index index of the vertices of the model 
	void deform_object_mvc(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position,
						CMap2::Attribute<uint32>* object_vertex_index)
	{

		foreach_cell(object, [&](Vertex v) -> bool {
			uint32 object_point_index = 
							value<uint32>(object, object_vertex_index, v);
			
			if (object_activation_cage_[object_point_index] == 27)
			{
				Vec3 new_position = {0.0, 0.0, 0.0};

				foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = 
					value<Vec3>(*control_cage_, 
								control_cage_vertex_position_, cv);
				uint32 cage_point_index = 
					value<uint32>(*control_cage_, 
								control_cage_vertex_index_, cv);

				new_position += 
					object_weights_.position_(object_point_index, 
											cage_point_index) * cage_point;

				return true;
			});
				value<Vec3>(object, object_vertex_position, v) = new_position;
			
			} else {
				
				Virtual_cube virtual_cube_target = 
						virtual_cubes_[object_activation_cage_[object_point_index]];  

				Vec3 new_position = {0.0, 0.0, 0.0};
				for (size_t p = 0; p < virtual_cube_target.points.size(); p++){
					const Point& cage_point = virtual_cube_target.points[p];
					uint32 cage_point_index = p;

					new_position += 
						object_weights_.position_(object_point_index, p) * 
									cage_point.position;

				}
				value<Vec3>(object, object_vertex_position, v) = new_position;
			 
			}

			return true; 
		}); 
	}

	/// @brief deform the influence area of the object with MVC deformation type
	/// rely object_fixed_data to handle case of point inside virtual cube
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model 
	/// @param object_vertex_index index of the vertices of the model 
	void deform_object_mvc_parallel(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position,
						CMap2::Attribute<uint32>* object_vertex_index)
	{

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			uint32 object_point_index = 
							value<uint32>(object, object_vertex_index, v);
			
			if (object_activation_cage_[object_point_index] == 27)
			{
				Vec3 new_position = {0.0, 0.0, 0.0};

				foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = 
					value<Vec3>(*control_cage_, 
								control_cage_vertex_position_, cv);
				uint32 cage_point_index = 
					value<uint32>(*control_cage_, 
								control_cage_vertex_index_, cv);

				new_position += 
					object_weights_.position_(object_point_index, 
											cage_point_index) * cage_point;

				return true;
			});
				value<Vec3>(object, object_vertex_position, v) = new_position;
			
			} else {
				
				Virtual_cube virtual_cube_target = 
						virtual_cubes_[object_activation_cage_[object_point_index]];  

				Vec3 new_position = {0.0, 0.0, 0.0};
				for (size_t p = 0; p < virtual_cube_target.points.size(); p++){
					const Point& cage_point = virtual_cube_target.points[p];
					uint32 cage_point_index = p;

					new_position += 
						object_weights_.position_(object_point_index, p) * 
									cage_point.position;

				}
				value<Vec3>(object, object_vertex_position, v) = new_position;
			 
			}

			return true; 
		}); 
	}

	void update_bind_handle_mvc(const std::string& graph_name, 
						const Vec3& handle_position)
	{
		if (local_handle_data_[graph_name].cage_index_ ==27 )
		{

			Eigen::VectorXd local_row_weights;
			compute_mvc_on_point_inside_cage(handle_position, local_row_weights);

			local_handle_data_[graph_name].weights_.position_ = local_row_weights; 

		} else {
			Virtual_cube virtual_cube_target = virtual_cubes_[local_handle_data_[graph_name].cage_index_]; 

			Eigen::VectorXd local_row_weights;
			compute_mvc_on_point_outside_cage(handle_position, 			local_row_weights, virtual_cube_target);

			local_handle_data_[graph_name].weights_.position_ = local_row_weights;
		}
	}


	// @brief bind handle with MVC deformation type
	void bind_handle_mvc(const std::string& graph_name, 
						const Vec3& handle_position){
 
		const double d_x = get_projection_on_direction(handle_position, local_x_direction_control_planes_.direction), 
			
		d_y = get_projection_on_direction(handle_position, local_y_direction_control_planes_.direction),
		
			d_z = get_projection_on_direction(handle_position, local_z_direction_control_planes_.direction); 

		const bool valid_x_dir = check_projection_in_area(d_x, local_x_direction_control_planes_.d_min, local_x_direction_control_planes_.d_max), 
					
			valid_y_dir = check_projection_in_area(d_y, local_y_direction_control_planes_.d_min, local_y_direction_control_planes_.d_max), 
					   
			valid_z_dir = check_projection_in_area(d_z, local_z_direction_control_planes_.d_min, local_z_direction_control_planes_.d_max); 

		
		if (valid_x_dir && valid_y_dir && valid_z_dir)
		{
			local_handle_data_[graph_name].cage_index_ = 27; 

			Eigen::VectorXd local_row_weights;
			compute_mvc_on_point_inside_cage(handle_position, local_row_weights);

			local_handle_data_[graph_name].weights_.position_ = local_row_weights; 
		}
		else
		{
			Vec3 projection_values = {d_x, d_y, d_z};  
			
			std::vector<bool> valid_values(3); 
			valid_values[0] = valid_x_dir, valid_values[1] = valid_y_dir, 
			valid_values[2] = valid_z_dir; 

			Vec3 max_area = {local_x_direction_control_planes_.d_max, local_y_direction_control_planes_.d_max, local_z_direction_control_planes_.d_max}; 

			std::size_t virtual_cube_index = get_index_virtual_cube(projection_values, valid_values, max_area); 

			local_handle_data_[graph_name].cage_index_ = virtual_cube_index; 

			Virtual_cube virtual_cube_target = virtual_cubes_[virtual_cube_index]; 

			Eigen::VectorXd local_row_weights;
			compute_mvc_on_point_outside_cage(handle_position, 			local_row_weights, virtual_cube_target);

			local_handle_data_[graph_name].weights_.position_ = local_row_weights;
		} 								
	}

	/// @brief deform handle with MVC deformation type 
	/// @param g handle 
	/// @param graph_name handle name  
	/// @param graph_vertex_position position of the handle vertex
	/// @param handle_vertex handle vertex
	void deform_handle_mvc(Graph& g, const std::string& graph_name, 
			std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
							const GraphVertex& handle_vertex)
	{
		VectorWeights handle_weights = local_handle_data_[graph_name].weights_;

		int handle_cage_index = local_handle_data_[graph_name].cage_index_; 

		if (handle_cage_index == 27)
		{
			Vec3 new_position = {0.0, 0.0, 0.0};

			foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*control_cage_, 
										control_cage_vertex_position_, cv);
				uint32 cage_point_index = value<uint32>(*control_cage_, 
											control_cage_vertex_index_, cv);

				new_position += 
					handle_weights.position_[cage_point_index] * cage_point;

			return true;
			});

			value<Vec3>(g, graph_vertex_position, handle_vertex) = new_position;
		} else {
			Virtual_cube virtual_cube_target =  
				virtual_cubes_[handle_cage_index];  

			Vec3 new_position = {0.0, 0.0, 0.0};
			for (size_t p = 0; p < virtual_cube_target.points.size(); p++){
				const Point& cage_point = virtual_cube_target.points[p];
				uint32 cage_point_index = p;

				new_position += 
					handle_weights.position_[p] * 
									cage_point.position;

				}
				value<Vec3>(g, graph_vertex_position, handle_vertex) = new_position;
		}
		
	}

	void bind_local_cage_mvc(MESH& local_cage, 
							const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{

		foreach_cell(local_cage, [&](Vertex v) -> bool {

			const Vec3& cage_position = value<Vec3>(local_cage, 
											local_cage_vertex_position, v);
			const uint32& cage_point_index = 
					value<uint32>(local_cage, local_cage_vertex_index, v);

			const double d_x = get_projection_on_direction(cage_position, local_x_direction_control_planes_.direction), 
			
			d_y = get_projection_on_direction(cage_position, local_y_direction_control_planes_.direction),
		
			d_z = get_projection_on_direction(cage_position, local_z_direction_control_planes_.direction); 

			const bool valid_x_dir = check_projection_in_area(d_x, local_x_direction_control_planes_.d_min, local_x_direction_control_planes_.d_max), 
					
			valid_y_dir = check_projection_in_area(d_y, local_y_direction_control_planes_.d_min, local_y_direction_control_planes_.d_max), 
					   
			valid_z_dir = check_projection_in_area(d_z, local_z_direction_control_planes_.d_min, local_z_direction_control_planes_.d_max); 
		
		if (valid_x_dir && valid_y_dir && valid_z_dir)
		{
			local_cage_data_[cage_name].cage_index_[cage_point_index] = 27; 

			Eigen::VectorXd local_row_weights;
			compute_mvc_on_point_inside_cage(cage_position, local_row_weights);

			local_cage_data_[cage_name].weights_.position_.row(cage_point_index) = local_row_weights; 

		}
		else
		{

			Vec3 projection_values = {d_x, d_y, d_z};  
			
			std::vector<bool> valid_values(3); 
			valid_values[0] = valid_x_dir, valid_values[1] = valid_y_dir, 
			valid_values[2] = valid_z_dir; 

			Vec3 max_area = {local_x_direction_control_planes_.d_max, local_y_direction_control_planes_.d_max, local_z_direction_control_planes_.d_max}; 

			std::size_t virtual_cube_index = get_index_virtual_cube(projection_values, valid_values, max_area); 

			local_cage_data_[cage_name].cage_index_[cage_point_index] = virtual_cube_index; 

			Virtual_cube virtual_cube_target = virtual_cubes_[virtual_cube_index]; 

			Eigen::VectorXd local_row_weights;
			compute_mvc_on_point_outside_cage(cage_position, 			local_row_weights, virtual_cube_target);

			local_cage_data_[cage_name].weights_.position_.row(cage_point_index) = local_row_weights;

		} 						

			return true; 
		}); 
	}

	void update_bind_local_cage_mvc(MESH& local_cage, 
							const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		foreach_cell(local_cage, [&](Vertex v) -> bool {
			const Vec3& cage_position = value<Vec3>(local_cage, 
											local_cage_vertex_position, v);
			const uint32& cage_point_index = 
					value<uint32>(local_cage, local_cage_vertex_index, v);

			const uint32 activation_index = local_cage_data_[cage_name].cage_index_[cage_point_index]; 

			if (activation_index == 27)
			{
				Eigen::VectorXd local_row_weights;
				compute_mvc_on_point_inside_cage(cage_position, local_row_weights);

				local_cage_data_[cage_name].weights_.position_.row(cage_point_index) = local_row_weights; 

			} 
			else {
				Virtual_cube virtual_cube_target = virtual_cubes_[activation_index]; 

				Eigen::VectorXd local_row_weights;
				compute_mvc_on_point_outside_cage(cage_position, 			local_row_weights, virtual_cube_target);

				local_cage_data_[cage_name].weights_.position_.row(cage_point_index) = local_row_weights;
			}

			return true; 
		}); 
	}

	void deform_local_cage_mvc(MESH& local_cage, 
							const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		MatrixWeights local_cage_weights = local_cage_data_[cage_name].weights_;

		 foreach_cell(local_cage, [&](Vertex v) -> bool {
			uint32 local_cage_point_index = value<uint32>(local_cage, 
												local_cage_vertex_index, v);

			uint32 activation_index = local_cage_data_[cage_name].cage_index_[local_cage_point_index]; 

			if (activation_index == 27)
			{
				Vec3 new_position = {0.0, 0.0, 0.0};

				foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
					const Vec3& cage_point = value<Vec3>(*control_cage_, 
										control_cage_vertex_position_, cv);

					uint32 cage_point_index = value<uint32>(*control_cage_, 
											control_cage_vertex_index_, cv);

					new_position += local_cage_data_[cage_name].weights_
					.position_(local_cage_point_index, cage_point_index) 
						* cage_point;

					return true;
				});

				value<Vec3>(local_cage, local_cage_vertex_position, v) = new_position;
			} 
			else {
				Virtual_cube virtual_cube_target =  virtual_cubes_[activation_index];  

				Vec3 new_position = {0.0, 0.0, 0.0};
				for (size_t p = 0; p < virtual_cube_target.points.size(); p++){
					const Point& cage_point = virtual_cube_target.points[p];
					uint32 cage_point_index = p;

					new_position += local_cage_data_[cage_name].weights_
					.position_(local_cage_point_index, cage_point_index) 
						* cage_point.position;

					}
					value<Vec3>(local_cage, local_cage_vertex_position, v) = new_position;
			}
			return true;
		});
	}
	

	/// @brief compute MVC on point inside the cage 
	/// state-of-the-art method 
	/// [Mean Value Coordinates for Closed Triangular Meshes, Ju et al. 2005]
	/// @param surface_point 
	/// @param result_weights 
	/// @return 
	bool compute_mvc_on_point_inside_cage(const Vec3& surface_point, 
										Eigen::VectorXd& result_weights)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		result_weights.resize(nbv_cage); 

		double sumWeights = 0.0;
		double epsilon = 0.00000001;

		Eigen::VectorXd w_control_cage_coords_;

		w_control_cage_coords_.resize(nbv_cage);
		w_control_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*control_cage_, 
										control_cage_vertex_position_, v);

			uint32 cage_point_index = value<uint32>(*control_cage_, 
											control_cage_vertex_index_, v);

			d[cage_point_index] = (surface_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				result_weights[cage_point_index] = 1.0; 
				return true;
			}

			u[cage_point_index] = 
				(cage_point - surface_point) / d[cage_point_index];

			return true;
		});

		double l[3], theta[3], w[3], c[3], s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			std::vector<uint32> triangle_index = {
				cage_triangles_[t].points[0].control_cage_index,
				cage_triangles_[t].points[1].control_cage_index,
				cage_triangles_[t].points[2].control_cage_index
			};

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = (u[triangle_index[(i + 1) % 3]] 
						- u[triangle_index[(i + 2) % 3]])
						.norm();
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				theta[i] = 2.0 * asin(l[i] / 2.0);
			}

			double h = (theta[0] + theta[1] + theta[2]) / 2.0;
			if (M_PI - h < epsilon)
			{
				for (std::size_t i = 0; i < 3; i++)
				{
					w[i] = sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_control_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_control_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_control_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / 
						(sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) 
						- 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * 
							sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || 
					fabs(s[2]) < epsilon)
			{
				continue; 
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = (theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] 
							- c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					   (2.0 * d[triangle_index[i]] 
					   	* sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_control_cage_coords_[triangle_index[0]] += w[0];
			w_control_cage_coords_[triangle_index[1]] += w[1];
			w_control_cage_coords_[triangle_index[2]] += w[2];
		}

		foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_index = value<uint32>(*control_cage_, 
											control_cage_vertex_index_, v);

			result_weights[cage_point_index] =
				w_control_cage_coords_[cage_point_index] / sumWeights;

			return true;
		});

		return false;
	}


	/// @brief compute MVC on point inside the cage 
	/// state-of-the-art method 
	/// [Mean Value Coordinates for Closed Triangular Meshes, Ju et al. 2005]
	/// @param surface_point 
	/// @param result_weights 
	/// @return 
	bool compute_mvc_on_point_inside_cage_parallel(const Vec3& surface_point, 
										Eigen::VectorXd& result_weights)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		result_weights.resize(nbv_cage); 

		double sumWeights = 0.0;
		double epsilon = 0.00000001;

		Eigen::VectorXd w_control_cage_coords_;

		w_control_cage_coords_.resize(nbv_cage);
		w_control_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*control_cage_, 
										control_cage_vertex_position_, v);

			uint32 cage_point_index = value<uint32>(*control_cage_, 
											control_cage_vertex_index_, v);

			d[cage_point_index] = (surface_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				result_weights[cage_point_index] = 1.0; 
				return true;
			}

			u[cage_point_index] = 
				(cage_point - surface_point) / d[cage_point_index];

			return true;
		});

		double l[3], theta[3], w[3], c[3], s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			std::vector<uint32> triangle_index = {
				cage_triangles_[t].points[0].control_cage_index,
				cage_triangles_[t].points[1].control_cage_index,
				cage_triangles_[t].points[2].control_cage_index
			};

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = (u[triangle_index[(i + 1) % 3]] 
						- u[triangle_index[(i + 2) % 3]])
						.norm();
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				theta[i] = 2.0 * asin(l[i] / 2.0);
			}

			double h = (theta[0] + theta[1] + theta[2]) / 2.0;
			if (M_PI - h < epsilon)
			{
				for (std::size_t i = 0; i < 3; i++)
				{
					w[i] = sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_control_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_control_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_control_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / 
						(sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) 
						- 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * 
							sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || 
					fabs(s[2]) < epsilon)
			{
				continue; 
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = (theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] 
							- c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					   (2.0 * d[triangle_index[i]] 
					   	* sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_control_cage_coords_[triangle_index[0]] += w[0];
			w_control_cage_coords_[triangle_index[1]] += w[1];
			w_control_cage_coords_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_index = value<uint32>(*control_cage_, 
											control_cage_vertex_index_, v);

			result_weights[cage_point_index] =
				w_control_cage_coords_[cage_point_index] / sumWeights;

			return true;
		});

		return false;
	}

	/**
	 * compute MVC on point outside cage 
	 * 	that is point inside a virtual cube 
	 * compute the classical mvc for the points of the virtual cube 
	 * 	that belongs to the local cage 
	 * 	set the other fixed part inside object_fixed_data
	*/

	/// @brief compute MVC on point outside cage (inside a virtual cube)
	/// separate the weights of the control cage points and 
	/// and the virtual ones 
	/// @param surface_point 
	/// @param surface_point_index 
	/// @param virtual_cube_target 
	/// @param object 
	/// @return 
	bool compute_mvc_on_point_outside_cage(const Vec3& surface_point, 
								Eigen::VectorXd& result_weights,
								const Virtual_cube& virtual_cube_target)
	{

		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		result_weights.resize(nbv_cage); 

		double sumWeights = 0.0;
		double epsilon = 0.00000001;

		Eigen::VectorXd w_control_cage_coords_;

		w_control_cage_coords_.resize(nbv_cage);
		w_control_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		std::size_t number_of_points = virtual_cube_target.points.size(); 

		for (std::size_t p = 0; p < number_of_points; p++)
		{
			const Point& cage_point = virtual_cube_target.points[p];

			uint32 cage_point_index = p;

			d[cage_point_index] = 
				(surface_point - cage_point.position).norm();
			if (d[cage_point_index] < epsilon)
			{
				result_weights[cage_point_index] = 1.0;
				return true;
			}

			u[cage_point_index] = 
				(cage_point.position - surface_point) / d[cage_point_index];
		} 

		double l[3], theta[3], w[3], c[3], s[3];

		for (std::size_t t = 0; t < virtual_cube_target.triangles.size(); t++)
		{
			std::vector<uint32> triangle_index = 
					virtual_cube_target.triangles[t].virtual_cage_indices;

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = 
					(u[triangle_index[(i + 1) % 3]] 
					- u[triangle_index[(i + 2) % 3]])
					.norm();
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				theta[i] = 2.0 * asin(l[i] / 2.0);
			}

			double h = (theta[0] + theta[1] + theta[2]) / 2.0;
			if (M_PI - h < epsilon)
			{
				for (std::size_t i = 0; i < 3; i++)
				{
					w[i] = sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_control_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_control_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_control_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = 
					(2.0 * sin(h) * sin(h - theta[i])) / 
					(sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) 
					- 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = 
				u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * 
						sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || 
										fabs(s[2]) < epsilon)
			{
				continue; 
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = 
					(theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] 
					- c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					(2.0 * d[triangle_index[i]] * 
					sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_control_cage_coords_[triangle_index[0]] += w[0];
			w_control_cage_coords_[triangle_index[1]] += w[1];
			w_control_cage_coords_[triangle_index[2]] += w[2];
		} 

		for (std::size_t p = 0; p < number_of_points; p++)
		{
			uint32 cage_point_index = p; 

			result_weights[cage_point_index] =
				w_control_cage_coords_[cage_point_index] / sumWeights;

		}

		return false;
	}


	/// @brief bind influence area of object with Green deformation type
	/// loop on each point of the influence area 
	/// check if point inside control cage => classical Green binding 
	/// otherwise: find corresponding virtual cube, compute Green inside it
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model 
	/// @param object_vertex_index index of the vertices of the model 
	void bind_object_green(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position, 
						CMap2::Attribute<uint32>* object_vertex_index)
	{
		
		/*parallel_foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = 
							value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_index = 
							value<uint32>(object, object_vertex_index, v);

			const double d_x = surface_point
							.dot(local_x_direction_control_planes_.direction),
						 
						d_y = surface_point
							.dot(local_y_direction_control_planes_.direction),
						
						d_z = surface_point
							.dot(local_z_direction_control_planes_.direction);

			const bool valid_x_dir = 
						(d_x <= local_x_direction_control_planes_.d_max &&
						d_x >= local_x_direction_control_planes_.d_min),

					   valid_y_dir = 
					   (d_y <= local_y_direction_control_planes_.d_max &&
						d_y >= local_y_direction_control_planes_.d_min),

					   valid_z_dir = 
					   (d_z <= local_z_direction_control_planes_.d_max &&
						d_z >= local_z_direction_control_planes_.d_min);

			if (valid_x_dir && valid_y_dir && valid_z_dir)
			{
				object_activation_cage_[surface_point_index] = 27; 
				compute_green_on_point_inside_cage(surface_point, 
												surface_point_index);
			}
			else
			{
				Virtual_cube virtual_cube_target = 
					find_virtual_cube_target(d_x, d_y, d_z, 
									valid_x_dir, valid_y_dir, valid_z_dir,
									surface_point_index);

				compute_green_on_point_outside_cage(surface_point, 
												surface_point_index, 
												virtual_cube_target);
			}
			return true;
		}); */
	}

	/// @brief deform the influence area of the object with Green deformation type
	/// rely object_fixed_data to handle case of point inside virtual cube
	/// @param object model to deform 
	/// @param object_vertex_position position of the vertices of the model 
	/// @param object_vertex_index index of the vertices of the model 
	void deform_object_green(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position,
						CMap2::Attribute<uint32>* object_vertex_index)
	{
		/*const std::size_t influence_area_length = 
											object_influence_area_.size();

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			Vertex v = object_influence_area_[i];*/
		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			uint32 object_point_index = 
							value<uint32>(object, object_vertex_index, v);

			
			if (object_activation_cage_[object_point_index] == 27)
			{
				Vec3 new_position = {0.0, 0.0, 0.0};

				foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = 
					value<Vec3>(*control_cage_, 
								control_cage_vertex_position_, cv);
				uint32 cage_point_index = 
					value<uint32>(*control_cage_, 
								control_cage_vertex_index_, cv);

				new_position += 
					object_weights_.position_(object_point_index, 
											cage_point_index) * cage_point;

				return true;
			});

			Vec3 new_normal = {0.0, 0.0, 0.0};

			const auto sqrt8 = sqrt(8);

			for (std::size_t t = 0; t < cage_triangles_.size(); t++)
			{

				std::vector<Vec3> triangle_position(3);
				for (std::size_t i = 0; i < 3; i++)
				{
					triangle_position[i] = cage_triangles_[t].points[i].position; 
				}

				const Vec3 t_normal = cage_triangles_[t].normal;
				const auto t_u0 = cage_triangles_[t].edges.first;
				const auto t_v0 = cage_triangles_[t].edges.second;

				const auto t_u1 = 
					triangle_position[1] - triangle_position[0];
				const auto t_v1 = 
					triangle_position[2] - triangle_position[1];

				const auto area_face = (t_u0.cross(t_v0)).norm() * 0.5;
				double t_sj =
					sqrt((t_u1.squaredNorm()) * (t_v0.squaredNorm()) 
							- 2.0 * (t_u1.dot(t_v1)) * (t_u0.dot(t_v0)) +
						 (t_v1.squaredNorm()) * (t_u0.squaredNorm())) /
							(sqrt8 * area_face);

				new_normal += 
					object_weights_.normal_(object_point_index, t) * t_sj * t_normal;
			}

			value<Vec3>(object, object_vertex_position, v) = 
												new_position + new_normal;
			
			} else {
				
				Virtual_cube virtual_cube_target = 
						virtual_cubes_[object_activation_cage_[object_point_index]];  

				Vec3 new_position = {0.0, 0.0, 0.0};
				for (size_t p = 0; p < virtual_cube_target.points.size(); p++){
					const Point& cage_point = virtual_cube_target.points[p];
					uint32 cage_point_index = p;

					new_position += 
						object_weights_.position_(object_point_index, p) * 
									cage_point.position;

				}

				Vec3 new_normal = {0.0, 0.0, 0.0};

				const auto sqrt8 = sqrt(8);

				for (std::size_t t = 0; t < virtual_cube_target.triangles.size(); t++)
				{
					Triangle local_triangle = virtual_cube_target.triangles[t]; 

					std::vector<Vec3> triangle_position(3);
					for (std::size_t i = 0; i < 3; i++)
					{
						triangle_position[i] = local_triangle.points[i].position;
					}

					const Vec3 t_normal = local_triangle.normal;
					const auto t_u0 = local_triangle.edges.first;
					const auto t_v0 = local_triangle.edges.second;

					const auto t_u1 = 
						triangle_position[1] - triangle_position[0];
					const auto t_v1 = 
						triangle_position[2] - triangle_position[1];

					const auto area_face = (t_u0.cross(t_v0)).norm() * 0.5;
					double t_sj =
						sqrt((t_u1.squaredNorm()) * (t_v0.squaredNorm()) 
							- 2.0 * (t_u1.dot(t_v1)) * (t_u0.dot(t_v0)) +
						 (t_v1.squaredNorm()) * (t_u0.squaredNorm())) /
							(sqrt8 * area_face);

					new_normal += 
							object_weights_.normal_(object_point_index, t) * t_sj * t_normal;
			}

			value<Vec3>(object, object_vertex_position, v) = 
													new_position + new_normal;
			 
			}

			return true; 
		}); 
	}


	/// @brief compute Green on point inside the cage 
	/// state-of-the-art method
	/// [Green coordinates, Lipman et al. 2008]
	/// @param surface_point 
	/// @param surface_point_index 
	void compute_green_on_point_inside_cage(const Vec3& surface_point, 
										const uint32& surface_point_index)
	{

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			std::vector<Vec3> triangle_position(3);

			for (std::size_t i = 0; i < 3; i++)
			{
				triangle_position[i] = cage_triangles_[t].points[i].position; 
			}

			const Vec3 t_normal = cage_triangles_[t].normal; 

			std::vector<Vec3> t_vj(3);
			for (std::size_t l = 0; l < 3; ++l)
			{
				t_vj[l] = triangle_position[l] - surface_point;
			}

			const Vec3 t_p_ = (t_vj[0].dot(t_normal)) * t_normal;

			Vec3 t_I = {0.0, 0.0, 0.0};
			std::vector<double> t_II(3);
			Vec3 t_s = {0.0, 0.0, 0.0};
			std::vector<Vec3> t_N(3);

			for (std::size_t l = 0; l < 3; ++l)
			{
				const Vec3 t_v0 = t_vj[l];
				const Vec3 t_v1 = t_vj[(l + 1) % 3];

				const auto t_vjpt = 
					((t_v0 - t_p_).cross((t_v1 - t_p_))).dot(t_normal);
				t_s[l] = t_vjpt < 0 ? -1.0 : 1.0;

				t_I[l] = modeling::GCTriInt2(t_p_, t_v0, t_v1);
				t_II[l] = modeling::GCTriInt2(NULL_VECTOR, t_v1, t_v0);
				t_N[l] = (t_v1.cross(t_v0)).normalized();
			}

			const double t_I_ = -abs(t_s.dot(t_I));

			object_weights_.normal_(surface_point_index, t) = -t_I_;

			Vec3 t_w = t_I_ * t_normal;

			for (std::size_t k = 0; k < 3; ++k)
			{
				t_w += (t_II[k] * t_N[k]);
			}

			if (t_w.norm() > DBL_EPSILON)
			{
				for (size_t l = 0; l < 3; l++)
				{
					const uint32 cage_point_index = 
								cage_triangles_[t].points[l].control_cage_index; 

					const Vec3 Nl = t_N[(l + 1) % 3];
					const double num = Nl.dot(t_w);
					const double denom = Nl.dot(t_vj[l]);

					object_weights_
						.position_(surface_point_index, cage_point_index) 
										+= num / denom;
				}
			}
		}
	}

	


	/// @brief compute Green on point outside the cage 
	/// TODO: fix the position weights 
	/// state-of-the-art method
	/// [Green coordinates, Lipman et al. 2008]
	/// @param surface_point 
	/// @param surface_point_index 
	void compute_green_on_point_outside_cage(const Vec3& surface_point, 
										const uint32& surface_point_index,
										const Virtual_cube& virtual_cube_target)
	{

		
		for (std::size_t t = 0; t < virtual_cube_target.triangles.size(); t++)
		{
			Triangle local_triangle = virtual_cube_target.triangles[t]; 

			std::vector<Vec3> triangle_position(3);

			for (std::size_t i = 0; i < 3; i++)
			{
				triangle_position[i] = local_triangle.points[i].position; 
			}

			const Vec3 t_normal = -local_triangle.normal; 

			std::vector<Vec3> t_vj(3);
			for (std::size_t l = 0; l < 3; ++l)
			{
				t_vj[l] = triangle_position[l] - surface_point;
			}

			const Vec3 t_p_ = (t_vj[0].dot(t_normal)) * t_normal;

			Vec3 t_I = {0.0, 0.0, 0.0};
			std::vector<double> t_II(3);
			Vec3 t_s = {0.0, 0.0, 0.0};
			std::vector<Vec3> t_N(3);

			for (std::size_t l = 0; l < 3; ++l)
			{
				const Vec3 t_v0 = t_vj[l];
				const Vec3 t_v1 = t_vj[(l + 1) % 3];

				const auto t_vjpt = 
					((t_v0 - t_p_).cross((t_v1 - t_p_))).dot(t_normal);
				t_s[l] = t_vjpt < 0 ? -1.0 : 1.0;

				t_I[l] = modeling::GCTriInt2(t_p_, t_v0, t_v1);
				t_II[l] = modeling::GCTriInt2(NULL_VECTOR, t_v1, t_v0);
				t_N[l] = (t_v1.cross(t_v0)).normalized();
			}

			const double t_I_ = -abs(t_s.dot(t_I));

			object_weights_.normal_(surface_point_index, t) = -t_I_;

			Vec3 t_w = t_I_ * t_normal;

			for (std::size_t k = 0; k < 3; ++k)
			{
				t_w += (t_II[k] * t_N[k]);
			}

			if (t_w.norm() > DBL_EPSILON)
			{
				for (size_t l = 0; l < 3; l++)
				{
					const uint32 cage_point_index = 
							local_triangle.virtual_cage_indices[l]; 

					const Vec3 Nl = t_N[(l + 1) % 3];
					const double num = Nl.dot(t_w);
					const double denom = Nl.dot(t_vj[l]);

					object_weights_
						.position_(surface_point_index, cage_point_index) 
										+= num / denom;

				}
				 
			}

		}
	}


	/// @brief check if a position is inside the influence cage 
	/// use the local_direction_control_plane structure
	/// to verify if a position of the model is between the planes 
	/// of the influence cage 
	/// @param surface_point 
	/// @return true if inside the influence cage, false otherwise 
	bool check_point_inside_influence_cage(const Vec3& surface_point)
	{
		const double d_x = surface_point
						.dot(local_x_direction_control_planes_.direction),
					 d_y = surface_point
					 	.dot(local_y_direction_control_planes_.direction),
					 d_z = surface_point
					 	.dot(local_z_direction_control_planes_.direction);

		const double influence_cage_d_x_min = 
			local_x_direction_control_planes_.d_min 
				- local_x_direction_control_planes_.d_gap;
		const double influence_cage_d_x_max =
			local_x_direction_control_planes_.d_max 
				+ local_x_direction_control_planes_.d_gap;

		const double influence_cage_d_y_min =
			local_y_direction_control_planes_.d_min 
				- local_y_direction_control_planes_.d_gap;
		const double influence_cage_d_y_max =
			local_y_direction_control_planes_.d_max 
				+ local_y_direction_control_planes_.d_gap;

		const double influence_cage_d_z_min =
			local_z_direction_control_planes_.d_min 
				- local_z_direction_control_planes_.d_gap;
		const double influence_cage_d_z_max =
			local_z_direction_control_planes_.d_max 
				+ local_z_direction_control_planes_.d_gap;

		const bool valid_x_dir = 
					(d_x <= influence_cage_d_x_max && 
						d_x >= influence_cage_d_x_min),

				   valid_y_dir = 
				   (d_y <= influence_cage_d_y_max && 
				   		d_y >= influence_cage_d_y_min),

				   valid_z_dir = 
				   (d_z <= influence_cage_d_z_max && 
				   		d_z >= influence_cage_d_z_min);

		return (valid_x_dir && valid_y_dir && valid_z_dir);
	}


	/// @brief find the intersections points between two faces
	/// @param face1 
	/// @param face2 
	/// @return common edge of the provided faces
	std::vector<Point> find_intersection_points_face(
			const std::pair<Triangle, Triangle>& face1,
			const std::pair<Triangle, Triangle>& face2)
	{

		std::vector<Point> intersect_points;
		double epsilon = 0.00000001;

		std::vector<Point> points_face_1(4);
		points_face_1[0] = face1.first.points[2];
		points_face_1[1] = face1.first.points[0];
		points_face_1[2] = face1.second.points[1];
		points_face_1[3] = face1.first.points[1];

		std::vector<Point> points_face_2(4);
		points_face_2[0] = face2.first.points[2];
		points_face_2[1] = face2.first.points[0];
		points_face_2[2] = face2.second.points[1];
		points_face_2[3] = face2.first.points[1];

		for (std::size_t i = 0; i < 4; i++)
		{
			Point target_point_face1 = points_face_1[i];
			for (std::size_t j = 0; j < 4; j++)
			{
				Point target_point_face2 = points_face_2[j];

				const double distance_p1_p2 = 
					(target_point_face1.position - target_point_face2.position)
					.norm();
				if (distance_p1_p2 < epsilon)
				{
					intersect_points.push_back(target_point_face1);
					continue;
				}
			}
		}

		return intersect_points;
	}

	/// @brief find intersection point between three faces 
	/// @param face1 
	/// @param face2 
	/// @param face3 
	/// @return common point intersecting the three provided faces
	Point find_intersection_point(const std::pair<Triangle, Triangle> face1, 
								const std::pair<Triangle, Triangle> face2,
								const std::pair<Triangle, Triangle> face3)
	{
		Point intersection_point;

		double epsilon = 0.00000001;

		std::vector<Point> points_face_1(4);
		points_face_1[0] = face1.first.points[2];
		points_face_1[1] = face1.first.points[0];
		points_face_1[2] = face1.second.points[1];
		points_face_1[3] = face1.first.points[1];

		std::vector<Point> points_face_2(4);
		points_face_2[0] = face2.first.points[2];
		points_face_2[1] = face2.first.points[0];
		points_face_2[2] = face2.second.points[1];
		points_face_2[3] = face2.first.points[1];

		std::vector<Point> points_face_3(4);
		points_face_3[0] = face3.first.points[2];
		points_face_3[1] = face3.first.points[0];
		points_face_3[2] = face3.second.points[1];
		points_face_3[3] = face3.first.points[1];

		for (std::size_t i = 0; i < 4; i++)
		{
			Point target_point_face1 = points_face_1[i];
			for (std::size_t j = 0; j < 4; j++)
			{
				Point target_point_face2 = points_face_2[j];
				const double distance_p1_p2 = 
					(target_point_face1.position 
						- target_point_face2.position)
					.norm();
				if (distance_p1_p2 < epsilon)
				{
					for (std::size_t k = 0; k < 4; k++)
					{
						Point target_point_face3 = points_face_3[k];
						const double distance_p1_p3 =
							(target_point_face1.position 
								- target_point_face3.position)
							.norm();
						if (distance_p1_p3 < epsilon)
						{
							return target_point_face1;
						}
					}
				}
			}
		}

		return intersection_point;
	}


	/// @brief get face from intersecting edge 
	/// create new face composed of the intersecting edge (in the vector)
	/// and this edge shifted by the vector
	/// @param intersect_points 
	/// @param shiftVector 
	/// @return Face 
	std::pair<Triangle, Triangle> get_face_from_intersecting_edge(
								const std::vector<Point>& intersect_points,
								const Vec3& shiftVector, const Vec3& unit_dir)
	{
		
		const Point face_point0 = intersect_points[0];
		const Point face_point1 = intersect_points[1];

		Point face_point2;
		face_point2.position = face_point1.position + shiftVector;
		face_point2.inside_control_cage = false;
		face_point2.shift_vector = {0.0, 0.0, 0.0}; 
		face_point2.shift_vector += face_point1.shift_vector; 
		face_point2.shift_vector += unit_dir; 
		face_point2.vertex = face_point1.vertex; 

		Point face_point3;
		face_point3.position = face_point0.position + shiftVector;
		face_point3.inside_control_cage = false;
		face_point3.shift_vector = {0.0, 0.0, 0.0};
		face_point3.shift_vector += face_point0.shift_vector; 
		face_point3.shift_vector += unit_dir;
		face_point3.vertex = face_point0.vertex; 

		Triangle local_triangle1, local_triangle2;
		local_triangle1.points = {face_point1, face_point3, face_point0};
		local_triangle1.normal = get_normal_from_triangle(local_triangle1); 
		local_triangle1.edges = get_edges_from_triangle(local_triangle1);

		local_triangle2.points = {face_point1, face_point2, face_point3};
		local_triangle2.normal = get_normal_from_triangle(local_triangle2);
		local_triangle2.edges = get_edges_from_triangle(local_triangle2);

		return std::make_pair(local_triangle1, local_triangle2);
	}


	/// @brief get face from intersecting vertex
	/// create new face composed of the intersection point
	/// and the shift of this point with the two provided directions
	/// @param intersection_point 
	/// @param shiftVector1 
	/// @param shiftVector2 
	/// @return Face 
	std::pair<Triangle, Triangle> get_face_from_intersecting_vertex(
											const Point& intersection_point,
											const Vec3& shiftVector1, 
											const Vec3& shiftVector2,
											const Vec3& unit_dir1, 
											const Vec3& unit_dir2)
	{
		const Point face_point0 = intersection_point;

		Point face_point1, face_point2, face_point3;
		face_point1.position = face_point0.position + shiftVector1;
		face_point1.inside_control_cage = false;
		face_point1.shift_vector = {0.0, 0.0, 0.0};
		face_point1.shift_vector += face_point0.shift_vector; 
		face_point1.shift_vector += unit_dir1; 
		face_point1.vertex = face_point0.vertex; 

		face_point2.position = 
						(face_point0.position + shiftVector1) + shiftVector2;
		face_point2.inside_control_cage = false;
		face_point2.shift_vector = {0.0, 0.0, 0.0};
		face_point2.shift_vector += face_point0.shift_vector; 
		face_point2.shift_vector += (unit_dir1 + unit_dir2); 
		face_point2.vertex = face_point0.vertex; 

		face_point3.position = face_point0.position + shiftVector2;
		face_point3.inside_control_cage = false;
		face_point3.shift_vector = {0.0, 0.0, 0.0};
		face_point3.shift_vector += face_point0.shift_vector; 
		face_point3.shift_vector += unit_dir2;
		face_point3.vertex = face_point0.vertex; 

		Triangle local_triangle1, local_triangle2;
		local_triangle1.points = {face_point1, face_point3, face_point0};
		local_triangle2.points = {face_point1, face_point2, face_point3};

		return std::make_pair(local_triangle1, local_triangle2);
	}


	Vec3 get_normal_from_triangle(const Triangle& local_triangle)
	{
		return (cgogn::geometry::normal(local_triangle.points[0].position, 
										local_triangle.points[1].position,
										local_triangle.points[2].position))
										.normalized();
	}

	std::pair <Vec3, Vec3> get_edges_from_triangle(const Triangle& local_triangle)
	{
		return std::make_pair(local_triangle.points[1].position 
								- local_triangle.points[0].position, 
							local_triangle.points[2].position 
								- local_triangle.points[1].position);
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
