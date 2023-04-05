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

#ifndef CGOGN_MODELING_GLOBAL_CAGE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_GLOBAL_CAGE_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cells_set.h>
#include <cgogn/core/types/cmap/cmap2.h>

#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>

#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <boost/synapse/connect.hpp>

namespace cgogn
{

namespace modeling
{

/**
 * @Class Global Cage Deformation Tool
 */
template <typename MESH>
class GlobalCageDeformationTool
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

	using Graph = IncidenceGraph;

	template <typename T>
	using GraphAttribute = typename mesh_traits<IncidenceGraph>::
													template Attribute<T>;
	using GraphVertex = IncidenceGraph::Vertex;
	using GraphEdge = IncidenceGraph::Edge;
	using GraphFace = IncidenceGraph::Face;


public:
	MESH* global_cage_;
	std::shared_ptr<Attribute<Vec3>> global_cage_vertex_position_;

	MatrixWeights object_weights_;
	std::unordered_map<std::string, MatrixWeights> local_cage_weights_; 
	std::unordered_map<std::string, MatrixWeights> local_axis_weights_; 
	std::unordered_map<std::string, VectorWeights> local_handle_weights_; 

	std::string deformation_type_;

	std::shared_ptr<boost::synapse::connection> 
										cage_attribute_update_connection_;

	GlobalCageDeformationTool() : global_cage_vertex_position_(nullptr)
	{
	}

	~GlobalCageDeformationTool()
	{
	}

	/// @brief create global cage space tool
	/// initialize the set of triangles from the square faces
	/// set the global cage start position used to reset the deformation
	/// @param m default mesh
	/// @param m_vertex_position position of the vertices of the default graph 
	/// @param bb_min minimum coordinates of the model bounding box
	/// @param bb_max maximum coordinates of the model bounding box
	void create_global_cage(MESH* m,
							CMap2::Attribute<Vec3>* m_vertex_position, 
							const Vec3& bb_min, const Vec3& bb_max)
	{
		global_cage_ = m;

		global_cage_vertices_ = modeling::create_bounding_box(*m, 
											m_vertex_position, bb_min, bb_max);

		global_cage_vertex_position_ = get_attribute<Vec3, Vertex>(*m, 
																"position");

		global_cage_vertex_index_ = add_attribute<uint32, Vertex>(*global_cage_, 
														"vertex_index");
		modeling::set_attribute_vertex_index(*global_cage_, 
											global_cage_vertex_index_.get());

		init_triangles(); 

		set_start_positions(); 
		
	}

	/// @brief update the global cage dimensions 
	/// call when the content spatial tool deformation goes outside the cage
	/// @param bb_min 
	/// @param bb_max 
	void update_global_cage(const Vec3& bb_min, const Vec3& bb_max)
	{
		modeling::update_bounding_box(*global_cage_, 
									global_cage_vertex_position_.get(), 
									global_cage_vertices_, 
									bb_min, bb_max);
	}

	/// @brief set deformation type
	/// choice between MVC and Green
	/// @param new_type MVC or Green
	void set_deformation_type(const std::string& new_type)
	{
		deformation_type_ = new_type; 
	}

	/// @brief reset global cage deformation
	/// allow to change binding deformation type 
	void reset_deformation()
	{
		foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);

			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

			value<Vec3>(*global_cage_, global_cage_vertex_position_, cv) = 		start_positions_[cage_point_index]; 

			return true;
		});
	}

	/// @brief initialize binding of the global cage on the model to deform
	/// call the corresponding binding function relative to chosen type 
	/// @param object model to deform 
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	void init_bind_object(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position,
						CMap2::Attribute<uint32>* object_vertex_index)
	{
		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		object_weights_.position_.resize(nbv_object, nbv_cage);
		object_weights_.position_.setZero();

		if (deformation_type_ == "MVC")
		{
			bind_object_mvc(object, object_vertex_position, 
									object_vertex_index); 
		}

		if (deformation_type_ == "Green")
		{
			uint32 nbt_cage = cage_triangles_.size();
			object_weights_.normal_.resize(nbv_object, nbt_cage);
			object_weights_.normal_.setZero();

			bind_object_green(object, object_vertex_position, 
										object_vertex_index);
		}
	}

	/// @brief update binding weights 
	/// call when the deformation type is changed (after reset)
	/// @param object model to deform 
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	void bind_object(MESH& object, 
					CMap2::Attribute<Vec3>* object_vertex_position,
					 CMap2::Attribute<uint32>* object_vertex_index)
	{
		object_weights_.position_.setZero();
		if (deformation_type_ == "MVC")
		{
			bind_object_mvc(object, object_vertex_position, 
										object_vertex_index); 
		}

		if (deformation_type_ == "Green")
		{
			if (object_weights_.normal_ == Eigen::MatrixXd{}){
				int32 nbv_object = nb_cells<Vertex>(object);
				int32 nbt_cage = cage_triangles_.size();
				object_weights_.normal_.resize(nbv_object, nbt_cage);
			}
			object_weights_.normal_.setZero();
			bind_object_green(object, object_vertex_position,
										 object_vertex_index);
		}
	}

	/// @brief update weights of local area of object 
	/// call when a content spatial tool is deformed 
	/// and there is a need to update its zone of influence weights with the ///global cage 
	/// call the specific deformation type function 
	/// @param object model to deform 
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	/// @param influence_area zone of influence of the content tool 
	void update_local_area_object(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position, 
						CMap2::Attribute<uint32>* object_vertex_index, 
						ui::CellsSet<MESH, Vertex>* influence_area) 
	{ 
		if (deformation_type_ == "MVC")
		{
			update_local_area_object_mvc(object, object_vertex_position, 
										object_vertex_index, influence_area); 
		}

		if (deformation_type_ == "Green")
		{
			update_local_area_object_green(object, object_vertex_position,
										object_vertex_index, influence_area);
		}
	}

	/// @brief deform the model 
	/// call the specific deformation type function 
	/// @param object model to deform 
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	void deform_object(MESH& object, 
					CMap2::Attribute<Vec3>* object_vertex_position,
					CMap2::Attribute<uint32>* object_vertex_index)
	{
		 
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

	/// @brief initialize binding for handle tool 
	/// @param graph_name name of the handle 
	/// @param handle_position position of the handle 
	void init_bind_handle(const std::string& graph_name, 
						const Vec3& handle_position)
	{

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		VectorWeights handle_weights; 
		 
		handle_weights.position_.resize(nbv_cage);
		handle_weights.position_.setZero();
		local_handle_weights_[graph_name] = handle_weights;

		if (deformation_type_ == "MVC"){
			bind_handle_mvc(graph_name, handle_position); 
		}

		if (deformation_type_ == "Green"){
			uint32 nbt_cage = cage_triangles_.size();

			local_handle_weights_[graph_name].normal_.resize(nbt_cage); 
			local_handle_weights_[graph_name].normal_.setZero();

			bind_handle_green(graph_name, handle_position); 
		}
	}

	/// @brief update binding for handle tool 
	/// used when deformation type is changed
	/// @param graph_name handle name 
	/// @param handle_position handle position 
	void bind_handle(const std::string& graph_name, 
					const Vec3& handle_position)
	{
		local_handle_weights_[graph_name].position_.setZero();
		if (deformation_type_ == "MVC"){
			bind_handle_mvc(graph_name, handle_position); 
		}

		if (deformation_type_ == "Green"){
			if (local_handle_weights_[graph_name].normal_ == Eigen::MatrixXd{}){
				uint32 nbt_cage = cage_triangles_.size();
				local_handle_weights_[graph_name].normal_.resize(nbt_cage); 
			} 

			local_handle_weights_[graph_name].normal_.setZero();
			bind_handle_green(graph_name, handle_position); 
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

		if (deformation_type_ == "Green")
		{
			deform_handle_green(g, graph_name, graph_vertex_position, 
								handle_vertex);
		}
	}


	/// @brief initialize binding for axis tool 
	/// @param g axis
	/// @param graph_name axis name  
	/// @param graph_vertex_position position of the axis vertices  
	/// @param axis_vertices axis vertices 
	void init_bind_axis(Graph& g, const std::string& graph_name, 
			std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
						const std::vector<Graph::Vertex>& axis_vertices)
	{
		uint32 nbv_axis = axis_vertices.size();
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		MatrixWeights axis_weights; 
		axis_weights.position_.resize(nbv_axis, nbv_cage);
		axis_weights.position_.setZero();
		local_axis_weights_[graph_name] = axis_weights;

		if (deformation_type_ == "MVC")
		{
			bind_axis_mvc(g, graph_name, graph_vertex_position, 
										axis_vertices); 
		}

		if (deformation_type_ == "Green")
		{
			uint32 nbt_cage = cage_triangles_.size();
			local_axis_weights_[graph_name].normal_
										.resize(nbv_axis, nbt_cage); 
			local_axis_weights_[graph_name].normal_.setZero();

			bind_axis_green(g, graph_name, graph_vertex_position,
										 axis_vertices);
		}
	}

	/// @brief update binding for axis tool 
	/// @param g axis 
	/// @param graph_name axis name  
	/// @param graph_vertex_position position of the axis vertices 
	/// @param axis_vertices axis vertices 
	void bind_axis(Graph& g, const std::string& graph_name, 
			std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
					const std::vector<Graph::Vertex>& axis_vertices)
	{
		local_axis_weights_[graph_name].position_.setZero();
		if (deformation_type_ == "MVC"){
			bind_axis_mvc(g, graph_name, graph_vertex_position, 
										axis_vertices); 
		}

		if (deformation_type_ == "Green"){
			if (local_axis_weights_[graph_name].normal_ == Eigen::MatrixXd{}){
				uint32 nbt_cage = cage_triangles_.size();
				uint32 nbv_axis = axis_vertices.size();
				local_axis_weights_[graph_name].normal_.resize(nbv_axis, nbt_cage); 
			} 

			local_axis_weights_[graph_name].normal_.setZero();
			bind_axis_green(g, graph_name, graph_vertex_position, 
										axis_vertices); 
		}
	}

	/// @brief deform the axis tool 
	/// call the corresponding deformation type function 
	/// @param g axis 
	/// @param graph_name axis name  
	/// @param graph_vertex_position position of the axis vertices 
	/// @param axis_vertices axis vertices
	void deform_axis(Graph& g, const std::string& graph_name, 
			std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position, 
				const std::vector<GraphVertex>& axis_vertices)
	{
		if (deformation_type_ == "MVC")
		{
			deform_axis_mvc(g, graph_name, graph_vertex_position, 
								axis_vertices); 
		}

		if (deformation_type_ == "Green")
		{
			deform_axis_green(g, graph_name, graph_vertex_position, 
								axis_vertices);
		}
	}

	/// @brief initialize binding for local cage 
	/// @param local_cage 
	/// @param cage_name local cage name 
	/// @param local_cage_vertex_position position of the vertices of the cage
	/// @param local_cage_vertex_index indices of the vertices of the cage
	void init_bind_local_cage(MESH& local_cage, 
							const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		uint32 nbv_local_cage = nb_cells<Vertex>(local_cage);
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		MatrixWeights local_cage_weights; 
		local_cage_weights.position_.resize(nbv_local_cage, nbv_cage);
		local_cage_weights.position_.setZero();
		local_cage_weights_[cage_name] = local_cage_weights;

		if (deformation_type_ == "MVC")
		{
			bind_local_cage_mvc(local_cage, cage_name, 
						local_cage_vertex_position, local_cage_vertex_index); 
		}

		if (deformation_type_ == "Green")
		{
			uint32 nbt_cage = cage_triangles_.size();
			local_cage_weights_[cage_name].normal_
										.resize(nbv_local_cage, nbt_cage); 
			local_cage_weights_[cage_name].normal_.setZero();

			bind_local_cage_green(local_cage, cage_name, 
						local_cage_vertex_position, local_cage_vertex_index);
		}
	}

	/// @brief update binding for local cage tool 
	/// used when deformation type is changed (after reset)
	/// @param local_cage 
	/// @param cage_name local cage name 
	/// @param local_cage_vertex_position positions of the vertices of the cage
	/// @param local_cage_vertex_index indices of the vertices of the cage
	void bind_local_cage(MESH& local_cage, const std::string& cage_name, 
	std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
	std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		local_cage_weights_[cage_name].position_.setZero();
		if (deformation_type_ == "MVC")
		{
			bind_local_cage_mvc(local_cage, cage_name, 
						local_cage_vertex_position, local_cage_vertex_index); 
		}

		if (deformation_type_ == "Green")
		{
			if (local_axis_weights_[cage_name].normal_ == Eigen::MatrixXd{}){
				uint32 nbt_cage = cage_triangles_.size();
				uint32 nbv_local_cage = nb_cells<Vertex>(local_cage);

				local_axis_weights_[cage_name].normal_.resize(nbv_local_cage, nbt_cage); 
			} 
			local_cage_weights_[cage_name].normal_.setZero();
			bind_local_cage_green(local_cage, cage_name, 
						local_cage_vertex_position, local_cage_vertex_index);
		}
	}

	/// @brief deform the local cage 
	/// call the corresponding deformation type function 
	/// @param local_cage 
	/// @param cage_name local cage name 
	/// @param local_cage_vertex_position position of the vertices of the cage 
	/// @param local_cage_vertex_index indices of the vertices of the cage
	void deform_local_cage(MESH& local_cage, const std::string& cage_name, 
	std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
	std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		 
		if (deformation_type_ == "MVC")
		{
			deform_local_cage_mvc(local_cage, cage_name, 
						local_cage_vertex_position, local_cage_vertex_index); 
		}

		if (deformation_type_ == "Green")
		{
			deform_local_cage_green(local_cage, cage_name, 
						local_cage_vertex_position, local_cage_vertex_index);
		}
	}



private:
	/**
	 * private structure triangle 
	 * different than local cage triangle
	 * no virtual cube here
	*/
	struct Triangle
	{
		std::vector<CMap2::Vertex> vertices_; 
		std::vector<Vec3> positions_;
		std::vector<uint32> indices_;
		Vec3 normal_;
		std::pair<Vec3, Vec3> edges_;
	}; 

	std::vector<Triangle> cage_triangles_;
	std::vector<CMap2::Vertex> global_cage_vertices_;
	std::shared_ptr<Attribute<uint32>> global_cage_vertex_index_;

	std::vector<Vec3> start_positions_; 
	
	/// @brief set start positions to reset the deformation
	///
	void set_start_positions()
	{
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		start_positions_.resize(nbv_cage);

		foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);

			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

			start_positions_[cage_point_index] = cage_point; 

			return true;
		});
	}

	
	/// @brief initialize the global cage triangles from the square faces
	void init_triangles()
	{
		foreach_cell(*global_cage_, [&](Face fc) -> bool {
			std::vector<CMap2::Vertex> face_vertices_ = 
										incident_vertices(*global_cage_, fc);

			std::vector<CMap2::Vertex> triangle1_vertices(3); 
			triangle1_vertices[0] = face_vertices_[1]; 
			triangle1_vertices[1] = face_vertices_[3]; 
			triangle1_vertices[2] = face_vertices_[0];

			std::vector<CMap2::Vertex> triangle2_vertices(3); 
			triangle2_vertices[0] = face_vertices_[1]; 
			triangle2_vertices[1] = face_vertices_[2]; 
			triangle2_vertices[2] = face_vertices_[3];


			std::vector<Vec3> triangle1_positions, triangle2_positions;
			std::vector<uint32> triangle1_indices, triangle2_indices;

			for (std::size_t i = 0; i < 3; i++){
				triangle1_positions.push_back(value<Vec3>(*global_cage_, 
						global_cage_vertex_position_, triangle1_vertices[i])); 
				
				triangle2_positions.push_back(value<Vec3>(*global_cage_, 
						global_cage_vertex_position_, triangle2_vertices[i])); 

				triangle1_indices.push_back(value<uint32>(*global_cage_, 
						global_cage_vertex_index_, triangle1_vertices[i]));
				
				triangle2_indices.push_back(value<uint32>(*global_cage_, 
						global_cage_vertex_index_, triangle2_vertices[i])); 
			}

			Triangle triangle1; 
			triangle1.vertices_ = triangle1_vertices; 
			triangle1.positions_ = triangle1_positions; 
			triangle1.indices_ = triangle1_indices; 
			triangle1.normal_ = 
				(geometry::normal(triangle1_positions[0], 
								triangle1_positions[1], 
								triangle1_positions[2]))
							.normalized(); 
			triangle1.edges_ = 
				std::make_pair(triangle1_positions[1] 
									- triangle1_positions[0], 
							triangle1_positions[2] 
									- triangle1_positions[1]);  
			cage_triangles_.push_back(triangle1); 

			Triangle triangle2; 
			triangle2.vertices_ = triangle2_vertices; 
			triangle2.positions_ = triangle2_positions;
			triangle2.indices_ = triangle2_indices; 
			triangle2.normal_ = 
				(geometry::normal(triangle2_positions[0], 
								triangle2_positions[1], 
								triangle2_positions[2]))
							.normalized();

			triangle2.edges_ = 
				std::make_pair(triangle2_positions[1] 
								- triangle2_positions[0], 
							triangle2_positions[2] 
								- triangle2_positions[1]);  
			cage_triangles_.push_back(triangle2); 

			return true;
		});
	}


	/// @brief bind object with MVC
	/// @param object model to deform 
	/// @param object_vertex_position positions of the vertices of the model  
	/// @param object_vertex_index indices of the vertices of the model 
	void bind_object_mvc(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position, 
						CMap2::Attribute<uint32>* object_vertex_index)
	{
		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, 
												object_vertex_position, v);
			const uint32& surface_point_index = 
							value<uint32>(object, object_vertex_index, v);

			compute_mvc_coordinates_on_point(surface_point, 
														surface_point_index);

			return true;
		});
	}

	/// @brief update the MVC weights of a local area of the model 
	/// useful when content spatial tools is deformed and then need to update the weights of its zone of influence
	/// @param object model to deform 
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	/// @param influence_area zone of influence
	void update_local_area_object_mvc(MESH& object, 
							CMap2::Attribute<Vec3>* object_vertex_position,
							CMap2::Attribute<uint32>* object_vertex_index,
							ui::CellsSet<MESH, Vertex>* influence_area)
	{
		influence_area->foreach_cell([&](Vertex v) {
			const Vec3& surface_point = value<Vec3>(object, 
												object_vertex_position, v);
			uint32 surface_point_index = value<uint32>(object, 
													object_vertex_index, v);

			compute_mvc_coordinates_on_point(surface_point, 
													surface_point_index);
		});
	}

 
	/// @brief deform model with MVC deformation type
	/// @param object model to deform 
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	void deform_object_mvc(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position,
						CMap2::Attribute<uint32>* object_vertex_index)
	{
		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			uint32 object_point_index = value<uint32>(object, 
													object_vertex_index, v);

			Vec3 new_position = {0.0, 0.0, 0.0};

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);

				uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

				new_position += object_weights_
					.position_(object_point_index, cage_point_index) 
															* cage_point;

				return true;
			});

			value<Vec3>(object, object_vertex_position, v) = new_position;
			return true;
		});
	}

	/// @brief bind model with Green deformation type 
	/// @param object model to deform
	/// @param object_vertex_position positions of the vertices of the model  
	/// @param object_vertex_index indices of the vertices of the model 
	void bind_object_green(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position, 
						CMap2::Attribute<uint32>* object_vertex_index)
	{

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, 
												object_vertex_position, v);
			const uint32& surface_point_index = 
								value<uint32>(object, object_vertex_index, v);

			compute_green_coordinates_on_point(surface_point, 
														surface_point_index);

			return true;
		});
	}


	/// @brief update the Green weights of a local area of the model 
	/// useful when content spatial tools is deformed and then need to update the weights of its zone of influence 
	/// @param object model to deform
	/// @param object_vertex_position positions of the vertices of the model  
	/// @param object_vertex_index indices of the vertices of the model 
	/// @param influence_area zone of influence 
	void update_local_area_object_green(MESH& object, 
							CMap2::Attribute<Vec3>* object_vertex_position,
							CMap2::Attribute<uint32>* object_vertex_index,
							ui::CellsSet<MESH, Vertex>* influence_area)
	{
		influence_area->foreach_cell([&](Vertex v) {
			const Vec3& surface_point = value<Vec3>(object, 
												object_vertex_position, v);
			uint32 surface_point_index = value<uint32>(object, 
													object_vertex_index, v);

			compute_green_coordinates_on_point(surface_point, 
														surface_point_index);
		});
	}

	/// @brief deform model with Green deformation model 
	/// @param object model to deform 
	/// @param object_vertex_position positions of the vertices of the model 
	/// @param object_vertex_index indices of the vertices of the model 
	void deform_object_green(MESH& object, 
							CMap2::Attribute<Vec3>* object_vertex_position,
							CMap2::Attribute<uint32>* object_vertex_index)
	{

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(object, object_vertex_index, v);

			Vec3 new_position = {0.0, 0.0, 0.0};

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);

				uint32 cage_point_idx = value<uint32>(*global_cage_, 
										global_cage_vertex_index_, cv);

				new_position += object_weights_
							.position_(vidx, cage_point_idx) * cage_point;

				return true;
			});

			Vec3 new_normal = {0.0, 0.0, 0.0};

			const auto sqrt8 = sqrt(8);

			for (std::size_t t = 0; t < cage_triangles_.size(); t++)
			{

				std::vector<Vec3> triangle_position(3);
				for (std::size_t i = 0; i < 3; i++)
				{
					triangle_position[i] =
						value<Vec3>(*global_cage_, 
									global_cage_vertex_position_, 
									cage_triangles_[t].vertices_[i]);
				}

				const Vec3 t_normal = cage_triangles_[t].normal_;
				const auto t_u0 = cage_triangles_[t].edges_.first;
				const auto t_v0 = cage_triangles_[t].edges_.second;

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
					object_weights_.normal_(vidx, t) * t_sj * t_normal;
			}

			value<Vec3>(object, object_vertex_position, v) = 
												new_position + new_normal;

			return true;
		});
	}

	/// @brief bind handle with MVC deformation type
	void bind_handle_mvc(const std::string& graph_name, 
						const Vec3& handle_position){

		compute_mvc_coordinates_on_handle(graph_name, handle_position);								
	}

	/// @brief bind handle with Green deformation type
	void bind_handle_green(const std::string& graph_name, 
												const Vec3& handle_position){

		compute_green_coordinates_on_handle(graph_name, handle_position);								
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
		VectorWeights handle_weights = local_handle_weights_[graph_name]; 
		Vec3 new_position = {0.0, 0.0, 0.0};

		foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);
			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

			new_position += 
					handle_weights.position_[cage_point_index] * cage_point;

			return true;
		});

		value<Vec3>(g, graph_vertex_position, handle_vertex) = new_position;
	}

	
	/// @brief deform handle with Green deformation type 
	/// @param g handle 
	/// @param graph_name handle name  
	/// @param graph_vertex_position position of the handle vertex
	/// @param handle_vertex handle vertex
	void deform_handle_green(Graph& g, const std::string& graph_name, 
			std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
							 const GraphVertex& handle_vertex)
	{
		VectorWeights handle_weights = local_handle_weights_[graph_name];
		Vec3 new_position = {0.0, 0.0, 0.0};

		foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);
			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

			new_position += 
					handle_weights.position_[cage_point_index] * cage_point;

			return true;
		});

		Vec3 new_normal = {0.0, 0.0, 0.0};
		const auto sqrt8 = sqrt(8);

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			std::vector<Vec3> triangle_position(3);
			for (std::size_t i = 0; i < 3; i++)
			{
				triangle_position[i] = value<Vec3>(*global_cage_, 
											global_cage_vertex_position_, 
											cage_triangles_[t].vertices_[i]);
			}

			const Vec3 t_normal = cage_triangles_[t].normal_;
			const auto t_u0 = cage_triangles_[t].edges_.first;
			const auto t_v0 = cage_triangles_[t].edges_.second;

			const auto t_u1 = triangle_position[1] - triangle_position[0];
			const auto t_v1 = triangle_position[2] - triangle_position[1];

			const auto area_face = (t_u0.cross(t_v0)).norm() * 0.5;
			double t_sj = 
				sqrt((t_u1.squaredNorm()) * (t_v0.squaredNorm()) 
					- 2.0 * (t_u1.dot(t_v1)) * (t_u0.dot(t_v0)) +
					(t_v1.squaredNorm()) * (t_u0.squaredNorm())) /
						  (sqrt8 * area_face);

			new_normal += handle_weights.normal_[t] * t_sj * t_normal;
		}

		value<Vec3>(g, graph_vertex_position, handle_vertex) = 
												new_position + new_normal;
	}

	/// @brief bind axis with MVC deformation type
	/// @param g axis
	/// @param graph_name axis name 
	/// @param graph_vertex_position positions of the axis vertices 
	/// @param axis_vertices axis vertices
	void bind_axis_mvc(Graph& g, const std::string& graph_name, 
				std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
						const std::vector<Graph::Vertex>& axis_vertices)
	{
		uint32 nbv_axis = axis_vertices.size();

		for (uint32 i = 0; i < nbv_axis; i++)
		{
			Graph::Vertex v = axis_vertices[i];
			const Vec3 axis_point = value<Vec3>(g, 
												graph_vertex_position, v);
			const uint32 axis_point_index = v.index_;

			compute_mvc_coordinates_on_axis_point(graph_name, 
											axis_point, axis_point_index);
		} 
	}

	/// @brief bind axis with Green deformation type 
	/// @param g axis 
	/// @param graph_name axis name  
	/// @param graph_vertex_position positions of the axis vertices 
	/// @param axis_vertices 
	void bind_axis_green(Graph& g, const std::string& graph_name, 
				std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
						const std::vector<Graph::Vertex>& axis_vertices)
	{
		uint32 nbv_axis = axis_vertices.size();

		for (uint32 i = 0; i < nbv_axis; i++)
		{
			Graph::Vertex v = axis_vertices[i];
			const Vec3& axis_point = value<Vec3>(g, 
													graph_vertex_position, v);
			const uint32& axis_point_index = v.index_;

			compute_green_coordinates_on_axis_point(graph_name, 
												axis_point, axis_point_index);
		 
		}
	}

	/// @brief deform axis with MVC deformation type
	/// @param g axis
	/// @param graph_name axis name 
	/// @param graph_vertex_position positions of the axis vertices 
	/// @param axis_vertices
	void deform_axis_mvc(Graph& g, const std::string& graph_name, 
				std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
						const std::vector<GraphVertex>& axis_vertices)
	{
		MatrixWeights axis_weights = local_axis_weights_[graph_name]; 

		const std::size_t axis_vertices_length = axis_vertices.size();
		for (std::size_t i = 0; i < axis_vertices_length; i++)
		{
			GraphVertex axis_vertex = axis_vertices[i];
			uint32 axis_vertex_index = axis_vertex.index_;

			Vec3 new_position = {0.0, 0.0, 0.0};

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);

				uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

				new_position += axis_weights
								.position_(axis_vertex_index, cage_point_index) 
								* cage_point;

				return true;
			});

			value<Vec3>(g, graph_vertex_position, axis_vertex) = 
															new_position;
		}
	}

	/// @brief deform axis with Green deformation type
	/// @param g axis
	/// @param graph_name axis name 
	/// @param graph_vertex_position positions of the axis vertices 
	/// @param axis_vertices
	void deform_axis_green(Graph& g, const std::string& graph_name, 
			std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
						   const std::vector<GraphVertex>& axis_vertices)
	{
		MatrixWeights axis_weights = local_axis_weights_[graph_name]; 

		const std::size_t axis_vertices_length = axis_vertices.size();
		for (std::size_t i = 0; i < axis_vertices_length; i++)
		{
			GraphVertex axis_vertex = axis_vertices[i];
			uint32 axis_vertex_index = axis_vertex.index_;

			Vec3 new_position = {0.0, 0.0, 0.0};

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);

				uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

				new_position += axis_weights
								.position_(axis_vertex_index, cage_point_index) 
								* cage_point;

				return true;
			});

			Vec3 new_normal = {0.0, 0.0, 0.0};
			const auto sqrt8 = sqrt(8);

			for (std::size_t t = 0; t < cage_triangles_.size(); t++)
			{
				std::vector<Vec3> triangle_position(3);
				for (std::size_t i = 0; i < 3; i++)
				{
					triangle_position[i] = value<Vec3>(*global_cage_,
											global_cage_vertex_position_, 
											cage_triangles_[t].vertices_[i]);
				}

				const Vec3 t_normal = cage_triangles_[t].normal_;
				const auto t_u0 = cage_triangles_[t].edges_.first;
				const auto t_v0 = cage_triangles_[t].edges_.second;

				const auto t_u1 = 
								triangle_position[1] - triangle_position[0];
				const auto t_v1 = 
								triangle_position[2] - triangle_position[1];

				const auto area_face = (t_u0.cross(t_v0)).norm() * 0.5;
				double t_sj = sqrt((t_u1.squaredNorm()) * (t_v0.squaredNorm()) 
								- 2.0 * (t_u1.dot(t_v1)) * (t_u0.dot(t_v0)) +
						 	(t_v1.squaredNorm()) * (t_u0.squaredNorm())) /
									(sqrt8 * area_face);

				new_normal +=  
					axis_weights.normal_(axis_vertex_index, t) * t_sj * t_normal;
			}

			value<Vec3>(g, graph_vertex_position, axis_vertex) = 
													new_position + new_normal;
		}
	}

	/// @brief bind local cage with MVC deformation type 
	/// @param local_cage 
	/// @param cage_name local cage name
	/// @param local_cage_vertex_position positions of the local cage vertices
	/// @param local_cage_vertex_index indices of the local cage vertices
	void bind_local_cage_mvc(MESH& local_cage, 
							const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{

		parallel_foreach_cell(local_cage, [&](Vertex v) -> bool {
			const Vec3& local_cage_point = value<Vec3>(local_cage, 
											local_cage_vertex_position, v);
			const uint32& local_cage_point_index = 
					value<uint32>(local_cage, local_cage_vertex_index, v);

			compute_mvc_coordinates_on_local_cage_point(cage_name, 
								local_cage_point, local_cage_point_index);

			return true;
		});

	}

	/// @brief bind local cage with Green deformation type 
	/// @param local_cage 
	/// @param cage_name local cage name
	/// @param local_cage_vertex_position positions of the local cage vertices
	/// @param local_cage_vertex_index indices of the local cage vertices
	void bind_local_cage_green(MESH& local_cage, 
							const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		parallel_foreach_cell(local_cage, [&](Vertex v) -> bool {
			const Vec3& local_cage_point = value<Vec3>(local_cage, 
											local_cage_vertex_position, v);
			const uint32& local_cage_point_index = 
					value<uint32>(local_cage, local_cage_vertex_index, v);

			compute_green_coordinates_on_local_cage_point(cage_name, 
								local_cage_point, local_cage_point_index);

			return true;
		});
	}

	/// @brief deform local cage with MVC deformation type 
	/// @param local_cage 
	/// @param cage_name local cage name
	/// @param local_cage_vertex_position positions of the local cage vertices
	/// @param local_cage_vertex_index indices of the local cage vertices
	void deform_local_cage_mvc(MESH& local_cage, 
							const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		MatrixWeights local_cage_weights = local_cage_weights_[cage_name]; 

		parallel_foreach_cell(local_cage, [&](Vertex v) -> bool {
			uint32 local_cage_point_index = value<uint32>(local_cage, 
												local_cage_vertex_index, v);

			Vec3 new_position = {0.0, 0.0, 0.0};

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);

				uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

				new_position += local_cage_weights
					.position_(local_cage_point_index, cage_point_index) 
						* cage_point;

				return true;
			});

			value<Vec3>(local_cage, local_cage_vertex_position, v) = 
															new_position;
			return true;
		});
	}

	/// @brief deform local cage with Green deformation type 
	/// @param local_cage 
	/// @param cage_name local cage name
	/// @param local_cage_vertex_position positions of the local cage vertices
	/// @param local_cage_vertex_index indices of the local cage vertices
	void deform_local_cage_green(MESH& local_cage, 
								const std::string& cage_name, 
			std::shared_ptr<Attribute<Vec3>>& local_cage_vertex_position, 
			std::shared_ptr<Attribute<uint32>>& local_cage_vertex_index)
	{
		MatrixWeights local_cage_weights = local_cage_weights_[cage_name]; 

		parallel_foreach_cell(local_cage, [&](Vertex v) -> bool {
			uint32 vidx = 
				value<uint32>(local_cage, local_cage_vertex_index, v);

			Vec3 new_position = {0.0, 0.0, 0.0};

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, cv);

				uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, cv);

				new_position += local_cage_weights
					.position_(vidx, cage_point_index) * cage_point;

				return true;
			});

			Vec3 new_normal = {0.0, 0.0, 0.0};

			const auto sqrt8 = sqrt(8);

			for (std::size_t t = 0; t < cage_triangles_.size(); t++)
			{

				std::vector<Vec3> triangle_position(3);
				for (std::size_t i = 0; i < 3; i++)
				{
					triangle_position[i] = value<Vec3>(*global_cage_, 
												global_cage_vertex_position_, 
											cage_triangles_[t].vertices_[i]);
				}

				const Vec3 t_normal = cage_triangles_[t].normal_;
				const auto t_u0 = cage_triangles_[t].edges_.first;
				const auto t_v0 = cage_triangles_[t].edges_.second;

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
					local_cage_weights.normal_(vidx, t) * t_sj * t_normal;
			}

			value<Vec3>(local_cage, local_cage_vertex_position, v) = 
												new_position + new_normal;

			return true;
		});
	}
	
	/// @brief compute MVC coordinate on a point of the object 
	/// state-of-the-art method 
	/// [Mean Value Coordinates for Closed Triangular Meshes, Ju et al. 2005]
	/// @param surface_point 
	/// @param surface_point_index 
	/// @return boolean 
	bool compute_mvc_coordinates_on_point(const Vec3& surface_point, 
										const uint32& surface_point_index)
	{

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_global_cage_weights_;

		w_global_cage_weights_.resize(nbv_cage);
		w_global_cage_weights_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, 
											global_cage_vertex_position_, v);

			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, v);

			d[cage_point_index] = (surface_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				object_weights_
					.position_(surface_point_index, cage_point_index) = 1.0;
				return true;
			}

			u[cage_point_index] = 
				(cage_point - surface_point) / d[cage_point_index];

			return true;
		});

		double l[3], theta[3], w[3], c[3], s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{

			std::vector<uint32> triangle_index = cage_triangles_[t].indices_; 

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
					w[i] = 
						sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_global_cage_weights_[triangle_index[0]] = w[0] / sumWeights;
				w_global_cage_weights_[triangle_index[1]] = w[1] / sumWeights;
				w_global_cage_weights_[triangle_index[2]] = w[2] / sumWeights;

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
				w[i] = (theta[i] - c[(i + 1) % 3] 
						* theta[(i + 2) % 3] - c[(i + 2) % 3] 
						* theta[(i + 1) % 3]) /
					   	(2.0 * d[triangle_index[i]] 
						* sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_global_cage_weights_[triangle_index[0]] += w[0];
			w_global_cage_weights_[triangle_index[1]] += w[1];
			w_global_cage_weights_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, v);

			object_weights_.position_(surface_point_index, cage_point_index) =
				w_global_cage_weights_[cage_point_index] / sumWeights;

			return true;
		});

		return false;
	}

	/// @brief compute MVC coordinate on the handle 
	/// state-of-the-art method 
	/// [Mean Value Coordinates for Closed Triangular Meshes, Ju et al. 2005]
	/// @param graph_name handle name 
	/// @param handle_point handle position 
	/// @return boolean 
	bool compute_mvc_coordinates_on_handle(const std::string& graph_name, 
											const Vec3& handle_point)
	{

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_global_cage_weights_;

		w_global_cage_weights_.resize(nbv_cage);
		w_global_cage_weights_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, 
										global_cage_vertex_position_, v);

			uint32 cage_point_idx = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, v);

			d[cage_point_idx] = (handle_point - cage_point).norm();
			if (d[cage_point_idx] < epsilon)
			{
				local_handle_weights_[graph_name]
									.position_[cage_point_idx] = 1.0;
				return true;
			}

			u[cage_point_idx] = 
						(cage_point - handle_point) / d[cage_point_idx];

			return true;
		});

		double l[3], theta[3], w[3], c[3], s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			std::vector<uint32> triangle_index = 
											cage_triangles_[t].indices_; 

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
					w[i] = 
						sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_global_cage_weights_[triangle_index[0]] = w[0] / sumWeights;
				w_global_cage_weights_[triangle_index[1]] = w[1] / sumWeights;
				w_global_cage_weights_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / (sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) - 1.0;
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
					(2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) 
					* s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_global_cage_weights_[triangle_index[0]] += w[0];
			w_global_cage_weights_[triangle_index[1]] += w[1];
			w_global_cage_weights_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, v);

			local_handle_weights_[graph_name]
				.position_[cage_point_index] = 
					w_global_cage_weights_[cage_point_index] / sumWeights;

			return true;
		});

		return false;
	}

	/// @brief compute MVC coordinates on a point of the axis
	/// state-of-the-art method 
	/// [Mean Value Coordinates for Closed Triangular Meshes, Ju et al. 2005]
	/// @param graph_name axis name 
	/// @param axis_point position of the target point on the axis
	/// @param axis_point_index index of the target point on the axis 
	/// @return boolean 
	bool compute_mvc_coordinates_on_axis_point(const std::string& graph_name, 
										const Vec3& axis_point, 
										const uint32& axis_point_index)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_global_cage_weights_;

		w_global_cage_weights_.resize(nbv_cage);
		w_global_cage_weights_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, 
											global_cage_vertex_position_, v);

			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, v);

			d[cage_point_index] = (axis_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				local_axis_weights_[graph_name]
					.position_(axis_point_index, cage_point_index) = 1.0;
				return true;
			}

			u[cage_point_index] = 
				(cage_point - axis_point) / d[cage_point_index];

			return true;
		});

		double l[3], theta[3], w[3], c[3], s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{

			std::vector<uint32> triangle_index = 
											cage_triangles_[t].indices_; 

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
					w[i] = 
						sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_global_cage_weights_[triangle_index[0]] = w[0] / sumWeights;
				w_global_cage_weights_[triangle_index[1]] = w[1] / sumWeights;
				w_global_cage_weights_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = 
					(2.0 * sin(h) * sin(h - theta[i])) / 
						(sin(theta[(i + 1) % 3]) * 
					sin(theta[(i + 2) % 3])) 
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
				w[i] = (theta[i] - c[(i + 1) % 3] * 
					theta[(i + 2) % 3] - c[(i + 2) % 3] * 
						theta[(i + 1) % 3]) /
					(2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) *
					 	s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_global_cage_weights_[triangle_index[0]] += w[0];
			w_global_cage_weights_[triangle_index[1]] += w[1];
			w_global_cage_weights_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, v);

			local_axis_weights_[graph_name]
				.position_(axis_point_index, cage_point_index) = 
					w_global_cage_weights_[cage_point_index] / sumWeights;
			return true;
		});

		return false;
	}

	/// @brief compute MVC coordinates on a point of the local cage
	/// state-of-the-art method 
	/// [Mean Value Coordinates for Closed Triangular Meshes, Ju et al. 2005]
	/// @param graph_name cage name 
	/// @param local_cage_point position of the target point on the cage
	/// @param local_cage_point_index index of the target point on the cage 
	/// @return boolean 
	bool compute_mvc_coordinates_on_local_cage_point(
				const std::string& cage_name, 
				const Vec3& local_cage_point, 
				const uint32& local_cage_point_index)
	{

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_global_cage_weights_;

		w_global_cage_weights_.resize(nbv_cage);
		w_global_cage_weights_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, 
											global_cage_vertex_position_, v);

			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, v);

			d[cage_point_index] = (local_cage_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				local_cage_weights_[cage_name]
					.position_(local_cage_point_index, cage_point_index) = 1.0;
				return true;
			}

			u[cage_point_index] = 
				(cage_point - local_cage_point) / d[cage_point_index];

			return true;
		});

		double l[3];
		double theta[3];
		double w[3];
		double c[3];
		double s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{

			std::vector<uint32> triangle_index = 
				cage_triangles_[t].indices_; 

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
					w[i] = 
						sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_global_cage_weights_[triangle_index[0]] = w[0] / sumWeights;
				w_global_cage_weights_[triangle_index[1]] = w[1] / sumWeights;
				w_global_cage_weights_[triangle_index[2]] = w[2] / sumWeights;

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
				s[i] = sign_Basis_u0u1u2 * sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
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
					(2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) * 
						s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_global_cage_weights_[triangle_index[0]] += w[0];
			w_global_cage_weights_[triangle_index[1]] += w[1];
			w_global_cage_weights_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_index = value<uint32>(*global_cage_, 
											global_cage_vertex_index_, v);

			local_cage_weights_[cage_name]
				.position_(local_cage_point_index, cage_point_index) =
					w_global_cage_weights_[cage_point_index] / sumWeights;

			return true;
		});

		return false;
	}

	/// @brief compute Green coordinates on a point of the object
	/// state-of-the-art method 
	/// [Green coordinates, Lipman et al. 2008]
	/// @param surface_point position of the target point on the object
	/// @param surface_point_index index of the target point on the object 
	/// @return boolean 
	void compute_green_coordinates_on_point(const Vec3& surface_point, 
										const uint32& surface_point_index)
	{

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			const std::vector<Vec3> triangle_position = 
											cage_triangles_[t].positions_; 

			const Vec3 t_normal = cage_triangles_[t].normal_; 

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
								cage_triangles_[t].indices_[l]; 

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

	/// @brief compute Green coordinates on the target handle
	/// state-of-the-art method 
	/// [Green coordinates, Lipman et al. 2008]
	/// @param graph_name handle name 
	/// @param handle_point position of the handle 
	/// @return boolean 
	void compute_green_coordinates_on_handle(const std::string& graph_name, 
											const Vec3& handle_point)
	{

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			const std::vector<Vec3> triangle_position = 
											cage_triangles_[t].positions_; 

			const Vec3 t_normal = cage_triangles_[t].normal_; 

			std::vector<Vec3> t_vj(3);
			for (std::size_t l = 0; l < 3; ++l)
			{
				t_vj[l] = triangle_position[l] - handle_point;
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

			const auto t_I_ = -abs(t_s.dot(t_I));

			local_handle_weights_[graph_name].normal_[t] = -t_I_;

			Vec3 t_w = t_I_ * t_normal;

			for (std::size_t k = 0; k < 3; ++k)
			{
				t_w += (t_II[k] * t_N[k]);
			}

			if (t_w.norm() > DBL_EPSILON)
			{
				for (size_t l = 0; l < 3; l++)
				{
					const uint32 cage_vertex_idx = 
						cage_triangles_[t].indices_[l]; 

					const Vec3 Nl = t_N[(l + 1) % 3];
					const double num = Nl.dot(t_w);
					const double denom = Nl.dot(t_vj[l]);

					local_handle_weights_[graph_name]
						.position_[cage_vertex_idx] 
							+= num / denom;
				}
			}
		}
	}

	/// @brief compute Green coordinates on a point of the axis
	/// state-of-the-art method 
	/// [Green coordinates, Lipman et al. 2008]
	/// @param graph_name axis name 
	/// @param axis_point position of the target point on the axis
	/// @param axis_point_index index of the target point on the axis 
	/// @return boolean 
	void compute_green_coordinates_on_axis_point(
						const std::string& graph_name, 
						const Vec3& axis_point, 
						const uint32& axis_point_index)
	{

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			const std::vector<Vec3> triangle_position = 
											cage_triangles_[t].positions_; 

			const Vec3 t_normal = cage_triangles_[t].normal_; 

			std::vector<Vec3> t_vj(3);
			for (std::size_t l = 0; l < 3; ++l)
			{
				t_vj[l] = triangle_position[l] - axis_point;
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

			const auto t_I_ = -abs(t_s.dot(t_I));

			local_axis_weights_[graph_name]
							.normal_(axis_point_index, t) = -t_I_;

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
									cage_triangles_[t].indices_[l]; 

					const Vec3 Nl = t_N[(l + 1) % 3];
					const double num = Nl.dot(t_w);
					const double denom = Nl.dot(t_vj[l]);

					local_axis_weights_[graph_name]
						.position_(axis_point_index, cage_point_index) 
							+= num / denom;
				}
			}
		}
	}

	/// @brief compute Green coordinates on a point of the local cage
	/// state-of-the-art method 
	/// [Green coordinates, Lipman et al. 2008]
	/// @param cage_name local cage name 
	/// @param local_cage_point position of the target point on the cage
	/// @param local_cage_point_index index of the target point on the cage 
	/// @return boolean 
	void compute_green_coordinates_on_local_cage_point(
						const std::string& cage_name, 
						const Vec3& local_cage_point, 
						const uint32& local_cage_point_index)
	{

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			const std::vector<Vec3> triangle_position = 
									cage_triangles_[t].positions_; 

			const Vec3 t_normal = cage_triangles_[t].normal_; 

			std::vector<Vec3> t_vj(3);
			for (std::size_t l = 0; l < 3; ++l)
			{
				t_vj[l] = triangle_position[l] - local_cage_point;
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

			const auto t_I_ = -abs(t_s.dot(t_I));

			local_cage_weights_[cage_name]
				.normal_(local_cage_point_index, t) = -t_I_;

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
						cage_triangles_[t].indices_[l]; 

					const Vec3 Nl = t_N[(l + 1) % 3];
					const double num = Nl.dot(t_w);
					const double denom = Nl.dot(t_vj[l]);

					local_cage_weights_[cage_name]
						.position_(local_cage_point_index, cage_point_index) 
							+= num / denom;
				}
			}
		}
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_GLOBAL_CAGE_DEFORMATION_TOOL_H_
