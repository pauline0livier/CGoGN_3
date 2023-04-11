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

#include <cgogn/core/types/cells_set.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>
#include <cgogn/geometry/algos/laplacian.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
/**
 * @Class Handle Deformation Tool
 * Represents the specificities of the handle tool 
*/
class HandleDeformationTool
{
	using Graph = cgogn::IncidenceGraph;

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using MeshVertex = typename mesh_traits<MESH>::Vertex;
	using MeshFace = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Influence_area_vertex
	{
		MeshVertex vertex; 
		Vec3 local_translation; 
		Vec3 max_local_translation; 
		bool shared; 

		bool reset;  
	};

public:

	Graph* control_handle_;
	Eigen::VectorXd object_weights_;

	std::shared_ptr<Graph::Attribute<Vec3>> 
									control_handle_vertex_position_;
	std::shared_ptr<boost::synapse::connection> 
									handle_attribute_update_connection_;

	Eigen::VectorXf global_cage_weights_;
	Eigen::VectorXf global_cage_normal_weights_;

	std::unordered_map<uint32, Influence_area_vertex> object_influence_area_; 

	std::unordered_map<uint32, MeshVertex> shared_vertex_; 

	std::string deformation_type_;

	Scalar radius_of_influence_; 

	uint32 handle_mesh_vertex_index_;

	bool need_full_bind_;  

	HandleDeformationTool(): control_handle_vertex_position_(nullptr), need_full_bind_(false)
	{
	}

	~HandleDeformationTool()
	{
	}

	/// @brief create handle space tool
	/// @param g default graph 
	/// @param g_vertex_position position of the vertices of the default graph
	/// @param g_vertex_radius radius of the vertices of the default graph 
	/// @param center handle position 
	/// @param handle_name handle name
	void create_space_tool(Graph* g, 
						Graph::Attribute<Vec3>* g_vertex_position, 
						Graph::Attribute<Scalar>* g_vertex_radius,
						const Vec3& center, 
						const std::string& handle_name)
	{
		control_handle_ = g;
		handle_vertex_ = cgogn::modeling::create_handle(*g, g_vertex_position, 
											g_vertex_radius, center, Scalar(5));
											// 3 for low poly fox, 0.25 sphere

		control_handle_vertex_position_ = 
				cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		std::shared_ptr<Graph::Attribute<uint32>> g_vertex_index =
				cgogn::add_attribute<uint32, Graph::Vertex>(*control_handle_, 
															"vertex_index");
		cgogn::modeling::set_attribute_vertex_index_graph(*control_handle_, 
														g_vertex_index.get());

		handle_position_ = center;
		handle_name_ = handle_name; 

		start_position_ = center; 
	}

	/// @brief set object influence area
	/// create set of influence area vertex 
	void set_object_influence_area(MESH& object, 
			CMap2::Attribute<uint32>* object_vertex_index,
			const std::vector<MeshVertex>& influence_set)
	{
		for (std::size_t i = 0; i < influence_set.size(); i++)
		{
			MeshVertex v = influence_set[i]; 
			uint32 vertex_index = 
					value<uint32>(object, object_vertex_index, v); 

			Influence_area_vertex new_element; 
			new_element.vertex = influence_set[i];  
			new_element.local_translation = {0.0, 0.0, 0.0}; 
			new_element.max_local_translation = {0.0, 0.0, 0.0}; 
			new_element.shared = false; 
			new_element.reset = false; 
			
			object_influence_area_[vertex_index] = new_element; 
		}
	}

	/// @brief set deformation type
	/// choice between Round and Spike
	/// @param new_type 
	void set_deformation_type(const std::string new_type){
		deformation_type_ = new_type; 
	}

	/// @brief reset deformation
	void reset_deformation(MESH& object, 
					CMap2::Attribute<Vec3>* object_vertex_position,
					std::unordered_map<uint32,modeling::SharedVertexData>& activation_map, std::unordered_map<std::string, 
		std::shared_ptr<modeling::HandleDeformationTool<MESH>>>& handle_container)
	{
		value<Vec3>(*control_handle_, 
							control_handle_vertex_position_, handle_vertex_) = start_position_; 
		handle_position_ = get_handle_position(); 

		for ( const auto &myPair : object_influence_area_ ) {
			uint32 vertex_index = myPair.first;
			MeshVertex v = myPair.second.vertex; 

			if (!object_influence_area_[vertex_index].shared){
				value<Vec3>(object, object_vertex_position, v) -= object_influence_area_[vertex_index].max_local_translation;
			} else {
				if (activation_map[vertex_index].current_max_handle_.first == handle_name_){
					value<Vec3>(object, object_vertex_position, v) -= object_influence_area_[vertex_index].max_local_translation;

					std::string name_next_max_handle; 
					double current_max = 0; 
					Vec3 current_max_vector = {0.0, 0.0, 0.0}; 
					for ( const auto &otherPair : activation_map[vertex_index].handle_translation_ )
						{
							if (otherPair.first != handle_name_)
							{
								std::shared_ptr<modeling::HandleDeformationTool<MESH>> other_hdt = handle_container[otherPair.first]; 

								double local_max = other_hdt->object_influence_area_[vertex_index].max_local_translation.squaredNorm(); 

								if (local_max > current_max)
								{
									current_max = local_max; 
									current_max_vector = other_hdt->object_influence_area_[vertex_index].max_local_translation; 

									name_next_max_handle = otherPair.first; 
								}
							}
						}
						value<Vec3>(object, object_vertex_position, v) += current_max_vector; 

						activation_map[vertex_index].current_max_handle_ = make_pair(name_next_max_handle, current_max); 
				} 
			}
			object_influence_area_[vertex_index].max_local_translation = {0.0, 0.0, 0.0}; 

			
		}


	}

	/// @brief require full re-binding of the handle
	/// happen when the handle position is updated
	/// outside of the handle deformation
	void require_full_binding()
	{
		need_full_bind_ = true; 
	}


	/// @brief reset transformation
	Vec3 get_reset_transformation()
	{ 
		return start_position_ - value<Vec3>(*control_handle_, 
							control_handle_vertex_position_, handle_vertex_); 
	}

	/// @brief set geodesic distance of the points 
	/// belonging to the zone of influence 
	/// @param object 
	/// @param object_vertex_position 
	void set_geodesic_distance(MESH& object, 
			const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{
		geodesic_distance(object, object_vertex_position.get());

	}


	/// @brief update handle position 
	/// useful when handle is displaced by other spatial tools 
	/// @param new_position 
	void update_handle_position_variable()
	{

		handle_position_ = get_handle_position(); 
		start_position_ = get_handle_position(); 
	}

	/// @brief update handle position 
	/// used when other handles around 
	void update_handle_position(const Vec3& new_position)
	{ 
		handle_position_ = new_position; 
		value<Vec3>(*control_handle_, 
				control_handle_vertex_position_, handle_vertex_) = new_position; 

	}

	/// @brief 
	/// @return handle vertex 
	Graph::Vertex get_handle_vertex(){
		return handle_vertex_; 
	}

	/// @brief useful to compute the handle displacement 
	/// or the weights with other spatial tools 
	/// @return handle position 
	Vec3 get_handle_position(){
		return value<Vec3>(*control_handle_, 
							control_handle_vertex_position_, handle_vertex_);
	}

	/// @brief initialize binding for object
	/// @param object model to deform
	/// @param object_vertex_position position of the vertices of the model 
	void init_bind_object(MESH& object)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		object_weights_.resize(nbv_object);
		object_weights_.setZero();

		Scalar current_max = 0; 
		for ( const auto &myPair : object_influence_area_ ) {
			uint32 vertex_index = myPair.first;

			if (geodesic_distance_[vertex_index] > current_max){
				current_max = geodesic_distance_[vertex_index]; 
			}
		}
 
		radius_of_influence_ = current_max; 

		bind_object_round();
	}

	/// @brief binding for object
	/// @param object model to deform
	/// @param object_vertex_position position of the vertices of the model 
	void bind_object()
	{
		object_weights_.setZero();

		bind_object_round();
	}

	/// @brief rebind required if the handle is displaced outside of its deformation or if a vertex is moved by another tool 
	/// @param object 
	/// @param object_vertex_position 
	void rebind_object(MESH& object, 
						CMap2::Attribute<Vec3>* object_vertex_position)
	{
		geodesic_distance(object, object_vertex_position);
		object_weights_.setZero();

		bind_object_round();

	}

	void bind_isolated_vertex(const uint32& vertex_index)
	{
		object_weights_[vertex_index] = 
				exp(-(geodesic_distance_[vertex_index]*geodesic_distance_[vertex_index]) / 
								(radius_of_influence_*radius_of_influence_));
	}

	/// @brief deform the object 
	/// for each point of the zone of influence 
	/// new_position_i += weights_i*deformation
	/// @param object 
	/// @param object_vertex_position 
	/// @param object_vertex_index 
	void deform_object(MESH& object, 
					CMap2::Attribute<Vec3>* object_vertex_position, 
					CMap2::Attribute<uint32>* object_vertex_index )
	{

		if (need_full_bind_)
		{
			rebind_object(object, object_vertex_position); 
			need_full_bind_ = false;
		}

		const Vec3 new_deformation = get_handle_deformation();

		for ( const auto &myPair : object_influence_area_ ) {
			uint32 vertex_index = myPair.first;
			MeshVertex v = myPair.second.vertex; 

			const Vec3 new_transformation = 
						object_weights_[vertex_index] * new_deformation; 

			object_influence_area_[vertex_index].max_local_translation += 
					new_transformation; 
			
			object_influence_area_[vertex_index].local_translation = 
					new_transformation; 
 
			if (!object_influence_area_[vertex_index].shared)
			{
				value<Vec3>(object, object_vertex_position, v) += 
					new_transformation;
			} 
			
		}
	}

	/// @brief set handle mesh vertex 
	/// vertex of the object associated to the handle 
	/// useful to compute geodesic distance 
	/// @param m_v mesh vertex
	void set_handle_mesh_vertex(const MeshVertex& m_v)
	{
		handle_mesh_vertex_ = m_v;
	}

	const uint32& get_handle_mesh_vertex_index(MESH& object, 
					CMap2::Attribute<uint32>* object_vertex_index)
	{
		return value<uint32>(object, 
							object_vertex_index, handle_mesh_vertex_);
	}
	
	/// @brief compute handle deformation
	/// @return handle deformation
	const Vec3 get_handle_deformation()
	{
		const Vec3 handle_new_position = value<Vec3>(*control_handle_, 
							control_handle_vertex_position_, handle_vertex_);

		const Vec3 deformation = (handle_new_position - handle_position_);

		handle_position_ = handle_new_position;

		return deformation;
	}

	

private:
	Graph::Vertex handle_vertex_;
	MeshVertex handle_mesh_vertex_;
	Vec3 handle_position_; 

	Vec3 start_position_;  

	std::string handle_name_; 

	std::unordered_map<uint32, Scalar> geodesic_distance_; 

	/// @brief bind object with round deformation type 
	/// geodesic distance between object point and handle 
	/// attenuation function exp(-d²/d²_max), d = distance
	/// exp(-sqrt(geodesic_distance_[v]/max_dist)) for spike
	/// @param object 
	/// @param object_vertex_position 
	void bind_object_round()
	{
		for ( const auto &myPair : object_influence_area_ ) {
			uint32 vertex_index = myPair.first;

			object_weights_[vertex_index] = 
				exp(-(geodesic_distance_[vertex_index]*geodesic_distance_[vertex_index]) / 
								(radius_of_influence_*radius_of_influence_));
		}
	}

	/// @brief set geodesic distance between a mesh m and the handle position
	/// @param m mesh 
	/// @param m_vertex_position position of the vertices of the mesh  
	void geodesic_distance(MESH& m, const Attribute<Vec3>* m_vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> m_vertex_index = 
					cgogn::get_attribute<uint32, MeshVertex>(m, "vertex_index");

		uint32 nb_vertices = nb_cells<MeshVertex>(m);

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lc =
			cgogn::geometry::cotan_operator_matrix(m, 
									m_vertex_index.get(), m_vertex_position);

		auto m_vertex_area = add_attribute<Scalar, MeshVertex>(m, 
														"__vertex_area");
		cgogn::geometry::compute_area<MeshVertex>(m, m_vertex_position, 
													m_vertex_area.get());

		Eigen::VectorXd A(nb_vertices);
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, m_vertex_index, v);
			A(vidx) = value<Scalar>(m, m_vertex_area, v);
			return true;
		});

		Eigen::VectorXd u0(nb_vertices);
		u0.setZero();
		uint32 vidx = value<uint32>(m, m_vertex_index, handle_mesh_vertex_);
		handle_mesh_vertex_index_ = vidx; 
		u0(vidx) = 1.0;

		Scalar h = cgogn::geometry::mean_edge_length(m, m_vertex_position);
		Scalar t = h * h;

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Am(A.asDiagonal());
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> 
													heat_solver(Am - t * Lc);
		Eigen::VectorXd u = heat_solver.solve(u0);

		auto m_vertex_heat = get_or_add_attribute<Scalar, MeshVertex>(m, 
															"__vertex_heat");
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, m_vertex_index, v);
			value<Scalar>(m, m_vertex_heat, v) = u(vidx);
			return true;
		});

		auto m_face_heat_gradient = get_or_add_attribute<Vec3, MeshFace>(m, 
													"__face_heat_gradient");

		parallel_foreach_cell(m, [&](MeshFace f) -> bool {
			Vec3 g(0, 0, 0);
			Vec3 n = cgogn::geometry::normal(m, f, m_vertex_position);
			Scalar a = cgogn::geometry::area(m, f, m_vertex_position);
			std::vector<MeshVertex> vertices = incident_vertices(m, f);
			for (uint32 i = 0; i < vertices.size(); ++i)
			{
				Vec3 e = 
					value<Vec3>(m, m_vertex_position, vertices[(i + 2) % vertices.size()]) 
					- value<Vec3>(m, m_vertex_position, vertices[(i + 1) % vertices.size()]);
				
				g += value<Scalar>(m, m_vertex_heat, vertices[i]) * n.cross(e);
			}
			g /= 2 * a;
			value<Vec3>(m, m_face_heat_gradient, f) = -1.0 * g.normalized();
			return true;
		});

		auto m_vertex_heat_gradient_div = 
			get_or_add_attribute<Scalar, MeshVertex>(m, 
										"__vertex_heat_gradient_div");

		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			Scalar d = 
				vertex_gradient_divergence(m, v, m_face_heat_gradient.get(), 
											m_vertex_position);
			value<Scalar>(m, m_vertex_heat_gradient_div, v) = d;
			return true;
		});

		Eigen::VectorXd b(nb_vertices);
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, m_vertex_index, v);
			b(vidx) = value<Scalar>(m, m_vertex_heat_gradient_div, v);
			return true;
		});

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> 
														poisson_solver(Lc);
		Eigen::VectorXd dist = poisson_solver.solve(b);

		Scalar min = dist.minCoeff();
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, m_vertex_index, v);

			geodesic_distance_[vidx] = dist(vidx) - min;
														
			return true;
		});

		remove_attribute<MeshVertex>(m, m_vertex_area);
		remove_attribute<MeshVertex>(m,m_vertex_heat_gradient_div);  
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_
