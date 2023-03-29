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

public:

	Graph* control_handle_;
	Eigen::VectorXd attenuation_;

	std::shared_ptr<Graph::Attribute<Vec3>> 
									control_handle_vertex_position_;
	std::shared_ptr<boost::synapse::connection> 
									handle_attribute_update_connection_;

	Eigen::VectorXf global_cage_weights_;
	Eigen::VectorXf global_cage_normal_weights_;

	std::vector<MeshVertex> object_influence_area_; 

	HandleDeformationTool(): control_handle_vertex_position_(nullptr)
	{
	}

	~HandleDeformationTool()
	{
	}

	/**
	 * create handle tool 
	 * set position and normal parameters
	*/
	void create_space_tool(Graph* g, 
						Graph::Attribute<Vec3>* vertex_position, 
						Graph::Attribute<Scalar>* vertex_radius,
						const Vec3& center, const Vec3& normal)
	{
		control_handle_ = g;
		handle_vertex_ = cgogn::modeling::create_handle(*g, vertex_position, 
											vertex_radius, center, Scalar(5));

		control_handle_vertex_position_ = 
				cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		std::shared_ptr<Graph::Attribute<uint32>> vertex_index =
				cgogn::add_attribute<uint32, Graph::Vertex>(*control_handle_, 
															"vertex_index");
		cgogn::modeling::set_attribute_vertex_index_graph(*control_handle_, 
														vertex_index.get());

		handle_normal_ = normal;
		handle_position_ = center;
	}

	/**
	 * set deformation type
	 * choice between Round and Spike
	*/
	void set_deformation_type(const std::string new_type){
		deformation_type_ = new_type; 
	}

	/**
	 * set geodesic distance 
	 * distance of points belonging to influence area to handle
	*/
	void set_geodesic_distance(MESH& object, 
			const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		geodesic_distance(object, vertex_position.get());
	}

	/**
	 * get handle vertex
	*/
	Graph::Vertex get_handle_vertex(){
		return handle_vertex_; 
	}

	/**
	 * get handle_position 
	 * useful to compute the handle_displacement 
	 * or the weights with other spatial tools
	*/
	Vec3 get_handle_position(){
		return value<Vec3>(*control_handle_, 
							control_handle_vertex_position_, handle_vertex_);
	}

	/**
	 * set attenuation or object weights for spike deformation
	*/
	void set_attenuation_spike(MESH& object, 
					const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		attenuation_.resize(nbv_object);
		attenuation_.setZero();

		compute_attenuation_spike(object, vertex_position);
	}

	/**
	 * set attenuation or object weights for round deformation
	*/
	void set_attenuation_round(MESH& object, 
					const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		attenuation_.resize(nbv_object);
		attenuation_.setZero();

		compute_attenuation_round(object, vertex_position);
	}

	/**
	 * deform object 
	*/
	void deform_object(MESH& object, 
					CMap2::Attribute<Vec3>* object_vertex_position, 
					CMap2::Attribute<uint32>* object_vertex_index )
	{
		const Vec3 new_deformation = get_handle_deformation();

		const std::size_t influence_area_length =
									 this->object_influence_area_.size(); 

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			uint32 vidx = value<uint32>(object, object_vertex_index, v);

			value<Vec3>(object, object_vertex_position, v) 
									+= attenuation_[vidx] * new_deformation;
		}
	}

	/**
	 * set handle mesh vertex 
	 * vertex on the object associated to the handle
	*/
	void set_handle_mesh_vertex(const MeshVertex& m_v)
	{
		handle_mesh_vertex_ = m_v;
	}

	/**
	 * get handle deformation
	*/
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
	Vec3 handle_normal_;

	std::string deformation_type_; 

	/**
	 * compute attenuation or object weights 
	 * for round deformation type 
	*/
	void compute_attenuation_round(MESH& object, 
			const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		std::shared_ptr<Attribute<Scalar>> vertex_geodesic_distance =
			cgogn::get_attribute<Scalar, MeshVertex>(object, 
												"geodesic_distance");

		Scalar max_dist = 0.0;
		std::vector<Vec2> attenuation_points;

		const std::size_t influence_area_length = 
									this->object_influence_area_.size(); 

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			uint32 surface_point_idx = 
							value<uint32>(object, object_vertex_index, v);

			Vec3 surface_point = value<Vec3>(object, vertex_position, v);

			Vec3 point_to_handle = (surface_point - handle_position_);
			Scalar i_dist = 
						value<Scalar>(object, vertex_geodesic_distance, v);

			if (i_dist > max_dist)
			{
				max_dist = i_dist;
			}

			attenuation_(surface_point_idx) = i_dist;
			
		}

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			uint32 surface_point_idx = 
					value<uint32>(object, object_vertex_index, v);

			attenuation_[surface_point_idx] = 
						exp(-(attenuation_[surface_point_idx]*
				attenuation_[surface_point_idx]) / (max_dist*max_dist));
			
		}
	}

	/**
	 * compute attenuation or object weights
	 * for spike deformation type 
	*/
	void compute_attenuation_spike(MESH& object, 
			const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		std::shared_ptr<Attribute<Scalar>> vertex_geodesic_distance =
			cgogn::get_attribute<Scalar, MeshVertex>(object, 
												"geodesic_distance");

		Scalar max_dist = 0.0;
		std::vector<Vec2> attenuation_points;

		const std::size_t influence_area_length = 
									this->object_influence_area_.size(); 

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			uint32 surface_point_idx = 
							value<uint32>(object, object_vertex_index, v);

			Vec3 surface_point = value<Vec3>(object, vertex_position, v);

			Vec3 point_to_handle = (surface_point - handle_position_);
			Scalar i_dist = 
						value<Scalar>(object, vertex_geodesic_distance, v);

			if (i_dist > max_dist)
			{
				max_dist = i_dist;
			}

			attenuation_(surface_point_idx) = i_dist;
			
		}

		for (std::size_t i = 0; i < influence_area_length; i++)
		{
			MeshVertex v = this->object_influence_area_[i]; 
			uint32 surface_point_idx = 
					value<uint32>(object, object_vertex_index, v);
			attenuation_[surface_point_idx] = 
					exp(-sqrt(attenuation_[surface_point_idx]/max_dist)); 
			
		}
	}

	/**
	 * get geodesic distance between the mesh and the handle position 
	*/
	void geodesic_distance(MESH& m, const Attribute<Vec3>* vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> vertex_index = cgogn::get_attribute<uint32, MeshVertex>(m, "vertex_index");

		std::shared_ptr<Attribute<Scalar>> vertex_geodesic_distance =
			cgogn::get_attribute<Scalar, MeshVertex>(m, "geodesic_distance");

		uint32 nb_vertices = nb_cells<MeshVertex>(m);

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lc =
			cgogn::geometry::cotan_operator_matrix(m, 
									vertex_index.get(), vertex_position);

		auto vertex_area = add_attribute<Scalar, MeshVertex>(m, 
														"__vertex_area");
		cgogn::geometry::compute_area<MeshVertex>(m, vertex_position, 
													vertex_area.get());

		Eigen::VectorXd A(nb_vertices);
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			A(vidx) = value<Scalar>(m, vertex_area, v);
			return true;
		});

		Eigen::VectorXd u0(nb_vertices);
		u0.setZero();
		uint32 vidx = value<uint32>(m, vertex_index, handle_mesh_vertex_);
		u0(vidx) = 1.0;

		Scalar h = cgogn::geometry::mean_edge_length(m, vertex_position);
		Scalar t = h * h;

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Am(A.asDiagonal());
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> 
													heat_solver(Am - t * Lc);
		Eigen::VectorXd u = heat_solver.solve(u0);

		auto vertex_heat = get_or_add_attribute<Scalar, MeshVertex>(m, 
															"__vertex_heat");
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			value<Scalar>(m, vertex_heat, v) = u(vidx);
			return true;
		});

		auto face_heat_gradient = get_or_add_attribute<Vec3, MeshFace>(m, 
													"__face_heat_gradient");

		parallel_foreach_cell(m, [&](MeshFace f) -> bool {
			Vec3 g(0, 0, 0);
			Vec3 n = cgogn::geometry::normal(m, f, vertex_position);
			Scalar a = cgogn::geometry::area(m, f, vertex_position);
			std::vector<MeshVertex> vertices = incident_vertices(m, f);
			for (uint32 i = 0; i < vertices.size(); ++i)
			{
				Vec3 e = 
					value<Vec3>(m, vertex_position, vertices[(i + 2) % vertices.size()]) 
					- value<Vec3>(m, vertex_position, vertices[(i + 1) % vertices.size()]);
				
				g += value<Scalar>(m, vertex_heat, vertices[i]) * n.cross(e);
			}
			g /= 2 * a;
			value<Vec3>(m, face_heat_gradient, f) = -1.0 * g.normalized();
			return true;
		});

		auto vertex_heat_gradient_div = 
			get_or_add_attribute<Scalar, MeshVertex>(m, 
										"__vertex_heat_gradient_div");
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			Scalar d = 
				vertex_gradient_divergence(m, v, face_heat_gradient.get(), 
											vertex_position);
			value<Scalar>(m, vertex_heat_gradient_div, v) = d;
			return true;
		});

		Eigen::VectorXd b(nb_vertices);
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			b(vidx) = value<Scalar>(m, vertex_heat_gradient_div, v);
			return true;
		});

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> 
														poisson_solver(Lc);
		Eigen::VectorXd dist = poisson_solver.solve(b);

		Scalar min = dist.minCoeff();
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			value<Scalar>(m, vertex_geodesic_distance, v) = 
														dist(vidx) - min;
			return true;
		});

		remove_attribute<MeshVertex>(m, vertex_area);
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_
