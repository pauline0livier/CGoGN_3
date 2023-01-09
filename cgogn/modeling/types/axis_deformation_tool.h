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
 * Inc., 51 Franklin Street, Fifth Floor, Bostraph_provider_->emit_connectivity_changed(*handle);
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_position.get());
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_radius.get());                                  *
 *******************************************************************************/

#ifndef CGOGN_MODELING_AXIS_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_AXIS_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cells_set.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>

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
	cgogn::ui::CellsSet<MESH, MeshVertex>* influence_area_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> axis_weights; 

	std::string deformation_type_;

	std::shared_ptr<Graph::Attribute<Vec3>> control_axis_vertex_position_;
	std::shared_ptr<boost::synapse::connection> axis_attribute_update_connection_;

	Eigen::MatrixXf global_cage_weights_;

	AxisDeformationTool() : control_axis_vertex_position_(nullptr), influence_area_(nullptr)
	{
	}

	~AxisDeformationTool()
	{
	}

	void create_space_tool(Graph* g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius,
						   const std::vector<Vec3>& vertex_coords)
	{
		control_axis_ = g;
		axis_skeleton_ = cgogn::modeling::create_axis(*g, vertex_position, vertex_radius, vertex_coords);

		control_axis_vertex_position_ = cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		std::shared_ptr<Graph::Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, Graph::Vertex>(*control_axis_, "vertex_index");

		cgogn::modeling::set_graph_attribute_vertex_index(*control_axis_, vertex_index.get());
	}

	void set_influence_area(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position,
							cgogn::ui::CellsSet<MESH, MeshVertex>* influence_set)
	{

		influence_set->foreach_cell([&](MeshVertex v) -> bool {
			influence_area_->select(v);
			return true;
		});
	}

	void set_geodesic_distance(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		geodesic_distance(object, vertex_position.get());
	}

	void set_binding_rigid(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		axis_weights_.resize(nbv_object, axis_skeleton_.size());
		axis_weights_.setZero();

		compute_weights(object, vertex_position);
	}

	void set_binding_loose(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		axis_weights_.resize(nbv_object, axis_skeleton_.size());
		axis_weights_.setZero();

		compute_weights(object, vertex_position);
	}

	void update_deformation_object(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			cgogn::get_attribute<uint32, MeshVertex>(object, "vertex_index");

		const Vec3 new_deformation = get_axis_deformation();

		influence_area_->foreach_cell([&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(object, object_vertex_index, v);

			value<Vec3>(object, object_vertex_position, v) += attenuation_[vidx] * new_deformation;

			return true;
		});
	}

	
private:
	std::vector<Graph::Vertex> axis_skeleton_;
	//Graph::Vertex axis_center_; 
	Vec3 axis_normal_;
	Eigen::Matrix3d local_frame_;
	Eigen::Matrix3d local_frame_inverse_;

	void compute_weights(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		std::shared_ptr<Attribute<Scalar>> vertex_geodesic_distance =
			cgogn::get_attribute<Scalar, MeshVertex>(object, "geodesic_distance");

		// float h = 0.0f;
		Scalar max_dist = 0.0;
		std::vector<Vec2> attenuation_points;

		influence_area_->foreach_cell([&](MeshVertex v) {
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			Vec3 surface_point = value<Vec3>(object, vertex_position, v);

			Vec3 point_to_handle = (surface_point - handle_position_);
			// double i_dist = point_to_handle.squaredNorm();
			// double i_dist = point_to_handle.norm();
			Scalar i_dist = value<Scalar>(object, vertex_geodesic_distance, v);
			// this->cage_influence_distance(surface_point_idx, nbf_cage, nbv_cage);

			if (i_dist > max_dist)
			{
				max_dist = i_dist;
			}

			// this->attenuation_(surface_point_idx) = (float)sin(0.5*M_PI * (i_dist ));
			attenuation_(surface_point_idx) = i_dist;
		});

		influence_area_->foreach_cell([&](MeshVertex v) {
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			// attenuation_(surface_point_idx) = 1.0 - pow((attenuation_(surface_point_idx) / max_dist),2);

			// attenuation_(surface_point_idx) = pow(1.0 - (attenuation_(surface_point_idx) / max_dist),5);

			attenuation_(surface_point_idx) = 1.0 - (attenuation_(surface_point_idx) / max_dist);
		});
	}



	void geodesic_distance(MESH& m, const Attribute<Vec3>* vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> vertex_index = cgogn::get_attribute<uint32, MeshVertex>(m, "vertex_index");

		std::shared_ptr<Attribute<Scalar>> vertex_geodesic_distance =
			cgogn::get_attribute<Scalar, MeshVertex>(m, "geodesic_distance");

		uint32 nb_vertices = nb_cells<MeshVertex>(m);

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lc =
			cgogn::geometry::cotan_operator_matrix(m, vertex_index.get(), vertex_position);

		auto vertex_area = add_attribute<Scalar, MeshVertex>(m, "__vertex_area");
		cgogn::geometry::compute_area<MeshVertex>(m, vertex_position, vertex_area.get());

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
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> heat_solver(Am - t * Lc);
		Eigen::VectorXd u = heat_solver.solve(u0);

		auto vertex_heat = get_or_add_attribute<Scalar, MeshVertex>(m, "__vertex_heat");
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			value<Scalar>(m, vertex_heat, v) = u(vidx);
			return true;
		});

		auto face_heat_gradient = get_or_add_attribute<Vec3, MeshFace>(m, "__face_heat_gradient");
		parallel_foreach_cell(m, [&](MeshFace f) -> bool {
			Vec3 g(0, 0, 0);
			Vec3 n = cgogn::geometry::normal(m, f, vertex_position);
			Scalar a = cgogn::geometry::area(m, f, vertex_position);
			std::vector<MeshVertex> vertices = incident_vertices(m, f);
			for (uint32 i = 0; i < vertices.size(); ++i)
			{
				Vec3 e = value<Vec3>(m, vertex_position, vertices[(i + 2) % vertices.size()]) -
						 value<Vec3>(m, vertex_position, vertices[(i + 1) % vertices.size()]);
				g += value<Scalar>(m, vertex_heat, vertices[i]) * n.cross(e);
			}
			g /= 2 * a;
			value<Vec3>(m, face_heat_gradient, f) = -1.0 * g.normalized();
			return true;
		});

		auto vertex_heat_gradient_div = get_or_add_attribute<Scalar, MeshVertex>(m, "__vertex_heat_gradient_div");
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			Scalar d = vertex_gradient_divergence(m, v, face_heat_gradient.get(), vertex_position);
			value<Scalar>(m, vertex_heat_gradient_div, v) = d;
			return true;
		});

		Eigen::VectorXd b(nb_vertices);
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			b(vidx) = value<Scalar>(m, vertex_heat_gradient_div, v);
			return true;
		});

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> poisson_solver(Lc);
		Eigen::VectorXd dist = poisson_solver.solve(b);

		Scalar min = dist.minCoeff();
		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			value<Scalar>(m, vertex_geodesic_distance, v) = dist(vidx) - min;
			return true;
		});

		remove_attribute<MeshVertex>(m, vertex_area);
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
