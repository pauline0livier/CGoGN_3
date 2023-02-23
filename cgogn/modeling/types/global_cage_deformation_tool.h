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

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <boost/synapse/connect.hpp>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
/**
 * @Class Global Cage Deformation Tool
*/
class GlobalCageDeformationTool
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

	using Graph = cgogn::IncidenceGraph;

	template <typename T>
	using GraphAttribute = typename mesh_traits<IncidenceGraph>::template Attribute<T>;
	using GraphVertex = IncidenceGraph::Vertex;
	using GraphEdge = IncidenceGraph::Edge;
	using GraphFace = IncidenceGraph::Face;

public:
	MESH* global_cage_;
	std::shared_ptr<Attribute<Vec3>> global_cage_vertex_position_;

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;

	std::vector<std::vector<CMap2::Vertex>> cage_triangles_;
	std::vector<Vec3> cage_triangles_normal_;
	std::vector<std::pair<Vec3, Vec3>> cage_triangles_edge_;

	std::string deformation_type_;

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> global_cage_coords_;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> global_cage_normal_coords_;

	GlobalCageDeformationTool() : global_cage_vertex_position_(nullptr)
	{
	}

	~GlobalCageDeformationTool()
	{
	}

	void create_global_cage(MESH* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max)
	{
		global_cage_ = m;

		vertices_ = cgogn::modeling::create_bounding_box(*m, vertex_position, bb_min, bb_max);

		global_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		std::shared_ptr<Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, Vertex>(*global_cage_, "vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*global_cage_, vertex_index.get());

		foreach_cell(*global_cage_, [&](Face fc) -> bool {
			std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*global_cage_, fc);

			// triangle 1
			std::vector<CMap2::Vertex> triangle1_vertex(3);
			triangle1_vertex[0] = face_vertices_[1];
			triangle1_vertex[1] = face_vertices_[3];
			triangle1_vertex[2] = face_vertices_[0];

			cage_triangles_.push_back(triangle1_vertex);

			std::vector<CMap2::Vertex> triangle2_vertex(3);
			triangle2_vertex[0] = face_vertices_[1];
			triangle2_vertex[1] = face_vertices_[2];
			triangle2_vertex[2] = face_vertices_[3];

			cage_triangles_.push_back(triangle2_vertex);

			return true;
		});
	}

	void update_global_cage(const Vec3& bb_min, const Vec3& bb_max)
	{
		cgogn::modeling::update_bounding_box(*global_cage_, global_cage_vertex_position_.get(), vertices_, bb_min,
											 bb_max);
	}

	void bind_mvc(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position)
	{
		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		global_cage_coords_.resize(nbv_object, nbv_cage);
		global_cage_coords_.setZero();

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			const uint32& surface_point_idx = value<uint32>(object, object_vertex_index, v);

			compute_mvc_coordinates_on_point(surface_point, surface_point_idx);

			return true;
		});
	}

	void update_local_mvc(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position,
						  cgogn::ui::CellsSet<MESH, Vertex>* influence_area)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		influence_area->foreach_cell([&](Vertex v) {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			compute_mvc_coordinates_on_point(surface_point, surface_point_idx);
		});
	}

	Eigen::VectorXf bind_mvc_handle(Graph& g, std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
									Graph::Vertex handle_vertex)
	{

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		std::shared_ptr<GraphAttribute<uint32>> graph_vertex_index =
			get_attribute<uint32, GraphVertex>(g, "vertex_index");

		Eigen::VectorXf handle_weights;

		handle_weights.resize(nbv_cage);
		handle_weights.setZero();

		const Vec3& graph_point = value<Vec3>(g, graph_vertex_position, handle_vertex);

		compute_mvc_coordinates_on_handle(graph_point, handle_weights);

		return handle_weights;
	}

	void bind_green(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);
		uint32 nbt_cage = cage_triangles_.size();

		global_cage_coords_.resize(nbv_object, nbv_cage);
		global_cage_coords_.setZero();

		global_cage_normal_coords_.resize(nbv_object, nbt_cage);
		global_cage_normal_coords_.setZero();

		cage_triangles_normal_.resize(nbt_cage);
		cage_triangles_edge_.resize(nbt_cage);

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);

			const uint32& surface_point_idx = value<uint32>(object, object_vertex_index, v);

			compute_green_coordinates_on_point(surface_point, surface_point_idx);

			return true;
		});
	}

	std::pair<Eigen::VectorXf, Eigen::VectorXf> bind_green_handle(Graph& g, std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
						   Graph::Vertex handle_vertex)
	{

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);
		uint32 nbt_cage = cage_triangles_.size();

		Eigen::VectorXf handle_weights;
		Eigen::VectorXf handle_normal_weights;

		handle_weights.resize(nbv_cage);
		handle_weights.setZero();

		handle_normal_weights.resize(nbt_cage);
		handle_normal_weights.setZero();

		// check if already computed or not => perhaps to compute by default
		cage_triangles_normal_.resize(nbt_cage);
		cage_triangles_edge_.resize(nbt_cage);

		const Vec3& graph_point = value<Vec3>(g, graph_vertex_position, handle_vertex);

		compute_green_coordinates_on_handle(graph_point, handle_weights, handle_normal_weights);

		return std::make_pair(handle_weights, handle_normal_weights);
	}

	void update_mvc(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(object, object_vertex_index, v);

			Vec3 new_pos_ = {0.0, 0.0, 0.0};

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, global_cage_vertex_position_, cv);

				uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, cv);

				new_pos_ += global_cage_coords_(vidx, cage_point_idx) * cage_point;

				return true;
			});

			value<Vec3>(object, object_vertex_position, v) = new_pos_;
			return true;
		});
	}

	void update_mvc_handle(Graph& g, std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position,
						   Eigen::VectorXf weights)
	{

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		foreach_cell(g, [&](GraphVertex v) -> bool {
			Vec3 new_pos_ = {0.0, 0.0, 0.0};

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, global_cage_vertex_position_, cv);

				uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, cv);

				new_pos_ += weights(cage_point_idx) * cage_point;

				return true;
			});

			value<Vec3>(g, graph_vertex_position, v) = new_pos_;
			return true;
		});
	}

	void update_green(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_face_index =
			cgogn::get_attribute<uint32, Face>(*global_cage_, "face_index");

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(object, object_vertex_index, v);

			Vec3 new_pos_update_ = {0.0, 0.0, 0.0};

			const auto sqrt8 = sqrt(8);

			foreach_cell(*global_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*global_cage_, global_cage_vertex_position_, cv);

				uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, cv);

				new_pos_update_ += global_cage_coords_(vidx, cage_point_idx) * cage_point;

				return true;
			});

			Vec3 new_norm_update_ = {0.0, 0.0, 0.0};

			for (std::size_t t = 0; t < cage_triangles_.size(); t++)
			{

				std::vector<Vec3> triangle_position(3);
				for (std::size_t i = 0; i < 3; i++)
				{
					triangle_position[i] =
						value<Vec3>(*global_cage_, global_cage_vertex_position_, cage_triangles_[t][i]);
				}

				const Vec3 t_normal = cage_triangles_normal_[t];
				const auto t_u0 = cage_triangles_edge_[t].first;
				const auto t_v0 = cage_triangles_edge_[t].second;

				const auto t_u1 = triangle_position[1] - triangle_position[0];
				const auto t_v1 = triangle_position[2] - triangle_position[1];

				const auto area_face = (t_u0.cross(t_v0)).norm() * 0.5;
				double t_sj =
					sqrt((t_u1.squaredNorm()) * (t_v0.squaredNorm()) - 2.0 * (t_u1.dot(t_v1)) * (t_u0.dot(t_v0)) +
						 (t_v1.squaredNorm()) * (t_u0.squaredNorm())) /
					(sqrt8 * area_face);

				new_norm_update_ += global_cage_normal_coords_(vidx, t) * t_sj * t_normal;
			}

			value<Vec3>(object, object_vertex_position, v) = new_pos_update_ + new_norm_update_;

			return true;
		});
	}

private:
	std::vector<CMap2::Vertex> vertices_;

	bool compute_mvc_coordinates_on_point(const Vec3& surface_point, const uint32& surface_point_idx)
	{

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_global_cage_coords_;

		w_global_cage_coords_.resize(nbv_cage);
		w_global_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, global_cage_vertex_position_, v);

			uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, v);

			d[cage_point_idx] = (surface_point - cage_point).norm();
			if (d[cage_point_idx] < epsilon)
			{
				global_cage_coords_(surface_point_idx, cage_point_idx) = 1.0;
				return true;
			}

			u[cage_point_idx] = (cage_point - surface_point) / d[cage_point_idx];

			return true;
		});

		double l[3];
		double theta[3];
		double w[3];
		double c[3];
		double s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{

			std::vector<uint32> triangle_index(3);
			for (std::size_t i = 0; i < 3; i++)
			{
				triangle_index[i] = value<uint32>(*global_cage_, cage_vertex_index, cage_triangles_[t][i]);
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = (u[triangle_index[(i + 1) % 3]] - u[triangle_index[(i + 2) % 3]]).norm();
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
				w_global_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_global_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_global_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / (sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) - 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || fabs(s[2]) < epsilon)
			{
				continue; // eta is on the same plane, outside t  ->  ignore triangle t :
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = (theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] - c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					   (2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_global_cage_coords_[triangle_index[0]] += w[0];
			w_global_cage_coords_[triangle_index[1]] += w[1];
			w_global_cage_coords_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, v);

			global_cage_coords_(surface_point_idx, cage_point_idx) = w_global_cage_coords_[cage_point_idx] / sumWeights;

			return true;
		});

		return false;
	}

	bool compute_mvc_coordinates_on_handle(const Vec3& handle_point, Eigen::VectorXf& handle_weights)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_global_cage_coords_;

		w_global_cage_coords_.resize(nbv_cage);
		w_global_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*global_cage_, global_cage_vertex_position_, v);

			uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, v);

			d[cage_point_idx] = (handle_point - cage_point).norm();
			if (d[cage_point_idx] < epsilon)
			{
				handle_weights[cage_point_idx] = 1.0;
				return true;
			}

			u[cage_point_idx] = (cage_point - handle_point) / d[cage_point_idx];

			return true;
		});

		double l[3];
		double theta[3];
		double w[3];
		double c[3];
		double s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{

			std::vector<uint32> triangle_index(3);
			for (std::size_t i = 0; i < 3; i++)
			{
				triangle_index[i] = value<uint32>(*global_cage_, cage_vertex_index, cage_triangles_[t][i]);
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = (u[triangle_index[(i + 1) % 3]] - u[triangle_index[(i + 2) % 3]]).norm();
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
				w_global_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_global_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_global_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / (sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) - 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || fabs(s[2]) < epsilon)
			{
				continue; // eta is on the same plane, outside t  ->  ignore triangle t :
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = (theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] - c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					   (2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_global_cage_coords_[triangle_index[0]] += w[0];
			w_global_cage_coords_[triangle_index[1]] += w[1];
			w_global_cage_coords_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*global_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, v);

			handle_weights[cage_point_idx] = w_global_cage_coords_[cage_point_idx] / sumWeights;

			return true;
		});

		return false;
	}

	// GREEN
	void compute_green_coordinates_on_point(const Vec3& surface_point, const uint32& surface_point_idx)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		const Vec3 NULL_VECTOR(0.0, 0.0, 0);

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			std::vector<Vec3> triangle_position(3);
			for (std::size_t i = 0; i < 3; i++)
			{
				triangle_position[i] = value<Vec3>(*global_cage_, global_cage_vertex_position_, cage_triangles_[t][i]);
			}

			const Vec3 t_normal =
				(cgogn::geometry::normal(triangle_position[0], triangle_position[1], triangle_position[2]))
					.normalized();
			cage_triangles_normal_[t] = t_normal;

			cage_triangles_edge_[t].first = triangle_position[1] - triangle_position[0];
			cage_triangles_edge_[t].second = triangle_position[2] - triangle_position[1];

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

				const auto t_vjpt = ((t_v0 - t_p_).cross((t_v1 - t_p_))).dot(t_normal);
				t_s[l] = t_vjpt < 0 ? -1.0 : 1.0;

				t_I[l] = cgogn::modeling::GCTriInt2(t_p_, t_v0, t_v1);
				t_II[l] = cgogn::modeling::GCTriInt2(NULL_VECTOR, t_v1, t_v0);
				t_N[l] = (t_v1.cross(t_v0)).normalized();
			}

			const auto t_I_ = -abs(t_s.dot(t_I));

			global_cage_normal_coords_(surface_point_idx, t) = -t_I_;

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
						value<uint32>(*global_cage_, cage_vertex_index, cage_triangles_[t][l]);

					const auto Nl = t_N[(l + 1) % 3];
					const auto num = Nl.dot(t_w);
					const auto denom = Nl.dot(t_vj[l]);

					global_cage_coords_(surface_point_idx, cage_vertex_idx) =
						global_cage_coords_(surface_point_idx, cage_vertex_idx) + num / denom;
				}
			}
		}
	}

	void compute_green_coordinates_on_handle(const Vec3& handle_point, Eigen::VectorXf& handle_weights,
											 Eigen::VectorXf& handle_normal_weights)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		const Vec3 NULL_VECTOR(0.0, 0.0, 0);

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{
			std::vector<Vec3> triangle_position(3);
			for (std::size_t i = 0; i < 3; i++)
			{
				triangle_position[i] = value<Vec3>(*global_cage_, global_cage_vertex_position_, cage_triangles_[t][i]);
			}

			const Vec3 t_normal =
				(cgogn::geometry::normal(triangle_position[0], triangle_position[1], triangle_position[2]))
					.normalized();
			cage_triangles_normal_[t] = t_normal;

			cage_triangles_edge_[t].first = triangle_position[1] - triangle_position[0];
			cage_triangles_edge_[t].second = triangle_position[2] - triangle_position[1];

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

				const auto t_vjpt = ((t_v0 - t_p_).cross((t_v1 - t_p_))).dot(t_normal);
				t_s[l] = t_vjpt < 0 ? -1.0 : 1.0;

				t_I[l] = cgogn::modeling::GCTriInt2(t_p_, t_v0, t_v1);
				t_II[l] = cgogn::modeling::GCTriInt2(NULL_VECTOR, t_v1, t_v0);
				t_N[l] = (t_v1.cross(t_v0)).normalized();
			}

			const auto t_I_ = -abs(t_s.dot(t_I));

			handle_normal_weights[t] = -t_I_;

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
						value<uint32>(*global_cage_, cage_vertex_index, cage_triangles_[t][l]);

					const auto Nl = t_N[(l + 1) % 3];
					const auto num = Nl.dot(t_w);
					const auto denom = Nl.dot(t_vj[l]);

					handle_weights[cage_vertex_idx] = handle_weights[cage_vertex_idx] + num / denom;
				}
			}
		}
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_GLOBAL_CAGE_DEFORMATION_TOOL_H_
