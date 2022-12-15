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

	std::string binding_type_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_cage_coords_;
	Eigen::Matrix<Vec2, Eigen::Dynamic, Eigen::Dynamic> global_cage_normal_coords_;

	GlobalCageDeformationTool() : global_cage_vertex_position_(nullptr)
	{
	}

	~GlobalCageDeformationTool()
	{
	}

	void create_global_cage(MESH* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max)
	{
		global_cage_ = m;
		cgogn::modeling::create_bounding_box(*m, vertex_position, bb_min, bb_max);

		global_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		std::shared_ptr<Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, Vertex>(*global_cage_, "vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*global_cage_, vertex_index.get());

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, Vertex>(*global_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*global_cage_, marked_vertices.get());
	}

	void bind_mvc(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position)
	{

		uint32 nbv_object = nb_cells<Vertex>(object);

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked =
			get_attribute<bool, Vertex>(*global_cage_, "marked_vertices");

		global_cage_coords_.resize(nbv_object, nbv_cage);

		foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			DartMarker dm(*global_cage_);

			float sumMVC = 0.0;
			for (Dart d = global_cage_->begin(), end = global_cage_->end(); d != end; d = global_cage_->next(d))
			{
				Vertex cage_vertex = CMap2::Vertex(d);

				bool vc_marked = value<bool>(*global_cage_, cage_vertex_marked, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{

					const Vec3& cage_point = value<Vec3>(*global_cage_, global_cage_vertex_position_, cage_vertex);
					uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, cage_vertex);

					float mvc_value = cgogn::modeling::compute_mvc(surface_point, d, *global_cage_, cage_point,
																   global_cage_vertex_position_.get());

					global_cage_coords_(surface_point_idx, cage_point_idx) = mvc_value;

					dm.mark(d);

					value<bool>(*global_cage_, cage_vertex_marked, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			float sum_lambda = 0.0;

			parallel_foreach_cell(*global_cage_, [&](Vertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(*global_cage_, cage_vertex_index, vc);

				global_cage_coords_(surface_point_idx, cage_point_idx2) =
					global_cage_coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				sum_lambda += global_cage_coords_(surface_point_idx, cage_point_idx2);

				value<bool>(*global_cage_, cage_vertex_marked, vc) = false;

				return true;
			});

			return true;
		});
	}

	Eigen::VectorXf bind_mvc_handle(Graph& g, std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position){

		Eigen::VectorXf handle_weights; 
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		handle_weights.resize(nbv_cage);
		handle_weights.setZero(); 
		
		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked =
			get_attribute<bool, Vertex>(*global_cage_, "marked_vertices");

		foreach_cell(g, [&](GraphVertex v) -> bool {
			const Vec3& graph_point = value<Vec3>(g, graph_vertex_position, v);
			//uint32 surface_point_idx = value<uint32>(*g, object_vertex_index, v);

			DartMarker dm(*global_cage_);

			float sumMVC = 0.0;
			for (Dart d = global_cage_->begin(), end = global_cage_->end(); d != end; d = global_cage_->next(d))
			{
				Vertex cage_vertex = CMap2::Vertex(d);

				bool vc_marked = value<bool>(*global_cage_, cage_vertex_marked, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{

					const Vec3& cage_point = value<Vec3>(*global_cage_, global_cage_vertex_position_, cage_vertex);
					uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, cage_vertex);

					float mvc_value = cgogn::modeling::compute_mvc(graph_point, d, *global_cage_, cage_point,
																   global_cage_vertex_position_.get());

					handle_weights[cage_point_idx] = mvc_value;

					dm.mark(d);

					value<bool>(*global_cage_, cage_vertex_marked, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			//float sum_lambda = 0.0;

			parallel_foreach_cell(*global_cage_, [&](Vertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(*global_cage_, cage_vertex_index, vc);

				handle_weights[cage_point_idx2] =
					handle_weights[cage_point_idx2] / sumMVC;

				//sum_lambda += global_cage_coords_(surface_point_idx, cage_point_idx2);

				value<bool>(*global_cage_, cage_vertex_marked, vc) = false;

				return true;
			});

			return true;
		});

		return handle_weights; 
	}

	void bind_green(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position)
	{

		const Vec3 NULL_VECTOR(0.0, 0.0, 0);

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked =
			get_attribute<bool, Vertex>(*global_cage_, "marked_vertices");

		std::shared_ptr<Attribute<uint32>> cage_face_index = add_attribute<uint32, Face>(*global_cage_, "face_index");
		uint32 nb_faces_cage = 0;
		foreach_cell(*global_cage_, [&](Face f) -> bool {
			value<uint32>(*global_cage_, cage_face_index, f) = nb_faces_cage++;
			return true;
		});

		std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_normal =
			add_attribute<std::vector<Vec3>, Face>(*global_cage_, "face_normal");

		std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_edge =
			add_attribute<std::vector<Vec3>, Face>(*global_cage_, "face_edge");

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);
		uint32 nbf_cage = nb_cells<Face>(*global_cage_); // Warning valid only for square face (1 face = 2 triangles)

		global_cage_coords_.resize(nbv_object, nbv_cage);
		global_cage_coords_.setZero();

		global_cage_normal_coords_.resize(nbv_object, nbf_cage);
		global_cage_normal_coords_.setZero();

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			foreach_cell(*global_cage_, [&](Face fc) -> bool {
				uint32 cage_face_idx = value<uint32>(*global_cage_, cage_face_index, fc);

				std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*global_cage_, fc);

				// triangle 1

				const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3], face_vertices_[0]};
				const std::vector<Vec3> t1_values = {
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle1[0]),
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle1[1]),
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle1[2])};

				Vec3 t1_normal = (cgogn::geometry::normal(t1_values[0], t1_values[1], t1_values[2])).normalized();

				// triangle 2
				const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2], face_vertices_[3]};
				const std::vector<Vec3> t2_values = {
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle2[0]),
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle2[1]),
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle2[2])};

				Vec3 t2_normal = (cgogn::geometry::normal(t2_values[0], t2_values[1], t2_values[2])).normalized();

				value<std::vector<Vec3>>(*global_cage_, cage_face_normal, fc) = {t1_normal, t2_normal};

				value<std::vector<Vec3>>(*global_cage_, cage_face_edge,
										 fc) = {t1_values[1] - t1_values[0], t1_values[2] - t1_values[1],
												t2_values[1] - t2_values[0], t2_values[2] - t2_values[1]};

				std::vector<Vec3> t1_vj(3);
				std::vector<Vec3> t2_vj(3);
				for (size_t l = 0; l < 3; ++l)
				{
					t1_vj[l] = t1_values[l] - surface_point;
					t2_vj[l] = t2_values[l] - surface_point;
				}

				const Vec3 t1_p_ = (t1_vj[0].dot(t1_normal)) * t1_normal;
				const Vec3 t2_p_ = (t2_vj[0].dot(t2_normal)) * t2_normal;

				Vec3 t1_I = {0.0, 0.0, 0.0};
				std::vector<double> t1_II(3);
				Vec3 t1_s = {0.0, 0.0, 0.0};
				std::vector<Vec3> t1_N(3);

				Vec3 t2_I = {0.0, 0.0, 0.0};
				std::vector<double> t2_II(3);
				Vec3 t2_s = {0.0, 0.0, 0.0};
				std::vector<Vec3> t2_N(3);

				for (size_t k = 0; k < 3; ++k)
				{
					const auto t1_v0 = t1_vj[k];
					const auto t1_v1 = t1_vj[(k + 1) % 3];

					const auto t1_vjpt = ((t1_v0 - t1_p_).cross((t1_v1 - t1_p_))).dot(t1_normal);
					t1_s[k] = t1_vjpt < 0 ? -1.0 : 1.0;

					t1_I[k] = cgogn::modeling::GCTriInt2(t1_p_, t1_v0, t1_v1);
					t1_II[k] = cgogn::modeling::GCTriInt2(NULL_VECTOR, t1_v1, t1_v0);
					t1_N[k] = (t1_v1.cross(t1_v0)).normalized();

					const auto t2_v0 = t2_vj[k];
					const auto t2_v1 = t2_vj[(k + 1) % 3];

					const auto t2_vjpt = ((t2_v0 - t2_p_).cross((t2_v1 - t2_p_))).dot(t2_normal);
					t2_s[k] = t2_vjpt < 0 ? -1.0 : 1.0;

					t2_I[k] = cgogn::modeling::GCTriInt2(t2_p_, t2_v0, t2_v1);
					t2_II[k] = cgogn::modeling::GCTriInt2(NULL_VECTOR, t2_v1, t2_v0);
					t2_N[k] = (t2_v1.cross(t2_v0)).normalized();
				}

				const auto t1_I_ = -abs(t1_s.dot(t1_I));
				const auto t2_I_ = -abs(t2_s.dot(t2_I));

				global_cage_normal_coords_(surface_point_idx, cage_face_idx) = {-t1_I_, -t2_I_};

				Vec3 t1_w = t1_I_ * t1_normal;
				Vec3 t2_w = t2_I_ * t2_normal;
				for (size_t a = 0; a < 3; ++a)
				{
					t1_w += (t1_II[a] * t1_N[a]);
					t2_w += (t2_II[a] * t2_N[a]);
				}

				if (t1_w.norm() > DBL_EPSILON)
				{
					for (size_t l = 0; l < 3; l++)
					{
						const uint32 cage_vertex_idx = value<uint32>(*global_cage_, cage_vertex_index, triangle1[l]);

						const auto Nl1 = t1_N[(l + 1) % 3];
						const auto num = Nl1.dot(t1_w);
						const auto denom = Nl1.dot(t1_vj[l]);

						global_cage_coords_(surface_point_idx, cage_vertex_idx) =
							global_cage_coords_(surface_point_idx, cage_vertex_idx) + num / denom;
					}
				}

				if (t2_w.norm() > DBL_EPSILON)
				{
					for (size_t l = 0; l < 3; l++)
					{
						const uint32 cage_vertex_idx = value<uint32>(*global_cage_, cage_vertex_index, triangle2[l]);

						const auto Nl1 = t2_N[(l + 1) % 3];
						const auto num = Nl1.dot(t2_w);
						const auto denom = Nl1.dot(t2_vj[l]);

						global_cage_coords_(surface_point_idx, cage_vertex_idx) =
							global_cage_coords_(surface_point_idx, cage_vertex_idx) + num / denom;
					}
				}

				return true;
			});

			return true;
		});
	}

	void update_mvc(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, Vertex>(object, "vertex_index");

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


	// TODO
	/*void update_mvc_handle(Graph& g, std::shared_ptr<GraphAttribute<Vec3>>& graph_vertex_position, weights)
	{

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		parallel_foreach_cell(g, [&](Vertex v) -> bool {

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
	}*/

	void update_green(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_face_index =
			cgogn::get_attribute<uint32, Face>(*global_cage_, "face_index");

		std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_normal =
			cgogn::get_attribute<std::vector<Vec3>, Face>(*global_cage_, "face_normal");

		std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_edge =
			cgogn::get_attribute<std::vector<Vec3>, Face>(*global_cage_, "face_edge");

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

			foreach_cell(*global_cage_, [&](Face cf) -> bool {
				uint32 cage_face_idx = value<uint32>(*global_cage_, cage_face_index, cf);

				std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*global_cage_, cf);

				const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3], face_vertices_[0]};

				const std::vector<Vec3> t1_values = {
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle1[0]),
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle1[1]),
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle1[2])};

				const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2], face_vertices_[3]};

				const std::vector<Vec3> t2_values = {
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle2[0]),
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle2[1]),
					value<Vec3>(*global_cage_, global_cage_vertex_position_, triangle2[2])};

				const auto t1_normal = value<std::vector<Vec3>>(*global_cage_, cage_face_normal, cf)[0];
				const auto t2_normal = value<std::vector<Vec3>>(*global_cage_, cage_face_normal, cf)[1];

				// update triangle 1
				const auto t1_u0 = value<std::vector<Vec3>>(*global_cage_, cage_face_edge, cf)[0];
				const auto t1_v0 = value<std::vector<Vec3>>(*global_cage_, cage_face_edge, cf)[1];

				const auto t1_u1 = t1_values[1] - t1_values[0];
				const auto t1_v1 = t1_values[2] - t1_values[1];

				const auto area_face = (t1_u0.cross(t1_v0)).norm() * 0.5;

				double t1_sj =
					sqrt((t1_u1.squaredNorm()) * (t1_v0.squaredNorm()) - 2.0 * (t1_u1.dot(t1_v1)) * (t1_u0.dot(t1_v0)) +
						 (t1_v1.squaredNorm()) * (t1_u0.squaredNorm())) /
					(sqrt8 * area_face);

				new_norm_update_ += global_cage_normal_coords_(vidx, cage_face_idx)[0] * t1_sj * t1_normal;

				// update triangle 2
				const auto t2_u0 = value<std::vector<Vec3>>(*global_cage_, cage_face_edge, cf)[2];
				const auto t2_v0 = value<std::vector<Vec3>>(*global_cage_, cage_face_edge, cf)[3];

				const auto t2_u1 = t2_values[1] - t2_values[0];
				const auto t2_v1 = t2_values[2] - t2_values[1];

				double t2_sj =
					sqrt((t2_u1.squaredNorm()) * (t2_v0.squaredNorm()) - 2.0 * (t2_u1.dot(t2_v1)) * (t2_u0.dot(t2_v0)) +
						 (t2_v1.squaredNorm()) * (t2_u0.squaredNorm())) /
					(sqrt8 * area_face);

				new_norm_update_ += global_cage_normal_coords_(vidx, cage_face_idx)[1] * t2_sj * t2_normal;

				return true;
			});

			value<Vec3>(object, object_vertex_position, v) = new_pos_update_ + new_norm_update_;

			return true;
		});
	}

private:
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
