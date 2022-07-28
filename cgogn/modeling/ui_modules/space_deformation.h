/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_CAGE_DEFORMATION_H_
#define CGOGN_MODULE_CAGE_DEFORMATION_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/geometry/ui_modules/surface_selectionPO.h>
#include <cgogn/modeling/ui_modules/surface_deformation.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <boost/synapse/connect.hpp>

#include <iostream>

#include <string>

namespace cgogn
{

namespace ui
{

using Vec2 = geometry::Vec2;
using Vec3 = geometry::Vec3;
using Mat3 = geometry::Mat3;
using Scalar = geometry::Scalar;

void create_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max)
{
	CMap2::Volume v = add_prism(m, 4);
	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi_1(m, f2))};

	Vec3 center = (bb_min + bb_max) / Scalar(2);
	Vec3 bb_min_ = ((bb_min - center) * 1.5) + center;
	Vec3 bb_max_ = ((bb_max - center) * 1.5) + center;

	value<Vec3>(m, vertex_position, vertices[0]) = bb_min_;
	value<Vec3>(m, vertex_position, vertices[1]) = {bb_min_[0], bb_max_[1], bb_min_[2]};
	value<Vec3>(m, vertex_position, vertices[2]) = {bb_max_[0], bb_max_[1], bb_min_[2]};
	value<Vec3>(m, vertex_position, vertices[3]) = {bb_max_[0], bb_min_[1], bb_min_[2]};

	value<Vec3>(m, vertex_position, vertices[4]) = {bb_min_[0], bb_max_[1], bb_max_[2]};
	value<Vec3>(m, vertex_position, vertices[5]) = {bb_min_[0], bb_min_[1], bb_max_[2]};
	value<Vec3>(m, vertex_position, vertices[6]) = {bb_max_[0], bb_min_[1], bb_max_[2]};
	value<Vec3>(m, vertex_position, vertices[7]) = {bb_max_[0], bb_max_[1], bb_max_[2]};
}

// Github SuperBoubek QMVC https://github.com/superboubek/QMVC/blob/master/coordinates/mvc/mvc.h
double getAngleBetweenUnitVectors(const Vec3& a, const Vec3& b)
{
	return 2.0 * asin((a - b).norm() / 2.0);
}

float compute_mvc(const Vec3& surface_point, Dart vertex, CMap2& cage, const Vec3& cage_point,
				  CMap2::Attribute<Vec3>* cage_position)
{

	double r = (cage_point - surface_point).norm();

	if (r == 0)
	{
		std::cerr << "Alert" << std::endl;
	}

	double sumU(0.);

	Dart it = vertex;

	do
	{
		Vec3 vi = value<Vec3>(cage, cage_position, CMap2::Vertex(it));
		Vec3 vj = value<Vec3>(cage, cage_position, CMap2::Vertex(phi1(cage, it)));
		Vec3 vk = value<Vec3>(cage, cage_position, CMap2::Vertex(phi_1(cage, it)));

		// double Bjk = cgogn::geometry::angle((vj - surface_point), (vk - surface_point));
		// double Bij = cgogn::geometry::angle((vi - surface_point), (vj - surface_point));
		// double Bki = cgogn::geometry::angle((vk - surface_point), (vi - surface_point));

		Vec3 ei = (vi - surface_point).normalized();
		Vec3 ej = (vj - surface_point).normalized();
		Vec3 ek = (vk - surface_point).normalized();

		double Bjk = getAngleBetweenUnitVectors(ej, ek);
		double Bij = getAngleBetweenUnitVectors(ei, ej);
		double Bki = getAngleBetweenUnitVectors(ek, ei);

		Vec3 eiej = ei.cross(ej);
		Vec3 ejek = ej.cross(ek);
		Vec3 ekei = ek.cross(ei);

		Vec3 nij = eiej.normalized();
		Vec3 njk = ejek.normalized();
		Vec3 nki = ekei.normalized();

		double ui = (Bjk + (Bij * (nij.dot(njk))) + (Bki * (nki.dot(njk)))) / (2.f * ei.dot(njk));

		sumU += ui;

		it = phi<2, 1>(cage, it);
	} while (it != vertex);

	return (1.0f / r) * sumU;
}

const double GCTriInt(const Vec3& p, const Vec3& v1, const Vec3& v2, const Vec3& nu)
{
	const Vec3 v2_v1 = v2 - v1;
	const Vec3 p_v1 = p - v1;
	const Vec3 v1_p = v1 - p;
	const Vec3 v2_p = v2 - p;

	const Vec3 p_nu = p - nu;

	const double alpha = std::acos((v2_v1.dot(p_v1)) / (v2_v1.norm() * p_v1.norm()));
	const double beta = std::acos((v1_p.dot(v2_p)) / (v1_p.norm() * v2_p.norm()));

	const double lambda = p_v1.squaredNorm() * std::sin(alpha) * std::sin(alpha);

	const double c = p_nu.squaredNorm();

	const double sqrt_c = sqrt(c);
	const double sqrt_lambda = sqrt(lambda);

	const std::array<double, 2> theta = {M_PI - alpha, M_PI - alpha - beta};
	std::array<double, 2> I;

	for (size_t i = 0; i < 2; ++i)
	{
		const double S = std::sin(theta[i]);
		const double C = std::cos(theta[i]);

		const double SS = S * S;

		const auto sign = S < 0 ? -1.0 : 1.0;

		const auto tan_part = 2 * sqrt_c * std::atan2(sqrt_c * C, sqrt(lambda + SS * c));

		const auto log_part =
			log((2 * sqrt_lambda * SS / (std::pow(1.0 - C, 2)) *
				 (1.0 - ((2 * c * C) / ((c * (1 + C) + lambda + sqrt((lambda * lambda) + (lambda * c * SS))))))));

		I[i] = (-sign * 0.5) * (tan_part + (sqrt_lambda * log_part));
	}

	return (-1.0 / (4.0 * M_PI)) * abs(I[0] - I[1] - sqrt_c * beta);
}

// https://github.com/blaisebundle/green_cage_deformer/blob/master/src/green_cage_deformer.cc
namespace
{
const Vec3 NULL_VECTOR(0.0, 0.0, 0);
constexpr auto TOLERANCE = 1e-6;
} // namespace

namespace MathConstants
{
constexpr auto PI = 3.14159265358979323846;
constexpr auto ONE_OVER_FOUR_PI = 0.07957747154594767;
constexpr auto SQRT8 = 2.828427124746190097603;
}; // namespace MathConstants

const double GCTriInt2(const Vec3& p, const Vec3& v1, const Vec3& v2)
{
	const auto v2_v1 = v2 - v1;
	const auto p_v1 = p - v1;
	const auto v1_p = v1 - p;
	const auto v2_p = v2 - p;

	const auto alpha =
		std::acos(std::min(std::max(((v2_v1).dot((p_v1))) / ((v2_v1).norm() * (p_v1).norm()), -1.0), 1.0));
	if (abs(alpha - MathConstants::PI) < TOLERANCE || abs(alpha) < TOLERANCE)
	{
		return 0.0;
	}

	const auto beta =
		std::acos(std::min(std::max((((v1_p).dot((v2_p)))) / ((v1_p).norm() * (v2_p).norm()), -1.0), 1.0));
	const auto lambda = p_v1.squaredNorm() * std::sin(alpha) * std::sin(alpha);
	const auto c = p.squaredNorm();

	const auto sqrt_c = sqrt(c);
	const auto sqrt_lambda = sqrt(lambda);

	const std::array<double, 2> theta{MathConstants::PI - alpha, MathConstants::PI - alpha - beta};
	std::array<double, 2> I;

	for (size_t i = 0; i < 2; ++i)
	{
		const auto S = std::sin(theta[i]);
		const auto C = std::cos(theta[i]);
		const auto sign = S < 0 ? -1.0 : 1.0;

		const auto SS = S * S;
		const auto half_sign = (-sign * 0.5);
		const auto tan_part = (2 * sqrt_c * std::atan2((sqrt_c * C), sqrt(lambda + (SS * c))));
		const auto log_part =
			log(((2 * sqrt_lambda * SS) / std::pow(1.0 - C, 2)) *
				(1.0 - ((2 * c * C) / ((c * (1 + C) + lambda + sqrt((lambda * lambda) + (lambda * c * SS)))))));

		I[i] = half_sign * (tan_part + (sqrt_lambda * log_part));
	}

	return -MathConstants::ONE_OVER_FOUR_PI * abs(I[0] - I[1] - sqrt_c * beta);
}

template <typename MESH>

class SpaceDeformation : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "SpaceDeformation can only be used with meshes of dimension 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	struct Cage_data
	{

		Cage_data() : cage_vertex_position_(nullptr), influence_set_(nullptr), control_set_(nullptr), local_def(false)
		{
		}

		~Cage_data()
		{
		}

		// CGOGN_NOT_COPYABLE_NOR_MOVABLE(Weights);

		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> coords_;
		Eigen::Matrix<Vec2, Eigen::Dynamic, Eigen::Dynamic> n_coords_;

		std::shared_ptr<Attribute<Vec3>> cage_vertex_position_;

		CellsSet<MESH, Vertex>* influence_set_;
		CellsSet<MESH, Vertex>* control_set_;

		MESH* influence_cage_;

		bool local_def;

		std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;
	};

	struct Parameters
	{
		Parameters() : vertex_position_(nullptr), list_cage_(100, nullptr), new_cage_(false), nb_cage(0)
		{
		}

		~Parameters()
		{
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		std::vector<MESH*> list_cage_;

		std::unique_ptr<CellCache<MESH>> working_cells_;

		// CellsSet<MESH, Vertex>* last_influence_set_;

		bool new_cage_;

		std::shared_ptr<boost::synapse::connection> cells_set_connection_;

		int nb_cage;
	};

public:
	SpaceDeformation(const App& app)
		: Module(app, "SpaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr)
	{
	}
	~SpaceDeformation()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		p.vertex_position_ = vertex_position;
	}

	MESH* generate_global_cage(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position, int i)
	{
		// const std::string& m_name = mesh_provider_->mesh_name(m);
		Parameters& p = parameters_[&m];

		MESH* cage = mesh_provider_->add_mesh("cage" + std::to_string(i)); // (m_name + "_cage");
		std::shared_ptr<Attribute<Vec3>> cage_vertex_position = cgogn::add_attribute<Vec3, Vertex>(*cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*cage, cage_vertex_position);

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		create_box(*cage, cage_vertex_position.get(), md.bb_min_, md.bb_max_);

		mesh_provider_->emit_connectivity_changed(*cage);
		mesh_provider_->emit_attribute_changed(*cage, cage_vertex_position.get());

		Cage_data& cd = cage_data_[cage];

		p.list_cage_[i] = cage;
		// p.list_weights_.push_back(Weights());
		cd.cage_vertex_position_ = cage_vertex_position;

		ui::View* v1 = app_.current_view();

		surface_render_->set_vertex_position(*v1, *cage, cage_vertex_position);
		surface_render_->set_render_faces(*v1, *cage, false);

		return cage;
	}

	MESH* generate_local_cage(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position,
							  CellsSet<MESH, Vertex>* control_set, int i)
	{

		MESH* l_cage = mesh_provider_->add_mesh("cage" + std::to_string(i));

		std::shared_ptr<Attribute<Vec3>> l_cage_vertex_position =
			cgogn::add_attribute<Vec3, Vertex>(*l_cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*l_cage, l_cage_vertex_position);

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 l_min = {bb_max[0], bb_max[1], bb_max[2]};
		Vec3 l_max = {bb_min[0], bb_min[1], bb_min[2]};

		Parameters& p = parameters_[&m];

		control_set->foreach_cell([&](Vertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);

			for (size_t j = 0; j < 3; j++)
			{
				if (pos[j] < l_min[j])
				{
					l_min[j] = pos[j];
				}

				if (pos[j] > l_max[j])
				{
					l_max[j] = pos[j];
				}
			}
		});

		create_box(*l_cage, l_cage_vertex_position.get(), l_min, l_max);

		mesh_provider_->emit_connectivity_changed(*l_cage);
		mesh_provider_->emit_attribute_changed(*l_cage, l_cage_vertex_position.get());

		Cage_data& cd = cage_data_[l_cage];

		p.list_cage_[i] = l_cage;

		cd.cage_vertex_position_ = l_cage_vertex_position;
		cd.local_def = true;
		cd.control_set_ = control_set;

		View* v1 = app_.current_view();

		surface_render_->set_vertex_position(*v1, *l_cage, l_cage_vertex_position);
		surface_render_->set_render_faces(*v1, *l_cage, false);

		std::shared_ptr<Attribute<Vec3>> l_cage_vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*l_cage, "normal");

		surface_diff_pptes_->compute_normal(*l_cage, l_cage_vertex_position.get(), l_cage_vertex_normal.get());

		MESH* i_cage = mesh_provider_->add_mesh("i_cage" + std::to_string(i));
		std::shared_ptr<Attribute<Vec3>> i_cage_vertex_position =
			cgogn::add_attribute<Vec3, Vertex>(*i_cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*i_cage, i_cage_vertex_position);

		Vec3 i_center = (l_min + l_max) / Scalar(2);
		Vec3 bb_min_ = ((l_min - i_center) * 1.5) + i_center;
		Vec3 bb_max_ = ((l_max - i_center) * 1.5) + i_center;

		create_box(*i_cage, i_cage_vertex_position.get(), bb_min_, bb_max_);

		mesh_provider_->emit_connectivity_changed(*i_cage);
		mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

		cd.influence_cage_ = i_cage;

		// surface_render_->set_vertex_position(*v1, *i_cage, i_cage_vertex_position);
		// surface_render_->set_render_faces(*v1, *i_cage, false);

		MeshData<Mesh>& i_md = mesh_provider_->mesh_data(*i_cage);
		CellsSet<Mesh, Vertex>& i_set = i_md.add_cells_set<Vertex>();

		std::shared_ptr<Attribute<uint32>> i_cage_vertex_index =
			cgogn::add_attribute<uint32, Vertex>(*i_cage, "weight_index");
		uint32 nb_vertices_cage = 0;
		foreach_cell(*i_cage, [&](Vertex v) -> bool {
			value<uint32>(*i_cage, i_cage_vertex_index, v) = nb_vertices_cage++;
			return true;
		});

		std::shared_ptr<Attribute<bool>> i_cage_vertex_marked = cgogn::add_attribute<bool, Vertex>(*i_cage, "marked_vertex");
		parallel_foreach_cell(*i_cage, [&](Vertex v) -> bool {
			value<bool>(*i_cage, i_cage_vertex_marked, v) = false;
			return true;
		});

		std::shared_ptr<Attribute<Vec3>> m_vertex_position = p.vertex_position_; 

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(m, m_vertex_position, v);
			
			bool inside_cage = local_mvc_pt_surface(surface_point, *i_cage, i_cage_vertex_position ); 
			if (inside_cage)
			{
				i_set.select(v); 
			}

			return true; 
		});  

		//cd.influence_set_ = &i_set; 

		return l_cage;
	}

	void bind_object_mvc(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position, MESH& cage,
						 const std::shared_ptr<Attribute<Vec3>>& cage_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "weight_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index = add_attribute<uint32, Vertex>(cage, "weight_index");
		uint32 nb_vertices_cage = 0;
		foreach_cell(cage, [&](Vertex v) -> bool {
			value<uint32>(cage, cage_vertex_index, v) = nb_vertices_cage++;
			return true;
		});

		std::shared_ptr<Attribute<bool>> cage_vertex_marked = add_attribute<bool, Vertex>(cage, "marked_vertex");
		parallel_foreach_cell(cage, [&](Vertex v) -> bool {
			value<bool>(cage, cage_vertex_marked, v) = false;
			return true;
		});

		Parameters& p = parameters_[&object];
		Cage_data& cd = cage_data_[&cage];

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(cage);

		cd.coords_.resize(nbv_object, nbv_cage);

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			DartMarker dm(cage);
			float sumMVC = 0.0;
			for (Dart d = cage.begin(), end = cage.end(); d != end; d = cage.next(d))
			{
				Vertex cage_vertex = CMap2::Vertex(d);
				bool vc_marked = value<bool>(cage, cage_vertex_marked, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{

					const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cage_vertex);
					uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cage_vertex);

					float mvc_value = compute_mvc(surface_point, d, cage, cage_point, cage_vertex_position.get());

					cd.coords_(surface_point_idx, cage_point_idx) = mvc_value;
					dm.mark(d);

					value<bool>(cage, cage_vertex_marked, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			float sum_lambda = 0.0;

			parallel_foreach_cell(cage, [&](Vertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(cage, cage_vertex_index, vc);

				cd.coords_(surface_point_idx, cage_point_idx2) =
					cd.coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				// sum_lambda += p.list_weights_[i].coords_(surface_point_idx, cage_point_idx2);

				value<bool>(cage, cage_vertex_marked, vc) = false;

				return true;
			});

			// std::cout << "sum_lambda " << sum_lambda << std::endl;

			return true;
		});

		// std::cout << p.coords_ << std::endl;

		cd.cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				&cage, [&](Attribute<Vec3>* attribute) {
					if (cd.cage_vertex_position_.get() == attribute)
					{

						std::shared_ptr<Attribute<uint32>> object_vertex_index =
							cgogn::get_attribute<uint32, Vertex>(object, "weight_index");
						std::shared_ptr<Attribute<uint32>> cage_vertex_index =
							cgogn::get_attribute<uint32, Vertex>(cage, "weight_index");

						parallel_foreach_cell(object, [&](Vertex v) -> bool {
							uint32 vidx = value<uint32>(object, object_vertex_index, v);

							Vec3 new_pos_ = {0.0, 0.0, 0.0};

							foreach_cell(cage, [&](Vertex cv) -> bool {
								const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cv);
								uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cv);

								new_pos_ += cd.coords_(vidx, cage_point_idx) * cage_point;

								return true;
							});

							value<Vec3>(object, object_vertex_position, v) = new_pos_;
							return true;
						});

						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					}
				});
	}

	void bind_local_mvc(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position, MESH& cage,
						const std::shared_ptr<Attribute<Vec3>>& cage_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "weight_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index = add_attribute<uint32, Vertex>(cage, "weight_index");
		uint32 nb_vertices_cage = 0;
		foreach_cell(cage, [&](Vertex v) -> bool {
			value<uint32>(cage, cage_vertex_index, v) = nb_vertices_cage++;
			return true;
		});

		std::shared_ptr<Attribute<bool>> cage_vertex_marked = add_attribute<bool, Vertex>(cage, "marked_vertex");
		parallel_foreach_cell(cage, [&](Vertex v) -> bool {
			value<bool>(cage, cage_vertex_marked, v) = false;
			return true;
		});

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(cage);

		Parameters& p = parameters_[&object];
		Cage_data& cd = cage_data_[&cage];

		cd.coords_.resize(nbv_object, nbv_cage);

		cd.influence_set_->foreach_cell([&](Vertex v) {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			DartMarker dm(cage);
			float sumMVC = 0.0;
			for (Dart d = cage.begin(), end = cage.end(); d != end; d = cage.next(d))
			{
				Vertex cage_vertex = CMap2::Vertex(d);
				bool vc_marked = value<bool>(cage, cage_vertex_marked, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{

					const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cage_vertex);
					uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cage_vertex);

					float mvc_value = compute_mvc(surface_point, d, cage, cage_point, cage_vertex_position.get());

					cd.coords_(surface_point_idx, cage_point_idx) = mvc_value;
					dm.mark(d);

					value<bool>(cage, cage_vertex_marked, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			float sum_lambda = 0.0;

			parallel_foreach_cell(cage, [&](Vertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(cage, cage_vertex_index, vc);

				cd.coords_(surface_point_idx, cage_point_idx2) =
					cd.coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				// sum_lambda += p.list_weights_[i].coords_(surface_point_idx, cage_point_idx2);

				value<bool>(cage, cage_vertex_marked, vc) = false;

				return true;
			});
		});

		cd.cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				&cage, [&](Attribute<Vec3>* attribute) {
					if (cd.cage_vertex_position_.get() == attribute)
					{

						std::shared_ptr<Attribute<uint32>> object_vertex_index =
							cgogn::get_attribute<uint32, Vertex>(object, "weight_index");

						std::shared_ptr<Attribute<uint32>> cage_vertex_index =
							cgogn::get_attribute<uint32, Vertex>(cage, "weight_index");

						cd.influence_set_->foreach_cell([&](Vertex v) -> bool {
							uint32 vidx = value<uint32>(object, object_vertex_index, v);

							Vec3 new_pos_ = {0.0, 0.0, 0.0};

							foreach_cell(cage, [&](Vertex cv) -> bool {
								const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cv);
								uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cv);

								new_pos_ += cd.coords_(vidx, cage_point_idx) * cage_point;

								return true;
							});

							value<Vec3>(object, object_vertex_position, v) = new_pos_;
							return true;
						});

						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					}
				});
	}

	bool local_mvc_pt_surface(Vec3 pt, MESH& cage, const std::shared_ptr<Attribute<Vec3>>& cage_vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> i_cage_vertex_index = get_attribute<uint32, Vertex>(cage, "weight_index");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked = get_attribute<bool, Vertex>(cage, "marked_vertex");

		DartMarker dm(cage);
		float sumMVC = 0.0;

		bool checked = true; 
		for (Dart d = cage.begin(), end = cage.end(); d != end; d = cage.next(d))
		{
			Vertex cage_vertex = CMap2::Vertex(d);
			bool vc_marked = value<bool>(cage, cage_vertex_marked, cage_vertex);

			if (!dm.is_marked(d) && !vc_marked)
			{

				const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cage_vertex);

				float mvc_value = compute_mvc(pt, d, cage, cage_point, cage_vertex_position.get());

				if (mvc_value < 0)
				{
					checked = false;
					break;
				}

				dm.mark(d);

				value<bool>(cage, cage_vertex_marked, cage_vertex) = true;
			}
		}

		parallel_foreach_cell(cage, [&](Vertex vc) -> bool {
				value<bool>(cage, cage_vertex_marked, vc) = false;
				return true;
			}); 

			return checked;
		}

		void bind_object_green(MESH & object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position,
							   MESH& cage, const std::shared_ptr<Attribute<Vec3>>& cage_vertex_position)
		{

			std::shared_ptr<Attribute<uint32>> object_vertex_index =
				get_attribute<uint32, Vertex>(object, "weight_index");

			std::shared_ptr<Attribute<uint32>> cage_vertex_index = add_attribute<uint32, Vertex>(cage, "weight_index");
			uint32 nb_vertices_cage = 0;
			foreach_cell(cage, [&](Vertex v) -> bool {
				value<uint32>(cage, cage_vertex_index, v) = nb_vertices_cage++;
				return true;
			});

			std::shared_ptr<Attribute<uint32>> cage_face_index = add_attribute<uint32, Face>(cage, "face_index");
			uint32 nb_faces_cage = 0;
			foreach_cell(cage, [&](Face f) -> bool {
				value<uint32>(cage, cage_face_index, f) = nb_faces_cage++;
				return true;
			});

			std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_normal =
				add_attribute<std::vector<Vec3>, Face>(cage, "face_normal");

			std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_edge =
				add_attribute<std::vector<Vec3>, Face>(cage, "face_edge");

			Parameters& p = parameters_[&object];
			Cage_data& cd = cage_data_[&cage];

			uint32 nbv_object = nb_cells<Vertex>(object);
			uint32 nbv_cage = nb_cells<Vertex>(cage);

			uint32 nbf_cage = nb_cells<Face>(cage); // Warning valid only for square face (1 face = 2 triangles)

			cd.coords_.resize(nbv_object, nbv_cage);
			cd.coords_.setZero();

			cd.n_coords_.resize(nbv_object, nbf_cage);
			cd.n_coords_.setZero();

			parallel_foreach_cell(object, [&](Vertex v) -> bool {
				const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
				uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

				foreach_cell(cage, [&](Face fc) -> bool {
					uint32 cage_face_idx = value<uint32>(cage, cage_face_index, fc);

					// Dart d1 = fc.dart;

					std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(cage, fc);

					// triangle 1
					// const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[0], face_vertices_[3],
					// face_vertices_[1]};
					const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3],
																  face_vertices_[0]};
					const std::vector<Vec3> t1_values = {value<Vec3>(cage, cage_vertex_position, triangle1[0]),
														 value<Vec3>(cage, cage_vertex_position, triangle1[1]),
														 value<Vec3>(cage, cage_vertex_position, triangle1[2])};

					Vec3 t1_normal = (cgogn::geometry::normal(t1_values[0], t1_values[1], t1_values[2])).normalized();

					const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2],
																  face_vertices_[3]};
					const std::vector<Vec3> t2_values = {value<Vec3>(cage, cage_vertex_position, triangle2[0]),
														 value<Vec3>(cage, cage_vertex_position, triangle2[1]),
														 value<Vec3>(cage, cage_vertex_position, triangle2[2])};

					Vec3 t2_normal = (cgogn::geometry::normal(t2_values[0], t2_values[1], t2_values[2])).normalized();

					value<std::vector<Vec3>>(cage, cage_face_normal, fc) = {t1_normal, t2_normal};

					value<std::vector<Vec3>>(cage, cage_face_edge,
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
						// I[l] = GCTriInt(p_, v0, v1, {0.0f, 0.0f, 0.0f});
						// II[l] = GCTriInt({0.0f, 0.0f, 0.0f}, v1, v0, {0.0f, 0.0f, 0.0f});
						t1_I[k] = GCTriInt2(t1_p_, t1_v0, t1_v1);
						t1_II[k] = GCTriInt2(NULL_VECTOR, t1_v1, t1_v0);
						t1_N[k] = (t1_v1.cross(t1_v0)).normalized();

						const auto t2_v0 = t2_vj[k];
						const auto t2_v1 = t2_vj[(k + 1) % 3];

						const auto t2_vjpt = ((t2_v0 - t2_p_).cross((t2_v1 - t2_p_))).dot(t2_normal);
						t2_s[k] = t2_vjpt < 0 ? -1.0 : 1.0;
						// I[l] = GCTriInt(p_, v0, v1, {0.0f, 0.0f, 0.0f});
						// II[l] = GCTriInt({0.0f, 0.0f, 0.0f}, v1, v0, {0.0f, 0.0f, 0.0f});
						t2_I[k] = GCTriInt2(t2_p_, t2_v0, t2_v1);
						t2_II[k] = GCTriInt2(NULL_VECTOR, t2_v1, t2_v0);
						t2_N[k] = (t2_v1.cross(t2_v0)).normalized();
					}

					const auto t1_I_ = -abs(t1_s.dot(t1_I));
					const auto t2_I_ = -abs(t2_s.dot(t2_I));

					cd.n_coords_(surface_point_idx, cage_face_idx) = {-t1_I_, -t2_I_};

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
							const uint32 cage_vertex_idx = value<uint32>(cage, cage_vertex_index, triangle1[l]);

							const auto Nl1 = t1_N[(l + 1) % 3];
							const auto num = Nl1.dot(t1_w);
							const auto denom = Nl1.dot(t1_vj[l]);

							cd.coords_(surface_point_idx, cage_vertex_idx) =
								cd.coords_(surface_point_idx, cage_vertex_idx) + num / denom;
						}
					}

					if (t2_w.norm() > DBL_EPSILON)
					{
						for (size_t l = 0; l < 3; l++)
						{
							const uint32 cage_vertex_idx = value<uint32>(cage, cage_vertex_index, triangle2[l]);

							const auto Nl1 = t2_N[(l + 1) % 3];
							const auto num = Nl1.dot(t2_w);
							const auto denom = Nl1.dot(t2_vj[l]);

							cd.coords_(surface_point_idx, cage_vertex_idx) =
								cd.coords_(surface_point_idx, cage_vertex_idx) + num / denom;
						}
					}

					return true;
				});

				return true;
			});

			cd.cage_attribute_update_connection_ =
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					&cage, [&](Attribute<Vec3>* attribute) {
						if (cd.cage_vertex_position_.get() == attribute)
						{

							std::shared_ptr<Attribute<uint32>> object_vertex_index =
								cgogn::get_attribute<uint32, Vertex>(object, "weight_index");

							std::shared_ptr<Attribute<uint32>> cage_vertex_index =
								cgogn::get_attribute<uint32, Vertex>(cage, "weight_index");

							std::shared_ptr<Attribute<uint32>> cage_face_index =
								cgogn::get_attribute<uint32, Face>(cage, "face_index");

							std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_normal =
								cgogn::get_attribute<std::vector<Vec3>, Face>(cage, "face_normal");

							std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_edge =
								cgogn::get_attribute<std::vector<Vec3>, Face>(cage, "face_edge");

							parallel_foreach_cell(object, [&](Vertex v) -> bool {
								uint32 vidx = value<uint32>(object, object_vertex_index, v);

								Vec3 new_pos_update_ = {0.0, 0.0, 0.0};

								const auto sqrt8 = sqrt(8);

								foreach_cell(cage, [&](Vertex cv) -> bool {
									const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cv);
									uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cv);

									new_pos_update_ += cd.coords_(vidx, cage_point_idx) * cage_point;

									return true;
								});

								Vec3 new_norm_update_ = {0.0, 0.0, 0.0};
								foreach_cell(cage, [&](Face cf) -> bool {
									uint32 cage_face_idx = value<uint32>(cage, cage_face_index, cf);

									std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(cage, cf);

									const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3],
																				  face_vertices_[0]};
									const std::vector<Vec3> t1_values = {
										value<Vec3>(cage, cage_vertex_position, triangle1[0]),
										value<Vec3>(cage, cage_vertex_position, triangle1[1]),
										value<Vec3>(cage, cage_vertex_position, triangle1[2])};

									const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2],
																				  face_vertices_[3]};
									const std::vector<Vec3> t2_values = {
										value<Vec3>(cage, cage_vertex_position, triangle2[0]),
										value<Vec3>(cage, cage_vertex_position, triangle2[1]),
										value<Vec3>(cage, cage_vertex_position, triangle2[2])};

									const auto t1_normal = value<std::vector<Vec3>>(cage, cage_face_normal, cf)[0];
									const auto t2_normal = value<std::vector<Vec3>>(cage, cage_face_normal, cf)[1];

									// update triangle 1
									const auto t1_u0 = value<std::vector<Vec3>>(cage, cage_face_edge, cf)[0];
									const auto t1_v0 = value<std::vector<Vec3>>(cage, cage_face_edge, cf)[1];

									const auto t1_u1 = t1_values[1] - t1_values[0];
									const auto t1_v1 = t1_values[2] - t1_values[1];

									const auto area_face = (t1_u0.cross(t1_v0)).norm() * 0.5;

									double t1_sj = sqrt((t1_u1.squaredNorm()) * (t1_v0.squaredNorm()) -
														2.0 * (t1_u1.dot(t1_v1)) * (t1_u0.dot(t1_v0)) +
														(t1_v1.squaredNorm()) * (t1_u0.squaredNorm())) /
												   (sqrt8 * area_face);

									new_norm_update_ += cd.n_coords_(vidx, cage_face_idx)[0] * t1_sj * t1_normal;

									// update triangle 2
									const auto t2_u0 = value<std::vector<Vec3>>(cage, cage_face_edge, cf)[2];
									const auto t2_v0 = value<std::vector<Vec3>>(cage, cage_face_edge, cf)[3];

									const auto t2_u1 = t2_values[1] - t2_values[0];
									const auto t2_v1 = t2_values[2] - t2_values[1];

									double t2_sj = sqrt((t2_u1.squaredNorm()) * (t2_v0.squaredNorm()) -
														2.0 * (t2_u1.dot(t2_v1)) * (t2_u0.dot(t2_v0)) +
														(t2_v1.squaredNorm()) * (t2_u0.squaredNorm())) /
												   (sqrt8 * area_face);

									new_norm_update_ += cd.n_coords_(vidx, cage_face_idx)[1] * t2_sj * t2_normal;

									return true;
								});

								value<Vec3>(object, object_vertex_position, v) = new_pos_update_ + new_norm_update_;
								return true;
							});

							mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
						}
					});
		}

		void bind_local_green(MESH & object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position, MESH& cage,
							  const std::shared_ptr<Attribute<Vec3>>& cage_vertex_position)
		{

			std::shared_ptr<Attribute<uint32>> object_vertex_index =
				get_attribute<uint32, Vertex>(object, "weight_index");

			std::shared_ptr<Attribute<uint32>> cage_vertex_index = add_attribute<uint32, Vertex>(cage, "weight_index");
			uint32 nb_vertices_cage = 0;
			foreach_cell(cage, [&](Vertex v) -> bool {
				value<uint32>(cage, cage_vertex_index, v) = nb_vertices_cage++;
				return true;
			});

			std::shared_ptr<Attribute<uint32>> cage_face_index = add_attribute<uint32, Face>(cage, "face_index");
			uint32 nb_faces_cage = 0;
			foreach_cell(cage, [&](Face f) -> bool {
				value<uint32>(cage, cage_face_index, f) = nb_faces_cage++;
				return true;
			});

			std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_normal =
				add_attribute<std::vector<Vec3>, Face>(cage, "face_normal");

			std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_edge =
				add_attribute<std::vector<Vec3>, Face>(cage, "face_edge");

			Parameters& p = parameters_[&object];
			Cage_data& cd = cage_data_[&cage];

			uint32 nbv_object = nb_cells<Vertex>(object);
			uint32 nbv_cage = nb_cells<Vertex>(cage);

			uint32 nbf_cage = nb_cells<Face>(cage); // Warning valid only for square face (1 face = 2 triangles)

			cd.coords_.resize(nbv_object, nbv_cage);
			cd.coords_.setZero();

			cd.n_coords_.resize(nbv_object, nbf_cage);
			cd.n_coords_.setZero();

			cd.influence_set_->foreach_cell([&](Vertex v) {
				const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
				uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

				foreach_cell(cage, [&](Face fc) -> bool {
					uint32 cage_face_idx = value<uint32>(cage, cage_face_index, fc);

					// Dart d1 = fc.dart;

					std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(cage, fc);

					const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3],
																  face_vertices_[0]};
					const std::vector<Vec3> t1_values = {value<Vec3>(cage, cage_vertex_position, triangle1[0]),
														 value<Vec3>(cage, cage_vertex_position, triangle1[1]),
														 value<Vec3>(cage, cage_vertex_position, triangle1[2])};

					Vec3 t1_normal = (cgogn::geometry::normal(t1_values[0], t1_values[1], t1_values[2])).normalized();

					const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2],
																  face_vertices_[3]};
					const std::vector<Vec3> t2_values = {value<Vec3>(cage, cage_vertex_position, triangle2[0]),
														 value<Vec3>(cage, cage_vertex_position, triangle2[1]),
														 value<Vec3>(cage, cage_vertex_position, triangle2[2])};

					Vec3 t2_normal = (cgogn::geometry::normal(t2_values[0], t2_values[1], t2_values[2])).normalized();

					value<std::vector<Vec3>>(cage, cage_face_normal, fc) = {t1_normal, t2_normal};

					value<std::vector<Vec3>>(cage, cage_face_edge,
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
						// I[l] = GCTriInt(p_, v0, v1, {0.0f, 0.0f, 0.0f});
						// II[l] = GCTriInt({0.0f, 0.0f, 0.0f}, v1, v0, {0.0f, 0.0f, 0.0f});
						t1_I[k] = GCTriInt2(t1_p_, t1_v0, t1_v1);
						t1_II[k] = GCTriInt2(NULL_VECTOR, t1_v1, t1_v0);
						t1_N[k] = (t1_v1.cross(t1_v0)).normalized();

						const auto t2_v0 = t2_vj[k];
						const auto t2_v1 = t2_vj[(k + 1) % 3];

						const auto t2_vjpt = ((t2_v0 - t2_p_).cross((t2_v1 - t2_p_))).dot(t2_normal);
						t2_s[k] = t2_vjpt < 0 ? -1.0 : 1.0;
						// I[l] = GCTriInt(p_, v0, v1, {0.0f, 0.0f, 0.0f});
						// II[l] = GCTriInt({0.0f, 0.0f, 0.0f}, v1, v0, {0.0f, 0.0f, 0.0f});
						t2_I[k] = GCTriInt2(t2_p_, t2_v0, t2_v1);
						t2_II[k] = GCTriInt2(NULL_VECTOR, t2_v1, t2_v0);
						t2_N[k] = (t2_v1.cross(t2_v0)).normalized();
					}

					const auto t1_I_ = -abs(t1_s.dot(t1_I));
					const auto t2_I_ = -abs(t2_s.dot(t2_I));

					cd.n_coords_(surface_point_idx, cage_face_idx) = {-t1_I_, -t2_I_};

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
							const uint32 cage_vertex_idx = value<uint32>(cage, cage_vertex_index, triangle1[l]);

							const auto Nl1 = t1_N[(l + 1) % 3];
							const auto num = Nl1.dot(t1_w);
							const auto denom = Nl1.dot(t1_vj[l]);

							cd.coords_(surface_point_idx, cage_vertex_idx) =
								cd.coords_(surface_point_idx, cage_vertex_idx) + num / denom;
						}
					}

					if (t2_w.norm() > DBL_EPSILON)
					{
						for (size_t l = 0; l < 3; l++)
						{
							const uint32 cage_vertex_idx = value<uint32>(cage, cage_vertex_index, triangle2[l]);

							const auto Nl1 = t2_N[(l + 1) % 3];
							const auto num = Nl1.dot(t2_w);
							const auto denom = Nl1.dot(t2_vj[l]);

							cd.coords_(surface_point_idx, cage_vertex_idx) =
								cd.coords_(surface_point_idx, cage_vertex_idx) + num / denom;
						}
					}

					return true;
				});
			});

			cd.cage_attribute_update_connection_ =
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					&cage, [&](Attribute<Vec3>* attribute) {
						if (cd.cage_vertex_position_.get() == attribute)
						{

							std::shared_ptr<Attribute<uint32>> object_vertex_index =
								cgogn::get_attribute<uint32, Vertex>(object, "weight_index");

							std::shared_ptr<Attribute<uint32>> cage_vertex_index =
								cgogn::get_attribute<uint32, Vertex>(cage, "weight_index");

							std::shared_ptr<Attribute<uint32>> cage_face_index =
								cgogn::get_attribute<uint32, Face>(cage, "face_index");

							std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_normal =
								cgogn::get_attribute<std::vector<Vec3>, Face>(cage, "face_normal");

							std::shared_ptr<Attribute<std::vector<Vec3>>> cage_face_edge =
								cgogn::get_attribute<std::vector<Vec3>, Face>(cage, "face_edge");

							cd.influence_set_->foreach_cell([&](Vertex v) {
								uint32 vidx = value<uint32>(object, object_vertex_index, v);

								Vec3 new_pos_update_ = {0.0, 0.0, 0.0};

								const auto sqrt8 = sqrt(8);

								foreach_cell(cage, [&](Vertex cv) -> bool {
									const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cv);
									uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cv);

									new_pos_update_ += cd.coords_(vidx, cage_point_idx) * cage_point;

									return true;
								});

								Vec3 new_norm_update_ = {0.0, 0.0, 0.0};
								foreach_cell(cage, [&](Face cf) -> bool {
									uint32 cage_face_idx = value<uint32>(cage, cage_face_index, cf);

									std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(cage, cf);

									const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3],
																				  face_vertices_[0]};
									const std::vector<Vec3> t1_values = {
										value<Vec3>(cage, cage_vertex_position, triangle1[0]),
										value<Vec3>(cage, cage_vertex_position, triangle1[1]),
										value<Vec3>(cage, cage_vertex_position, triangle1[2])};

									const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2],
																				  face_vertices_[3]};
									const std::vector<Vec3> t2_values = {
										value<Vec3>(cage, cage_vertex_position, triangle2[0]),
										value<Vec3>(cage, cage_vertex_position, triangle2[1]),
										value<Vec3>(cage, cage_vertex_position, triangle2[2])};

									const auto t1_normal = value<std::vector<Vec3>>(cage, cage_face_normal, cf)[0];
									const auto t2_normal = value<std::vector<Vec3>>(cage, cage_face_normal, cf)[1];

									// update triangle 1
									const auto t1_u0 = value<std::vector<Vec3>>(cage, cage_face_edge, cf)[0];
									const auto t1_v0 = value<std::vector<Vec3>>(cage, cage_face_edge, cf)[1];

									const auto t1_u1 = t1_values[1] - t1_values[0];
									const auto t1_v1 = t1_values[2] - t1_values[1];

									const auto area_face = (t1_u0.cross(t1_v0)).norm() * 0.5;

									double t1_sj = sqrt((t1_u1.squaredNorm()) * (t1_v0.squaredNorm()) -
														2.0 * (t1_u1.dot(t1_v1)) * (t1_u0.dot(t1_v0)) +
														(t1_v1.squaredNorm()) * (t1_u0.squaredNorm())) /
												   (sqrt8 * area_face);

									new_norm_update_ += cd.n_coords_(vidx, cage_face_idx)[0] * t1_sj * t1_normal;

									// update triangle 2
									const auto t2_u0 = value<std::vector<Vec3>>(cage, cage_face_edge, cf)[2];
									const auto t2_v0 = value<std::vector<Vec3>>(cage, cage_face_edge, cf)[3];

									const auto t2_u1 = t2_values[1] - t2_values[0];
									const auto t2_v1 = t2_values[2] - t2_values[1];

									double t2_sj = sqrt((t2_u1.squaredNorm()) * (t2_v0.squaredNorm()) -
														2.0 * (t2_u1.dot(t2_v1)) * (t2_u0.dot(t2_v0)) +
														(t2_v1.squaredNorm()) * (t2_u0.squaredNorm())) /
												   (sqrt8 * area_face);

									new_norm_update_ += cd.n_coords_(vidx, cage_face_idx)[1] * t2_sj * t2_normal;

									return true;
								});

								value<Vec3>(object, object_vertex_position, v) = new_pos_update_ + new_norm_update_;
							});

							mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
						}
					});
		}

	protected:
		void init() override
		{
			mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
				app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
			mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
			connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
				mesh_provider_, this, &SpaceDeformation<MESH>::init_mesh));

			surface_render_ = static_cast<ui::SurfaceRender<MESH>*>(
				app_.module("SurfaceRender (" + std::string{mesh_traits<MESH>::name} + ")"));

			surface_diff_pptes_ = static_cast<ui::SurfaceDifferentialProperties<MESH>*>(
				app_.module("SurfaceDifferentialProperties (" + std::string{mesh_traits<MESH>::name} + ")"));

			surface_selection_ = static_cast<ui::SurfaceSelectionPO<MESH>*>(
				app_.module("SurfaceSelectionPO (" + std::string{mesh_traits<MESH>::name} + ")"));

			surface_deformation_ = static_cast<ui::SurfaceDeformation<MESH>*>(
				app_.module("SurfaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));
		}

		void left_panel() override
		{
			imgui_mesh_selector(mesh_provider_, selected_mesh_, "Object", [&](MESH& m) {
				selected_mesh_ = &m;
				mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
			});

			if (selected_mesh_)
			{
				MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);
				Parameters& p = parameters_[selected_mesh_];

				imgui_combo_attribute<Vertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
													[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
														set_vertex_position(*selected_mesh_, attribute);
													});
				if (p.vertex_position_)
				{
					ImGui::Separator();
					ImGui::Text("Global");
					if (ImGui::Button("Generate global cage"))
					{
						generate_global_cage(*selected_mesh_, p.vertex_position_, p.nb_cage);

						p.nb_cage++;
					}

					ImGui::Separator();
					ImGui::Text("Local");
					CellsSet<MESH, Vertex>* control_set = nullptr;
					// CellsSet<MESH, Vertex>* influence_set = nullptr;

					imgui_combo_cells_set(md, control_set, "Control ",
										  [&](CellsSet<MESH, Vertex>* cs) { control_set = cs; });

					/* imgui_combo_cells_set(md, influence_set, "Influence ",
											[&](CellsSet<MESH, Vertex>* cs) { influence_set = cs; });*/

					bool newCage = false;

					/* if (influence_set && influence_set->size() > 0)
					{


						if (p.new_cage_){
							MESH* l_cage = p.list_cage_[p.nb_cage-1];

							Cage_data& cd = cage_data_[l_cage];
							cd.influence_set_ = influence_set;

							p.new_cage_ = false;
						} else {
							p.last_influence_set_ = influence_set;
						}
					}*/

					if (control_set && control_set->size() > 0)
					{
						MESH* l_cage = generate_local_cage(*selected_mesh_, p.vertex_position_, control_set, p.nb_cage);
						p.nb_cage++;

						p.new_cage_ = true;

						/* if (p.last_influence_set_ && p.last_influence_set_->size() > 0)
						{
							Cage_data& cd = cage_data_[l_cage];
							cd.influence_set_ = p.last_influence_set_;
						}*/
					}

					/*CellsSet<MESH, Vertex>* influence_set = nullptr;

					imgui_combo_cells_set(md, influence_set, "Influence ",
											[&](CellsSet<MESH, Vertex>* cs) { influence_set = cs; });

					if (influence_set && influence_set->size() > 0){
						std::cout << "ici " << newCage << std::endl;
						if (newCage){
							std::cout << "check new cage" << std::endl;


							cd.influence_set_ = influence_set;

						} else {
							p.last_influence_set_ = influence_set;
						}

					}*/

					ImGui::Separator();
					ImGui::Text("Binding");
					if (p.list_cage_.size() > 1)
					{
						imgui_mesh_selector(mesh_provider_, selected_cage_, "Cage", [&](MESH& m) {
							selected_cage_ = &m;
							mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
						});

						Cage_data& cd = cage_data_[selected_cage_];

						const std::string& cage_name = mesh_provider_->mesh_name(*selected_cage_);

						if (cage_name.length() > 0)
						{

							/*std::string c_name = cage_name.substr(4, cage_name.length());

							const int selected_index = std::stoi(c_name);*/

							// inspired from https://github.com/ocornut/imgui/issues/1658
							const char* items[] = {"MVC", "Green"};
							static const char* current_item = "MVC";
							ImGuiComboFlags flags = ImGuiComboFlags_NoArrowButton;

							ImGuiStyle& style = ImGui::GetStyle();
							float w = ImGui::CalcItemWidth();
							float spacing = style.ItemInnerSpacing.x;
							float button_sz = ImGui::GetFrameHeight();
							ImGui::PushItemWidth(w - spacing * 2.0f - button_sz * 2.0f);
							if (ImGui::BeginCombo("##custom combo", current_item, ImGuiComboFlags_NoArrowButton))
							{
								for (int n = 0; n < IM_ARRAYSIZE(items); n++)
								{
									bool is_selected = (current_item == items[n]);
									if (ImGui::Selectable(items[n], is_selected))
										current_item = items[n];
									if (is_selected)
										ImGui::SetItemDefaultFocus();
								}
								ImGui::EndCombo();
							}

							if (ImGui::Button("Bind object"))
							{

								if (current_item == "MVC")
								{
									if (!cd.local_def)
									{
										bind_object_mvc(*selected_mesh_, p.vertex_position_, *selected_cage_,
														cd.cage_vertex_position_);
									}
									else
									{
										bind_local_mvc(*selected_mesh_, p.vertex_position_, *selected_cage_,
													   cd.cage_vertex_position_);

										// p.last_influence_set_ = nullptr;
									}
								}

								else if (current_item == "Green")
								{
									if (!cd.local_def)
									{
										bind_object_green(*selected_mesh_, p.vertex_position_, *selected_cage_,
														  cd.cage_vertex_position_);
									}
									else
									{
										bind_local_green(*selected_mesh_, p.vertex_position_, *selected_cage_,
														 cd.cage_vertex_position_);

										// p.last_influence_set_ = nullptr;
									}
								}
								else
								{
									std::cout << "not available yet" << std::endl;
								}
							}
						}
					}
				}
			}
		}

	private:
		MESH* selected_mesh_;
		MESH* selected_cage_;
		std::unordered_map<const MESH*, Parameters> parameters_;

		std::unordered_map<const MESH*, Cage_data> cage_data_;
		std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
		MeshProvider<MESH>* mesh_provider_;
		SurfaceRender<MESH>* surface_render_;
		SurfaceDifferentialProperties<MESH>* surface_diff_pptes_;
		SurfaceSelectionPO<MESH>* surface_selection_;
		SurfaceDeformation<MESH>* surface_deformation_;
	};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_CAGE_DEFORMATION_H_
