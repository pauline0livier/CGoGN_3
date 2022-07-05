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
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <boost/synapse/connect.hpp>

namespace cgogn
{

namespace ui
{

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
	value<Vec3>(m, vertex_position, vertices[0]) = bb_min;
	value<Vec3>(m, vertex_position, vertices[1]) = {bb_max[0]*1.5, bb_min[1] -1.5, bb_min[2]-1.5};
	value<Vec3>(m, vertex_position, vertices[2]) = {bb_max[0]*1.5, bb_max[1]*1.5, bb_min[2]-1.5};
	value<Vec3>(m, vertex_position, vertices[3]) = {bb_min[0]-1.5, bb_max[1]*1.5, bb_min[2]-1.5};
	value<Vec3>(m, vertex_position, vertices[4]) = {bb_max[0]*1.5, bb_min[1]-1.5, bb_max[2]*1.5};
	value<Vec3>(m, vertex_position, vertices[5]) = {bb_min[0]-1.5, bb_min[1]-1.5, bb_max[2]*1.5};
	value<Vec3>(m, vertex_position, vertices[6]) = {bb_min[0]-1.5, bb_max[1]*1.5, bb_max[2]*1.5};
	value<Vec3>(m, vertex_position, vertices[7]) = {bb_max[0]*1.5, bb_max[1]*1.5, bb_max[2]*1.5};
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
		double Bij = getAngleBetweenUnitVectors(ei, ek);
		double Bki = getAngleBetweenUnitVectors(ej, ei);

		Vec3 eiej = ei.cross(ek);
		Vec3 ejek = ek.cross(ej);
		Vec3 ekei = ej.cross(ei);

		Vec3 nij = eiej.normalized();
		Vec3 njk = ejek.normalized();
		Vec3 nki = ekei.normalized();

		double ui = (Bjk + (Bij * (nij.dot(njk))) + (Bki * (nki.dot(njk)))) / (2.f * ei.dot(njk));

		sumU += ui;

		it = phi<2, 1>(cage, it);
	} while (it != vertex);

	return (1.0f / r) * sumU;
}

template <typename MESH>

class CageDeformation : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "CageDeformation can only be used with meshes of dimension 2");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	struct Parameters
	{
		Parameters() : vertex_position_(nullptr), cage_(nullptr), cage_vertex_position_(nullptr)
		{
		}

		~Parameters()
		{
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;
		MESH* cage_;
		std::shared_ptr<Attribute<Vec3>> cage_vertex_position_;
		// Eigen::MatrixXd coords_;
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> coords_;

		std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;
	};

public:
	CageDeformation(const App& app)
		: Module(app, "CageDeformation (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr)
	{
	}
	~CageDeformation()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
	}

	void initialize_mesh_data(MESH& m)
	{
		Parameters& p = parameters_[&m];
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		p.vertex_position_ = vertex_position;
	}

	MESH* generate_cage(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		const std::string& m_name = mesh_provider_->mesh_name(m);

		MESH* cage = mesh_provider_->add_mesh("cage"); // (m_name + "_cage");
		std::shared_ptr<Attribute<Vec3>> cage_vertex_position = cgogn::add_attribute<Vec3, Vertex>(*cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*cage, cage_vertex_position);

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		create_box(*cage, cage_vertex_position.get(), md.bb_min_, md.bb_max_);

		mesh_provider_->emit_connectivity_changed(*cage);
		mesh_provider_->emit_attribute_changed(*cage, cage_vertex_position.get());

		Parameters& p = parameters_[&m];
		p.cage_ = cage;
		p.cage_vertex_position_ = cage_vertex_position;

		return cage;
	}

	void bind_object(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position, MESH& cage,
					 const std::shared_ptr<Attribute<Vec3>>& cage_vertex_position)
	{

		// createWeightedMatrix(const MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position,
		// const MESH& cage, const std::shared_ptr<Attribute<Vec3>>& cage_vertex_position);

		Parameters& p = parameters_[&object];

		std::shared_ptr<Attribute<uint32>> object_vertex_index = add_attribute<uint32, Vertex>(object, "weight_index");
		uint32 nb_vertices = 0;
		foreach_cell(object, [&](Vertex v) -> bool {
			value<uint32>(object, object_vertex_index, v) = nb_vertices++;
			return true;
		});

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

		p.coords_.resize(nbv_object, nbv_cage);

		foreach_cell(object, [&](Vertex v) -> bool {
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

					p.coords_(surface_point_idx, cage_point_idx) = mvc_value;
					dm.mark(d);

					value<bool>(cage, cage_vertex_marked, cage_vertex) = true;

					sumMVC += mvc_value;
				}
				// const CELL c(d);
			}

			//std::cout << "sumMVC " << sumMVC << std::endl; 

			float sum_lambda = 0.0;

			parallel_foreach_cell(cage, [&](Vertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(cage, cage_vertex_index, vc);

				p.coords_(surface_point_idx, cage_point_idx2) = p.coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				sum_lambda += p.coords_(surface_point_idx, cage_point_idx2);

				//std::cout << "local lambda" << p.coords_(surface_point_idx, cage_point_idx2) << " idx " << cage_point_idx2<< std::endl; 

				value<bool>(cage, cage_vertex_marked, vc) = false;

				return true;
			});

			//std::cout << "sum_lambda " << sum_lambda << std::endl; 

			return true;
		});

		// std::cout << p.coords_ << std::endl;

		p.cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				&cage, [&](Attribute<Vec3>* attribute) {
					if (p.cage_vertex_position_.get() == attribute)
					{
						std::cout << "recompute object position" << std::endl;

						/*parallel_foreach_cell(object, [&](Vertex v) -> bool {
							const Vec3& val = value<Vec3>(m, vertex_attribute, v);
							uint32 vidx = value<uint32>(m, vertex_index, v);
							vattr(vidx, 0) = val[0];
							vattr(vidx, 1) = val[1];
							vattr(vidx, 2) = val[2];
							return true;
						});*/
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
			mesh_provider_, this, &CageDeformation<MESH>::init_mesh));
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
				if (!p.cage_)
				{
					if (ImGui::Button("Generate cage"))
						generate_cage(*selected_mesh_, p.vertex_position_);
				}
				else
				{
					if (ImGui::Button("Bind object"))
						bind_object(*selected_mesh_, p.vertex_position_, *p.cage_, p.cage_vertex_position_);
				}
			}
		}
	}

private:
	MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_CAGE_DEFORMATION_H_
