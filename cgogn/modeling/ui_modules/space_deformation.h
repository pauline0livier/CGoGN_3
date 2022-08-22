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
#include <cgogn/geometry/ui_modules/graph_selection.h>
#include <cgogn/modeling/ui_modules/surface_deformation.h>
#include <cgogn/modeling/ui_modules/graph_deformation.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>

#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/angle.h>


#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_cage.h>

#include <cgogn/modeling/types/cage_deformation_tool.h>
#include <cgogn/modeling/types/handle_deformation_tool.h>
#include <cgogn/modeling/types/axis_deformation_tool.h>

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
	

template <typename MESH, typename GRAPH>

class SpaceDeformation : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "SpaceDeformation can only be used with meshes of dimension 2");

	template <typename T>
	using GraphAttribute = typename mesh_traits<GRAPH>::template Attribute<T>;

	template <typename T>
	using MeshAttribute = typename mesh_traits<MESH>::template Attribute<T>;

	using GraphVertex = typename mesh_traits<GRAPH>::Vertex;
	using GraphEdge = typename mesh_traits<GRAPH>::Edge;

	using MeshVertex = typename mesh_traits<MESH>::Vertex;
	using MeshEdge = typename mesh_traits<MESH>::Edge;
	using MeshFace = typename mesh_traits<MESH>::Face;

	enum SelectionTool
	{
		Handle = 0,
		Axis,
		Cage
	};

	struct Parameters
	{
		Parameters() : vertex_position_(nullptr), new_cage_(false), nb_cage(0)
		{
		}

		~Parameters()
		{
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<MeshAttribute<Vec3>> vertex_position_;

		std::vector<MESH*> list_cage_;

		std::vector<Vec3> init_position_;

		std::shared_ptr<MeshAttribute<Vec3>> gammaColor;

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


	// creation deformation tools 
	MESH* generate_global_cage(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position, int i)
	{
		// const std::string& m_name = mesh_provider_->mesh_name(m);
		Parameters& p = parameters_[&m];

		MESH* cage = mesh_provider_->add_mesh("cage" + std::to_string(i)); // (m_name + "_cage");
		std::shared_ptr<MeshAttribute<Vec3>> cage_vertex_position =
			cgogn::add_attribute<Vec3, MeshVertex>(*cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*cage, cage_vertex_position);

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		cgogn::modeling::create_box(*cage, cage_vertex_position.get(), md.bb_min_, md.bb_max_);

		mesh_provider_->emit_connectivity_changed(*cage);
		mesh_provider_->emit_attribute_changed(*cage, cage_vertex_position.get());

		std::shared_ptr<MeshAttribute<uint32>> position_indices =
			cgogn::add_attribute<uint32, MeshVertex>(*cage, "position_indices");
		cgogn::modeling::set_attribute_position_indices(*cage, position_indices.get()); 

		// Cage_data& cd = cage_data_[cage];

		p.list_cage_.push_back(cage);
		// p.list_weights_.push_back(Weights());
		// cd.cage_vertex_position_ = cage_vertex_position;

		ui::View* v1 = app_.current_view();

		surface_render_->set_vertex_position(*v1, *cage, cage_vertex_position);
		surface_render_->set_render_faces(*v1, *cage, false);

		return cage;
	}

	MESH* generate_local_cage(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position,
							  CellsSet<MESH, MeshVertex>* control_set, int i)
	{

		MESH* l_cage = mesh_provider_->add_mesh("cage" + std::to_string(i));

		std::shared_ptr<MeshAttribute<Vec3>> l_cage_vertex_position =
			cgogn::add_attribute<Vec3, MeshVertex>(*l_cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*l_cage, l_cage_vertex_position);

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 l_min = {bb_max[0], bb_max[1], bb_max[2]};
		Vec3 l_max = {bb_min[0], bb_min[1], bb_min[2]};

		Parameters& p = parameters_[&m];

		control_set->foreach_cell([&](MeshVertex v) {
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

		cgogn::modeling::CageDeformationTool<MESH> cdt; 

		cdt.create_space_tool(l_cage, l_cage_vertex_position.get(), l_min, l_max); 

		mesh_provider_->emit_connectivity_changed(*l_cage);
		mesh_provider_->emit_attribute_changed(*l_cage, l_cage_vertex_position.get());

		// Cage_data& cd = cage_data_[l_cage];

		p.list_cage_.push_back(l_cage);

		// cd.cage_vertex_position_ = l_cage_vertex_position;
		// cd.local_def = true;
		// cd.control_set_ = control_set;

		View* v1 = app_.current_view();

		surface_render_->set_vertex_position(*v1, *l_cage, l_cage_vertex_position);
		surface_render_->set_render_faces(*v1, *l_cage, false);

		surface_selection_->set_vertex_position(*l_cage, l_cage_vertex_position);

		std::shared_ptr<MeshAttribute<Vec3>> l_cage_vertex_normal =
			cgogn::add_attribute<Vec3, MeshVertex>(*l_cage, "normal");

		surface_diff_pptes_->compute_normal(*l_cage, l_cage_vertex_position.get(), l_cage_vertex_normal.get());

		MESH* i_cage = mesh_provider_->add_mesh("i_cage" + std::to_string(i));

		std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
			cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "position");

		mesh_provider_->set_mesh_bb_vertex_position(*i_cage, i_cage_vertex_position);

		Vec3 i_center = (l_min + l_max) / Scalar(2);
		Vec3 bb_min_ = ((l_min - i_center) * 1.5f) + i_center;
		Vec3 bb_max_ = ((l_max - i_center) * 1.5f) + i_center;

		// cd.center_ctrl_cage_ = i_center;
		cdt.set_center_control_cage(i_center); 
		cdt.set_influence_cage(i_cage, i_cage_vertex_position.get(), bb_min_, bb_max_); 
		//cgogn::modeling::create_box(*i_cage, i_cage_vertex_position.get(), bb_min_, bb_max_);

		mesh_provider_->emit_connectivity_changed(*i_cage);
		mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

		// cd.influence_cage_ = i_cage;
		// cd.i_cage_vertex_position_ = i_cage_vertex_position;

		surface_render_->set_vertex_position(*v1, *i_cage, i_cage_vertex_position);
		surface_render_->set_render_faces(*v1, *i_cage, false);

		CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();
		cdt.influence_area_ = &i_set; 

		//  cd.smoothing_factor_ = 0.3;

		// computeCagesAttenuationFactor(*l_cage);

		//cdt.update_influence_area(m, vertex_position.get()); 

		mesh_provider_->emit_cells_set_changed(m, cdt.influence_area_);

		return l_cage;
	}

	GRAPH* generate_handle(const Vec3& center1, const Vec3& center2)
	{

		GRAPH* handle_ = graph_provider_->add_mesh("handle");

		auto handle_vertex_position_ = add_attribute<Vec3, GraphVertex>(*handle_, "position");
		auto handle_vertex_radius_ = add_attribute<Scalar, GraphVertex>(*handle_, "radius");

		cgogn::modeling::create_handle(*handle_, handle_vertex_position_.get(), handle_vertex_radius_.get(), center1,
									   center2);

		graph_provider_->emit_connectivity_changed(*handle_);
		graph_provider_->emit_attribute_changed(*handle_, handle_vertex_position_.get());
		graph_provider_->emit_attribute_changed(*handle_, handle_vertex_radius_.get());

		// Handle_data& hd = handle_data_[handle_];
		// hd.handle_vertex_position_ = handle_vertex_position_;

		View* v1 = app_.current_view();

		graph_render_->set_vertex_position(*v1, *handle_, handle_vertex_position_);

		graph_selection_->set_vertex_position(*handle_, handle_vertex_position_);

		return handle_;
	}

	GRAPH* generate_axis(const std::vector<Vec3>& vertices_positions)
	{

		GRAPH* axis_ = graph_provider_->add_mesh("axis");

		auto axis_vertex_position_ = add_attribute<Vec3, GraphVertex>(*axis_, "position");
		auto axis_vertex_radius_ = add_attribute<Scalar, GraphVertex>(*axis_, "radius");

		cgogn::modeling::create_axis(*axis_, axis_vertex_position_.get(), axis_vertex_radius_.get(),
									 vertices_positions);

		graph_provider_->emit_connectivity_changed(*axis_);
		graph_provider_->emit_attribute_changed(*axis_, axis_vertex_position_.get());
		graph_provider_->emit_attribute_changed(*axis_, axis_vertex_radius_.get());

		// Handle_data& hd = handle_data_[axis_];
		// hd.handle_vertex_position_ = axis_vertex_position_;

		View* v1 = app_.current_view();

		graph_render_->set_vertex_position(*v1, *axis_, axis_vertex_position_);

		graph_selection_->set_vertex_position(*axis_, axis_vertex_position_);

		return axis_;
	}



	/*float cageInfluenceDistance(const int object_point_idx, Cage_data& cd, const int& nbf_cage, const int& nbv_cage)
	{

		MESH* i_cage = cd.influence_cage_;

		std::shared_ptr<MeshAttribute<uint32>> i_cage_vertex_index = get_attribute<uint32, MeshVertex>(*i_cage, "weight_index");

		std::shared_ptr<MeshAttribute<uint32>> i_cage_face_index = get_attribute<uint32, MeshFace>(*i_cage, "face_index");

		float r = 1.0;
		foreach_cell(*i_cage, [&](MeshFace fc) -> bool {
			uint32 i_cage_face_idx = value<uint32>(*i_cage, i_cage_face_index, fc);

			std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*i_cage, fc);

			// triangle 1
			const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3], face_vertices_[0]};
			const std::vector<uint32> t1_index = {value<uint32>(*i_cage, i_cage_vertex_index, triangle1[0]),
												  value<uint32>(*i_cage, i_cage_vertex_index, triangle1[1]),
												  value<uint32>(*i_cage, i_cage_vertex_index, triangle1[2])};

			const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2], face_vertices_[3]};
			const std::vector<uint32> t2_index = {value<uint32>(*i_cage, i_cage_vertex_index, triangle2[0]),
												  value<uint32>(*i_cage, i_cage_vertex_index, triangle2[1]),
												  value<uint32>(*i_cage, i_cage_vertex_index, triangle2[2])};

			r *= (1.f - (cd.coords_(object_point_idx, t1_index[0]) + cd.coords_(object_point_idx, t1_index[1]) +
						 cd.coords_(object_point_idx, t1_index[2])));

			r *= (1.f - (cd.coords_(object_point_idx, t2_index[0]) + cd.coords_(object_point_idx, t2_index[1]) +
						 cd.coords_(object_point_idx, t2_index[2])));

			return true;
		});

		r /= std::pow((1.0 - (3.0 / nbv_cage)), nbf_cage);

		if (r > 1.0f)
		{
			r = 1.0f;
		}

		return r;
	}*/

	/*void compute_attenuation(MESH& object, Cage_data& cd)
	{

		MESH* i_cage = cd.influence_cage_;

		const float h = cd.smoothing_factor_;
		assert((h <= 1.0f && h >= 0.0f) || !"Cage's attenuation factor must be computed!");

		std::shared_ptr<MeshAttribute<uint32>> i_cage_face_index = add_attribute<uint32, MeshFace>(*i_cage, "face_index");
		uint32 nb_faces_cage = 0;
		foreach_cell(*i_cage, [&](MeshFace f) -> bool {
			value<uint32>(*i_cage, i_cage_face_index, f) = nb_faces_cage++;
			return true;
		});

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index = get_attribute<uint32, MeshVertex>(object, "weight_index");

		uint32 nbf_cage = 2 * nb_cells<MeshFace>(*i_cage);
		uint32 nbv_cage = nb_cells<MeshVertex>(*i_cage);

		cd.influence_set_->foreach_cell([&](MeshVertex v) {
			uint32 surface_point_idx = value<uint32>(object, object_vertex_igenerate_local_ca
			cd.attenuation_(surface_point_idx) =
				(i_dist < h) ? (0.5f * ((float)sin(M_PI * ((i_dist / h) - 0.5f))) + 0.5f) : 1.0f;
		});
	}*/

	/*void displayGammaColor(MESH& object)
	{

		std::shared_ptr<MeshAttribute<Vec3>> gamma_color = cgogn::add_attribute<Vec3, MeshVertex>(object, "color_gamma");
		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index = get_attribute<uint32, MeshVertex>(object, "weight_index");

		Parameters& p = parameters_[&object];

		foreach_cell(object, [&](MeshVertex v) -> bool {
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			float gamma_value = 0.f;

			for (int i = 0; i < p.list_cage_.size(); i++)
			{
				const MESH* current_cage = p.list_cage_[i];
				Cage_data& cd = cage_data_[current_cage];

				gamma_value += cd.attenuation_(surface_point_idx);
			}

			Vec3 color;
			if (gamma_value <= 0.f)
				color = {0.f, 0.f, 0.f};
			else
			{
				gamma_value = (gamma_value > 1.f) ? 1.f : gamma_value;
				color[0] = (1.f + std::sin((gamma_value - 0.5f) * M_PI)) / 2.f;
				color[1] = std::sin((gamma_value + 2.f) * M_PI);
				color[2] = (1.f + std::sin((gamma_value + 0.5f) * M_PI)) / 2.f;
			}

			value<Vec3>(object, gamma_color, v) = color;

			return true;
		});

		View* v1 = app_.current_view();
		surface_render_->set_vertex_color(*v1, object, gamma_color);
	}*/

	/*float cageNormalizedDistance(MESH& ctrl_cage, const CageVolList& l_ct)
	{

		Cage_data& cd = cage_data_[&ctrl_cage];

		float sum_coeff = 0.f;
		float ct_coeff, tmp;
		for (Vol_CageIndex cf : l_ct)
		{
			CageTool& cage = m_cageTools[cf.first];
			if (cage.matricesData.pendingTaskState < CageMatricesSet::COMPUTE_GAMMAS)
			{
				sum_coeff += (tmp = cage.matricesData.gamma(cf.second));
			}
			else
			{
				sum_coeff += (tmp =
								  // influenceDistance(ptMatrixIndex,cage,cagesMap)
							  cagePointGamma(cf.second, cage));
			}
			if (cage.cageVol.dart == ct.cageVol.dart)
				ct_coeff = tmp;
		}

		// return (pointGammaCage(ptMatrixIndex,ct,cagesMap)/sum_dinf);
		return ct_coeff / sum_coeff;
	}*/

	/*void computeCageMixFactor(MESH& object, MESH& ctrl_cage)
	{
		Cage_data& cd = cage_data_[&ctrl_cage];

		uint32 nbv_object = nb_cells<MeshVertex>(object);
		cd.mixing_factor_.resize(nbv_object);
		cd.mixing_factor_.setZero();

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index = get_attribute<uint32, MeshVertex>(object, "weight_index");
		std::shared_ptr<MeshAttribute<Vec3>> object_vertex_position = get_attribute<Vec3, MeshVertex>(object, "position");

		if (!cd.b_simple_mix)
		{
			cd.influence_set_->foreach_cell([&](MeshVertex v) {
				const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
				uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

				float dnorm = 1.0f; // cageNormalizedDistance(cage, m_objectVertexCagesList[cont_i]);

				cd.mixing_factor_(surface_point_idx) = dnorm;
			});
		}
		else
		{
			cd.influence_set_->foreach_cell([&](MeshVertex v) {
				const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
				uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

				// int n = m_objectVertexCagesList[cont_i].size();
				int n = 1; // only one cage so far

				cd.mixing_factor_(surface_point_idx) = 1.0f / ((float)n);
			});
		}
	}*/

	/*void computeWeights(MESH& ctrl_cage)
	{

		Cage_data& cd = cage_data_[&ctrl_cage];
		const MESH* i_cage = cd.influence_cage_;

		std::shared_ptr<MeshAttribute<uint32>> ctrl_cage_vertex_index =
			get_attribute<uint32, MeshVertex>(ctrl_cage, "weight_index");
		std::shared_ptr<MeshAttribute<Vec3>> ctrl_cage_vertex_position = get_attribute<Vec3, MeshVertex>(ctrl_cage, "position");

		std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position = get_attribute<Vec3, MeshVertex>(*i_cage, "position");
		std::shared_ptr<MeshAttribute<uint32>> i_cage_vertex_index = get_attribute<uint32, MeshVertex>(*i_cage, "weight_index");

		std::shared_ptr<MeshAttribute<bool>> i_cage_vertex_marked = get_attribute<bool, MeshVertex>(*i_cage, "marked_vertex");

		uint32 nbv_object = nb_cells<MeshVertex>(ctrl_cage);
		uint32 nbv_cage = nb_cells<MeshVertex>(*i_cage);

		cd.ctrl_cage_coords_.resize(nbv_object, nbv_cage);

		parallel_foreach_cell(ctrl_cage, [&](MeshVertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(ctrl_cage, ctrl_cage_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(ctrl_cage, ctrl_cage_vertex_index, v);

			DartMarker dm(*i_cage);
			float sumMVC = 0.0;
			for (Dart d = i_cage->begin(), end = i_cage->end(); d != end; d = i_cage->next(d))
			{
				MeshVertex cage_vertex = CMap2::Vertex(d);
				bool vc_marked = value<bool>(*i_cage, i_cage_vertex_marked, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{

					const Vec3& cage_point = value<Vec3>(*i_cage, i_cage_vertex_position, cage_vertex);
					uint32 cage_point_idx = value<uint32>(*i_cage, i_cage_vertex_index, cage_vertex);

					float r = (cage_point - surface_point).norm();
					Vec3 e = (cage_point - surface_point).normalized();

					float s = 0.0f;

					foreach_incident_face(*i_cage, cage_vertex, [&](MeshFace f) -> bool {
						Dart dart = f.dart;

						Vec3 p_prev = value<Vec3>(*i_cage, i_cage_vertex_position, CMap2::Vertex(phi_1(*i_cage, dart)));
						Vec3 p_next = value<Vec3>(*i_cage, i_cage_vertex_position, CMap2::Vertex(phi1(*i_cage, dart)));

						Vec3 e_prev = (p_prev - surface_point).normalized();
						Vec3 e_next = (p_next - surface_point).normalized();

						Vec3 n_pv = e_prev.cross(e);
						n_pv.normalize();
						Vec3 n_vn = e.cross(e_next);
						n_vn.normalize();
						Vec3 n_np = e_next.cross(e_prev);
						n_np.normalize();

						float a_pv = getAngleBetweenUnitVectors(e_prev, e);
						float a_vn = getAngleBetweenUnitVectors(e, e_next);
						float a_np = getAngleBetweenUnitVectors(e_next, e_prev);

						float u =
							(a_np + (a_vn * (n_vn.dot(n_np))) + (a_pv * (n_pv.dot(n_np)))) / ((2.0f * e.dot(n_np)));

						s += u;

						return true;
					});

					float res = (1 / r) * s;

					cd.ctrl_cage_coords_(surface_point_idx, cage_point_idx) = res;
					dm.mark(d);

					value<bool>(*i_cage, i_cage_vertex_marked, cage_vertex) = true;

					sumMVC += res;
				}
			}

			parallel_foreach_cell(*i_cage, [&](MeshVertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(*i_cage, i_cage_vertex_index, vc);

				cd.ctrl_cage_coords_(surface_point_idx, cage_point_idx2) =
					cd.ctrl_cage_coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				value<bool>(*i_cage, i_cage_vertex_marked, vc) = false;

				return true;
			});

			return true;
		});
	}*/

	/*void computeCagesAttenuationFactor(MESH& ctrl_cage)
	{
		Cage_data& cd = cage_data_[&ctrl_cage];

		// Dart d= ctrl_cage.dart;
		float h = 1.0f;
generate_local_ca
		uint32 nbf_cage = 2 * nb_cells<MeshFace>(ctrl_cage); // same for influence & ctrl cage
		uint32 nbv_cage = nb_cells<MeshVertex>(ctrl_cage);

		if (cd.m_hFactor < 0.f)
		{
			MESH* i_cage = cd.influence_cage_;

			std::shared_ptr<MeshAttribute<uint32>> ctrl_vertex_index =
				get_attribute<uint32, MeshVertex>(ctrl_cage, "weight_index");
			std::shared_ptr<MeshAttribute<Vec3>> ctrl_vertex_position = get_attribute<Vec3, MeshVertex>(ctrl_cage, "position");

			std::shared_ptr<MeshAttribute<Vec3>> i_vertex_position = get_attribute<Vec3, MeshVertex>(*i_cage, "position");
			std::shared_ptr<MeshAttribute<uint32>> i_vertex_index = get_attribute<uint32, MeshVertex>(*i_cage, "weight_index");

			std::shared_ptr<MeshAttribute<bool>> i_vertex_marked = get_attribute<bool, MeshVertex>(*i_cage, "marked_vertex");

			cd.ctrl_cage_coords_.resize(nbv_cage, nbv_cage);
			cd.ctrl_cage_coords_.setZero();

			parallel_foreach_cell(ctrl_cage, [&](MeshVertex v) -> bool {
				const Vec3& ctrl_point = value<Vec3>(ctrl_cage, ctrl_vertex_position, v);
				uint32 ctrl_point_idx = value<uint32>(ctrl_cage, ctrl_vertex_index, v);

				DartMarker dm(*i_cage);
				float sumMVC = 0.0;
				for (Dart d = i_cage->begin(), end = i_cage->end(); d != end; d = i_cage->next(d))
				{
					MeshVertex i_vertex = CMap2::Vertex(d);
					bool vi_marked = value<bool>(*i_cage, i_vertex_marked, i_vertex);

					if (!dm.is_marked(d) && !vi_marked)
					{

						const Vec3& i_point = value<Vec3>(*i_cage, i_vertex_position, i_vertex);
						uint32 i_point_idx = value<uint32>(*i_cage, i_vertex_index, i_vertex);

						float mvc_value = compute_mvc(ctrl_point, d, *i_cage, i_point, i_vertex_position.get());

						cd.ctrl_cage_coords_(ctrl_point_idx, i_point_idx) = mvc_value;
						dm.mark(d);

						value<bool>(*i_cage, i_vertex_marked, i_vertex) = true;

						sumMVC += mvc_value;
					}
				}

				parallel_foreach_cell(*i_cage, [&](MeshVertex vc) -> bool {
					uint32 i_point_idx = value<uint32>(*i_cage, i_vertex_index, vc);

					cd.ctrl_cage_coords_(ctrl_point_idx, i_point_idx) /= sumMVC;

					value<bool>(*i_cage, i_vertex_marked, vc) = false;

					return true;
				});

				return true;
			});

			std::shared_ptr<MeshAttribute<uint32>> i_face_index = get_attribute<uint32, MeshFace>(*i_cage, "face_index");

			parallel_foreach_cell(ctrl_cage, [&](MeshVertex v) -> bool {
				uint32 ctrl_point_idx = value<uint32>(ctrl_cage, ctrl_vertex_index, v);
				float r = 1.0f;

				foreach_cell(*i_cage, [&](MeshFace fc) -> bool {
					uint32 i_face_idx = value<uint32>(*i_cage, i_face_index, fc);

					std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*i_cage, fc);

					// triangle 1
					const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3],
																  face_vertices_[0]};
					const std::vector<uint32> t1_index = {value<uint32>(*i_cage, i_vertex_index, triangle1[0]),
														  value<uint32>(*i_cage, i_vertex_index, triangle1[1]),
														  value<uint32>(*i_cage, i_vertex_index, triangle1[2])};

					const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2],
																  face_vertices_[3]};

					const std::vector<uint32> t2_index = {value<uint32>(*i_cage, i_vertex_index, triangle2[0]),
														  value<uint32>(*i_cage, i_vertex_index, triangle2[1]),
														  value<uint32>(*i_cage, i_vertex_index, triangle2[2])};

					r *= (1.f - (cd.ctrl_cage_coords_(ctrl_point_idx, t1_index[0]) +
								 cd.ctrl_cage_coords_(ctrl_point_idx, t1_index[1]) +
								 cd.ctrl_cage_coords_(ctrl_point_idx, t1_index[2])));

					r *= (1.f - (cd.ctrl_cage_coords_(ctrl_point_idx, t2_index[0]) +
								 cd.ctrl_cage_coords_(ctrl_point_idx, t2_index[1]) +
								 cd.ctrl_cage_coords_(ctrl_point_idx, t2_index[2])));

					return true;
				});

				r /= std::pow((1.f - (3.f / nbv_cage)), nbf_cage);

				if (r < h)
					h = r;

				return true;
			});
		}
		else
		{
			h = (cd.m_hFactor > 1.0f) ? 1.0f : cd.m_hFactor;
		}

		h = (1.f - h) / 3.f;

		cd.smoothing_factor_ = (h < 0.f) ? 0.f : h;
	}*/

	/*void computeGlobalDeformationMatrix(MESH& object, MESH& ctrl_cage)
	{
		Cage_data& cd = cage_data_[&ctrl_cage];

		uint32 nbv_cage = nb_cells<MeshVertex>(ctrl_cage);
		uint32 nbv_object = nb_cells<MeshVertex>(object);

		cd.global_matrix_.resize(nbv_object, nbv_cage);
		cd.global_matrix_.setZero();

		cd.global_matrix_comp_.resize(nbv_object, nbv_cage);
		cd.global_matrix_comp_.setZero();

		// along with the vectors that will contained the cages vertices position (2D) and initial postions

		for (int dim = 0; dim < DIM; ++dim)
		{
			m_toolPos[dim] = Eigen::VectorXf::Zero(nbCageVertex);
			m_toolInitPos[dim] = Eigen::VectorXf::Zero(nbCageVertex);
		}


			CageTool& cage = m_cageTools[volCont_i];
			for (int dim = 0; dim < DIM; ++dim)
			{
				int index = cage.start_index;
				int cage_index = 0;
				Dart n = cage.cageVol.dart;
				//-- we fill the postions (and intial positions) cage vertices containing vector
				//-- accordingly
				for (MeshVertex vertex : verticesIncidentToVolume3(*m_cagesMap, cage.cageVol))
				{
					m_toolPos[dim](index) = m_cageInf[vertex.dart][dim];
					m_toolInitPos[dim](index) = cage.matricesData.initialPostions[dim](cage_index);

					++index;
					++cage_index;
				}
			}

			int innerMatIndex = 0;
			//-- for each space point contained by this the current cage
			for (unsigned int sp_index : cage.containedSpacePoints)
			{
				unsigned int globalMatrixIndex = m_verticesGlobalMatrixIndex[sp_index];
				//-- -- we compute its associated definitive coefficent,
				//-- -- complementary coefficent and mixing factors
				float mixFactor = cage.matricesData.mixingFactor(innerMatIndex);
				float gamma = cage.matricesData.gamma(innerMatIndex);

				//-- -- and fill accordingly the global defromation matrix
				//-- -- (and complementary matrix)
				for (int cageIndex = 0; cageIndex < cage.inf_nb_vertex; ++cageIndex)
				{
					m_globalMatrix(globalMatrixIndex, cage.start_index + cageIndex) =
						mixFactor * gamma * cage.matricesData.coeffs(innerMatIndex, cageIndex);

					m_globalMatrix_comp(globalMatrixIndex, cage.start_index + cageIndex) =
						mixFactor * (1.0f - gamma) * cage.matricesData.coeffs(innerMatIndex, cageIndex);
				}

				++innerMatIndex;
			}
		}*/

public:
	void set_selected_mesh(MESH& m)
	{
		selected_mesh_ = &m;
	}

	void set_vertex_position(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		p.vertex_position_ = vertex_position;

		uint32 nbv_m = nb_cells<MeshVertex>(m);
		p.init_position_.resize(nbv_m);

		std::shared_ptr<MeshAttribute<uint32>> vertex_index = get_attribute<uint32, MeshVertex>(m, "position_indices");

		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(m, vertex_position, v);
			uint32 surface_point_idx = value<uint32>(m, vertex_index, v);

			p.init_position_[surface_point_idx] = surface_point;
			return true;
		});
	}

	

	
	/*void bind_local_mvc(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position, MESH& cage,
						const std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position)
	{

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index = get_attribute<uint32, MeshVertex>(object, "weight_index");

		std::shared_ptr<MeshAttribute<uint32>> cage_vertex_index = add_attribute<uint32, MeshVertex>(cage, "weight_index");
		uint32 nb_vertices_cage = 0;
		foreach_cell(cage, [&](MeshVertex v) -> bool {
			value<uint32>(cage, cage_vertex_index, v) = nb_vertices_cage++;
			return true;
		});

		std::shared_ptr<MeshAttribute<bool>> cage_vertex_marked = add_attribute<bool, MeshVertex>(cage, "marked_vertex");
		parallel_foreach_cell(cage, [&](MeshVertex v) -> bool {
			value<bool>(cage, cage_vertex_marked, v) = false;
			return true;
		});

		uint32 nbv_object = nb_cells<MeshVertex>(object);
		uint32 nbv_cage = nb_cells<MeshVertex>(cage);

		Parameters& p = parameters_[&object];
		Cage_data& cd = cage_data_[&cage];

		cd.coords_.resize(nbv_object, nbv_cage);

		cd.influence_set_->foreach_cell([&](MeshVertex v) {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			DartMarker dm(cage);
			float sumMVC = 0.0;
			for (Dart d = cage.begin(), end = cage.end(); d != end; d = cage.next(d))
			{
				MeshVertex cage_vertex = CMap2::Vertex(d);
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

			parallel_foreach_cell(cage, [&](MeshVertex vc) -> bool {
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
				&cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cd.cage_vertex_position_.get() == attribute)
					{

						std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(object, "weight_index");

						std::shared_ptr<MeshAttribute<uint32>> cage_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(cage, "weight_index");

						cd.influence_set_->foreach_cell([&](MeshVertex v) -> bool {
							uint32 vidx = value<uint32>(object, object_vertex_index, v);

							Vec3 new_pos_ = {0.0, 0.0, 0.0};

							foreach_cell(cage, [&](MeshVertex cv) -> bool {
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
	}*/

	/*/

	/*void bind_influence_cage_mvc(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
								 MESH& ctrl_cage)
	{

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index = get_attribute<uint32, MeshVertex>(object, "weight_index");

		Cage_data& cd = cage_data_[&ctrl_cage];
		MESH* i_cage = cd.influence_cage_;

		std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position = get_attribute<Vec3, MeshVertex>(*i_cage, "position");

		std::shared_ptr<MeshAttribute<uint32>> i_cage_vertex_index = get_attribute<uint32, MeshVertex>(*i_cage, "weight_index");

		std::shared_ptr<MeshAttribute<bool>> i_cage_vertex_marked = get_attribute<bool, MeshVertex>(*i_cage, "marked_vertex");

		uint32 nbv_object = nb_cells<MeshVertex>(object);
		uint32 nbv_cage = nb_cells<MeshVertex>(*i_cage);

		Parameters& p = parameters_[&object];

		cd.coords_.resize(nbv_object, nbv_cage);
		cd.coords_.setZero();

		cd.influence_set_->foreach_cell([&](MeshVertex v) {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			DartMarker dm(*i_cage);
			float sumMVC = 0.0;

			for (Dart d = i_cage->begin(), end = i_cage->end(); d != end; d = i_cage->next(d))
			{
				MeshVertex cage_vertex = CMap2::Vertex(d);
				bool vc_marked = value<bool>(*i_cage, i_cage_vertex_marked, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{
					const Vec3& i_cage_point = value<Vec3>(*i_cage, i_cage_vertex_position, cage_vertex);
					uint32 i_cage_point_idx = value<uint32>(*i_cage, i_cage_vertex_index, cage_vertex);

					float mvc_value =
						compute_mvc(surface_point, d, *i_cage, i_cage_point, i_cage_vertex_position.get());

					cd.coords_(surface_point_idx, i_cage_point_idx) = mvc_value;

					dm.mark(d);

					value<bool>(*i_cage, i_cage_vertex_marked, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			parallel_foreach_cell(*i_cage, [&](MeshVertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(*i_cage, i_cage_vertex_index, vc);

				cd.coords_(surface_point_idx, cage_point_idx2) =
					cd.coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				value<bool>(*i_cage, i_cage_vertex_marked, vc) = false;

				return true;
			});
		});

		cd.attenuation_.resize(nbv_object);
		cd.attenuation_.setZero();

		compute_attenuation(object, cd);

		displayGammaColor(object);

		cd.cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				&ctrl_cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cd.cage_vertex_position_.get() == attribute)
					{
						MESH* i_cage = cd.influence_cage_;
						std::shared_ptr<MeshAttribute<Vec3>> cage_vertex_position =
							get_attribute<Vec3, MeshVertex>(ctrl_cage, "position");

						std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
							get_attribute<Vec3, MeshVertex>(*i_cage, "position");

						foreach_cell(ctrl_cage, [&](MeshVertex cv) -> bool {
							const Vec3& cage_point = value<Vec3>(ctrl_cage, cage_vertex_position, cv);

							value<Vec3>(*i_cage, i_cage_vertex_position, cv) =
								((cage_point - cd.center_ctrl_cage_) * 1.5) + cd.center_ctrl_cage_;

							return true;
						});

						mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

						std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(object, "weight_index");

						std::shared_ptr<MeshAttribute<uint32>> i_cage_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(*i_cage, "weight_index");

						cd.influence_set_->foreach_cell([&](MeshVertex v) -> bool {
							uint32 vidx = value<uint32>(object, object_vertex_index, v);

							float current_attenuation = cd.attenuation_(vidx);

							Vec3 new_pos_ = {0.0, 0.0, 0.0};

							foreach_cell(*i_cage, [&](MeshVertex cv) -> bool {
								const Vec3& i_cage_point = value<Vec3>(*i_cage, i_cage_vertex_position, cv);
								uint32 i_cage_point_idx = value<uint32>(*i_cage, i_cage_vertex_index, cv);

								new_pos_ += cd.coords_(vidx, i_cage_point_idx) * i_cage_point;

								return true;
							});

							const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
							new_pos_ = new_pos_ * current_attenuation;

							Vec3 current_pos = (1.0f - current_attenuation) * p.init_position_[vidx];

							value<Vec3>(object, object_vertex_position, v) = new_pos_ + current_pos;

							return true;
						});

						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					}
				});
	}*/

	/*void bind_object_green(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position, MESH& cage,
						   const std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position)
	{

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index = get_attribute<uint32, MeshVertex>(object, "weight_index");

		std::shared_ptr<MeshAttribute<uint32>> cage_vertex_index = add_attribute<uint32, MeshVertex>(cage, "weight_index");
		uint32 nb_vertices_cage = 0;
		foreach_cell(cage, [&](MeshVertex v) -> bool {
			value<uint32>(cage, cage_vertex_index, v) = nb_vertices_cage++;
			return true;
		});

		std::shared_ptr<MeshAttribute<uint32>> cage_face_index = add_attribute<uint32, MeshFace>(cage, "face_index");
		uint32 nb_faces_cage = 0;
		foreach_cell(cage, [&](MeshFace f) -> bool {
			value<uint32>(cage, cage_face_index, f) = nb_faces_cage++;
			return true;
		});

		std::shared_ptr<MeshAttribute<std::vector<Vec3>>> cage_face_normal =
			add_attribute<std::vector<Vec3>, MeshFace>(cage, "face_normal");

		std::shared_ptr<MeshAttribute<std::vector<Vec3>>> cage_face_edge =
			add_attribute<std::vector<Vec3>, MeshFace>(cage, "face_edge");

		Parameters& p = parameters_[&object];
		Cage_data& cd = cage_data_[&cage];

		uint32 nbv_object = nb_cells<MeshVertex>(object);
		uint32 nbv_cage = nb_cells<MeshVertex>(cage);

		uint32 nbf_cage = nb_cells<MeshFace>(cage); // Warning valid only for square face (1 face = 2 triangles)

		cd.coords_.resize(nbv_object, nbv_cage);
		cd.coords_.setZero();

		cd.n_coords_.resize(nbv_object, nbf_cage);
		cd.n_coords_.setZero();

		parallel_foreach_cell(object, [&](MeshVertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			foreach_cell(cage, [&](MeshFace fc) -> bool {
				uint32 cage_face_idx = value<uint32>(cage, cage_face_index, fc);

				// Dart d1 = fc.dart;

				std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(cage, fc);

				// triangle 1
				// const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[0], face_vertices_[3],
				// face_vertices_[1]};
				const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3], face_vertices_[0]};
				const std::vector<Vec3> t1_values = {value<Vec3>(cage, cage_vertex_position, triangle1[0]),
													 value<Vec3>(cage, cage_vertex_position, triangle1[1]),
													 value<Vec3>(cage, cage_vertex_position, triangle1[2])};

				Vec3 t1_normal = (cgogn::geometry::normal(t1_values[0], t1_values[1], t1_values[2])).normalized();

				const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2], face_vertices_[3]};
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
				&cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cd.cage_vertex_position_.get() == attribute)
					{

						std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(object, "weight_index");

						std::shared_ptr<MeshAttribute<uint32>> cage_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(cage, "weight_index");

						std::shared_ptr<MeshAttribute<uint32>> cage_face_index =
							cgogn::get_attribute<uint32, MeshFace>(cage, "face_index");

						std::shared_ptr<MeshAttribute<std::vector<Vec3>>> cage_face_normal =
							cgogn::get_attribute<std::vector<Vec3>, MeshFace>(cage, "face_normal");

						std::shared_ptr<MeshAttribute<std::vector<Vec3>>> cage_face_edge =
							cgogn::get_attribute<std::vector<Vec3>, MeshFace>(cage, "face_edge");

						parallel_foreach_cell(object, [&](MeshVertex v) -> bool {
							uint32 vidx = value<uint32>(object, object_vertex_index, v);

							Vec3 new_pos_update_ = {0.0, 0.0, 0.0};

							const auto sqrt8 = sqrt(8);

							foreach_cell(cage, [&](MeshVertex cv) -> bool {
								const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cv);
								uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cv);

								new_pos_update_ += cd.coords_(vidx, cage_point_idx) * cage_point;

								return true;
							});

							Vec3 new_norm_update_ = {0.0, 0.0, 0.0};
							foreach_cell(cage, [&](MeshFace cf) -> bool {
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
	}*/

	/*void bind_local_green(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position, MESH& cage,
						  const std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position)
	{

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index = get_attribute<uint32, MeshVertex>(object, "weight_index");

		std::shared_ptr<MeshAttribute<uint32>> cage_vertex_index = add_attribute<uint32, MeshVertex>(cage, "weight_index");
		uint32 nb_vertices_cage = 0;
		foreach_cell(cage, [&](MeshVertex v) -> bool {
			value<uint32>(cage, cage_vertex_index, v) = nb_vertices_cage++;
			return true;
		});

		std::shared_ptr<MeshAttribute<uint32>> cage_face_index = add_attribute<uint32, MeshFace>(cage, "face_index");
		uint32 nb_faces_cage = 0;
		foreach_cell(cage, [&](MeshFace f) -> bool {
			value<uint32>(cage, cage_face_index, f) = nb_faces_cage++;
			return true;
		});

		std::shared_ptr<MeshAttribute<std::vector<Vec3>>> cage_face_normal =
			add_attribute<std::vector<Vec3>, MeshFace>(cage, "face_normal");

		std::shared_ptr<MeshAttribute<std::vector<Vec3>>> cage_face_edge =
			add_attribute<std::vector<Vec3>, MeshFace>(cage, "face_edge");

		Parameters& p = parameters_[&object];
		Cage_data& cd = cage_data_[&cage];

		uint32 nbv_object = nb_cells<MeshVertex>(object);
		uint32 nbv_cage = nb_cells<MeshVertex>(cage);

		uint32 nbf_cage = nb_cells<MeshFace>(cage); // Warning valid only for square face (1 face = 2 triangles)

		cd.coords_.resize(nbv_object, nbv_cage);
		cd.coords_.setZero();

		cd.n_coords_.resize(nbv_object, nbf_cage);
		cd.n_coords_.setZero();

		cd.influence_set_->foreach_cell([&](MeshVertex v) {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			foreach_cell(cage, [&](MeshFace fc) -> bool {
				uint32 cage_face_idx = value<uint32>(cage, cage_face_index, fc);

				// Dart d1 = fc.dart;

				std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(cage, fc);

				const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3], face_vertices_[0]};
				const std::vector<Vec3> t1_values = {value<Vec3>(cage, cage_vertex_position, triangle1[0]),
													 value<Vec3>(cage, cage_vertex_position, triangle1[1]),
													 value<Vec3>(cage, cage_vertex_position, triangle1[2])};

				Vec3 t1_normal = (cgogn::geometry::normal(t1_values[0], t1_values[1], t1_values[2])).normalized();

				const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2], face_vertices_[3]};
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
				&cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cd.cage_vertex_position_.get() == attribute)
					{

						std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(object, "weight_index");

						std::shared_ptr<MeshAttribute<uint32>> cage_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(cage, "weight_index");

						std::shared_ptr<MeshAttribute<uint32>> cage_face_index =
							cgogn::get_attribute<uint32, MeshFace>(cage, "face_index");

						std::shared_ptr<MeshAttribute<std::vector<Vec3>>> cage_face_normal =
							cgogn::get_attribute<std::vector<Vec3>, MeshFace>(cage, "face_normal");

						std::shared_ptr<MeshAttribute<std::vector<Vec3>>> cage_face_edge =
							cgogn::get_attribute<std::vector<Vec3>, MeshFace>(cage, "face_edge");

						cd.influence_set_->foreach_cell([&](MeshVertex v) {
							uint32 vidx = value<uint32>(object, object_vertex_index, v);

							Vec3 new_pos_update_ = {0.0, 0.0, 0.0};

							const auto sqrt8 = sqrt(8);

							foreach_cell(cage, [&](MeshVertex cv) -> bool {
								const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cv);
								uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cv);

								new_pos_update_ += cd.coords_(vidx, cage_point_idx) * cage_point;

								return true;
							});

							Vec3 new_norm_update_ = {0.0, 0.0, 0.0};
							foreach_cell(cage, [&](MeshFace cf) -> bool {
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
	}*/


	

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &SpaceDeformation<MESH, GRAPH>::init_mesh));

		graph_provider_ = static_cast<ui::MeshProvider<GRAPH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<GRAPH>::name} + ")"));

		surface_render_ = static_cast<ui::SurfaceRender<MESH>*>(
			app_.module("SurfaceRender (" + std::string{mesh_traits<MESH>::name} + ")"));

		graph_render_ = static_cast<ui::SurfaceRender<GRAPH>*>(
			app_.module("SurfaceRender (" + std::string{mesh_traits<GRAPH>::name} + ")"));
		

		surface_diff_pptes_ = static_cast<ui::SurfaceDifferentialProperties<MESH>*>(
			app_.module("SurfaceDifferentialProperties (" + std::string{mesh_traits<MESH>::name} + ")"));

		surface_selection_ = static_cast<ui::SurfaceSelectionPO<MESH>*>(
			app_.module("SurfaceSelectionPO (" + std::string{mesh_traits<MESH>::name} + ")"));

		/*graph_selection_ = static_cast<ui::SurfaceSelectionPO<GRAPH>*>(
			app_.module("SurfaceSelectionPO (" + std::string{mesh_traits<GRAPH>::name} + ")"));*/

		graph_selection_ = static_cast<ui::GraphSelection<GRAPH>*>(
			app_.module("GraphSelection (" + std::string{mesh_traits<GRAPH>::name} + ")")); 

		surface_deformation_ = static_cast<ui::SurfaceDeformation<MESH>*>(
			app_.module("SurfaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));

		graph_deformation_ = static_cast<ui::GraphDeformation<GRAPH>*>(
			app_.module("GraphDeformation (" + std::string{mesh_traits<GRAPH>::name} + ")"));
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

			imgui_combo_attribute<MeshVertex, Vec3>(*selected_mesh_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<MeshAttribute<Vec3>>& attribute) {
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

				ImGui::RadioButton("Handle", reinterpret_cast<int*>(&selected_tool_), Handle);
				ImGui::SameLine();
				ImGui::RadioButton("Axis", reinterpret_cast<int*>(&selected_tool_), Axis);
				ImGui::SameLine();
				ImGui::RadioButton("Cage", reinterpret_cast<int*>(&selected_tool_), Cage);

				if (selected_tool_ == Cage)
				{
					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					imgui_combo_cells_set(md, control_set, "Control ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });

					bool newCage = false;

					if (control_set && control_set->size() > 0)
					{
						MESH* l_cage = generate_local_cage(*selected_mesh_, p.vertex_position_, control_set, p.nb_cage);

						p.nb_cage++;

						p.new_cage_ = true;
					}

					ImGui::Separator();
					ImGui::Text("Binding");

					/*if (p.list_cage_.size() > 0)
					{
						imgui_mesh_selector(mesh_provider_, selected_cage_, "Cage", [&](MESH& m) {
							selected_cage_ = &m;
							mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
						});

						Cage_data& cd = cage_data_[selected_cage_];

						const std::string& cage_name = mesh_provider_->mesh_name(*selected_cage_);

						if (cage_name.length() > 0)
						{

							
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
										bind_influence_cage_mvc(*selected_mesh_, p.vertex_position_, *selected_cage_);
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
									}
								}
								else
								{
									std::cout << "not available yet" << std::endl;
								}
							}
						}
					}*/
				}

				if (selected_tool_ == Handle)
				{
					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					imgui_combo_cells_set(md, control_set, "Control ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });

					if (control_set && control_set->size() > 0)
					{
						Vec3 center = {0.0, 0.0, 0.0};

						control_set->foreach_cell([&](MeshVertex v) {
							const Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_, v);

							center += pos;
						});

						center /= control_set->size();

						MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);
						const Vec3& bb_min = md.bb_min_;
						const Vec3& bb_max = md.bb_max_;

						Vec3 center1 = center + Vec3{bb_min[0], 0.0, 0.0}; 
						Vec3 center1bis = center + Vec3{bb_max[0], 0.0, 0.0}; 
						Vec3 center2 = (center1bis + center1)*0.5; 
						Vec3 center0 = (center1bis*0.25 + center1*0.75); 

						generate_handle(center0, center2);
						
					}
				}

				if (selected_tool_ == Axis)
				{
					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					imgui_combo_cells_set(md, control_set, "Control ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });

					if (control_set && control_set->size() > 0)
					{
						
						MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);
						const Vec3& bb_min = md.bb_min_;
						const Vec3& bb_max = md.bb_max_;

						Vec3 center = {0.0, 0.0, 0.0};

						Vec3 min_pos = bb_max; 
						Vec3 max_pos = bb_min;  

						control_set->foreach_cell([&](MeshVertex v) {
							const Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_, v);

							center += pos;

							for (int i = 0; i < 3; i++){
								if (pos[i] < min_pos[i]){
									min_pos[i] = pos[i]; 
								}

								if (pos[i] > max_pos[i]){
									max_pos[i] = pos[i]; 
								}
							}

							
						});

						center /= control_set->size();

						

						Vec3 center1 = center + Vec3{bb_min[0], 0.0, 0.0}; 
						Vec3 center1bis = center + Vec3{bb_max[0], 0.0, 0.0}; 
						Vec3 center2 = (center1bis + center1)*0.5; 
						Vec3 center0 = (center1bis*0.25 + center1*0.75); 

						double center_pos = (max_pos[1] + min_pos[1])/2.0; 

						Vec3 oneExtrem = {min_pos[0], center_pos, min_pos[2]}; 
						Vec3 otherExtrem = {min_pos[0], center_pos, max_pos[2]}; 
						std::vector<Vec3> setVertices; 
						setVertices.push_back(oneExtrem); 
						setVertices.push_back(center0); 
						setVertices.push_back(otherExtrem);

						generate_axis(setVertices);  
					}
				}
			}
		}
	}

private:
	MESH* selected_mesh_;
	MESH* selected_cage_;
	std::unordered_map<const MESH*, Parameters> parameters_;

	//std::unordered_map<const MESH*, Cage_data> cage_data_;
	//std::unordered_map<const IncidenceGraph*, Handle_data> handle_data_;

	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	MeshProvider<MESH>* mesh_provider_;
	MeshProvider<GRAPH>* graph_provider_;

	SurfaceRender<MESH>* surface_render_;
	SurfaceRender<GRAPH>* graph_render_;

	SurfaceDifferentialProperties<MESH>* surface_diff_pptes_;
	
	SurfaceSelectionPO<MESH>* surface_selection_;
	GraphSelection<GRAPH>* graph_selection_;

	SurfaceDeformation<MESH>* surface_deformation_;
	GraphDeformation<GRAPH>* graph_deformation_;

	SelectionTool selected_tool_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_CAGE_DEFORMATION_H_
