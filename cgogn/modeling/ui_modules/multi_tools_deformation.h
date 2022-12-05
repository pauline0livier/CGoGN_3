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

#ifndef CGOGN_MODULE_MULTI_TOOLS_DEFORMATION_H_
#define CGOGN_MODULE_MULTI_TOOLS_DEFORMATION_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/geometry/ui_modules/graph_selection.h>
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/geometry/ui_modules/surface_selectionPO.h>
#include <cgogn/modeling/ui_modules/graph_deformation.h>
#include <cgogn/modeling/ui_modules/surface_deformation.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/angle.h>

#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

#include <cgogn/modeling/types/axis_deformation_tool.h>
#include <cgogn/modeling/types/cage_deformation_tool.h>
#include <cgogn/modeling/types/global_cage_deformation_tool.h>
#include <cgogn/modeling/types/handle_deformation_tool.h>

#include <boost/synapse/connect.hpp>

#include <iostream>

#include <algorithm>
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

class MultiToolsDeformation : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension == 2,
				  "MultiToolsDeformation can only be used with meshes of dimension 2");

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

	enum SelectionMethod
	{
		SingleCell = 0,
		WithinSphere,
		ConnectedComponent
	};

	enum SelectingCell
	{
		VertexSelect = 0,
		EdgeSelect,
		FaceSelect
	};

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), new_cage_(false), nb_cage(0), nb_max_tool_(0), nb_tool_(0),
			  selection_method_(SingleCell), selecting_cell_(VertexSelect), selected_vertices_set_(nullptr),
			  vertex_scale_factor_(1.0), sphere_scale_factor_(10.0), object_update_(false), back_selection_(false)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_point_sprite_->set_vbos({&selected_vertices_vbo_});
		}

		~Parameters()
		{
		}

	public:
		void update_selected_vertices_vbo()
		{
			if (selected_vertices_set_)
			{
				std::vector<Vec3> selected_vertices_position;
				selected_vertices_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell([&](MeshVertex v) {
					selected_vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_, v));
				});
				rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);
			}
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		MESH* mesh_;

		std::shared_ptr<MeshAttribute<Vec3>> vertex_position_;

		std::vector<MESH*> list_cage_;

		std::vector<Vec3> init_position_;

		std::shared_ptr<MeshAttribute<Vec3>> gammaColor;

		bool new_cage_;

		std::shared_ptr<boost::synapse::connection> cells_set_connection_;

		int nb_cage;

		int nb_tool_;

		int nb_max_tool_;

		bool object_update_;

		bool back_selection_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;

		rendering::VBO selected_vertices_vbo_;

		CellsSet<MESH, MeshVertex>* selected_vertices_set_;

		SelectionMethod selection_method_;
		SelectingCell selecting_cell_;
	};

public:
	std::unordered_map<std::string, std::shared_ptr<cgogn::modeling::AxisDeformationTool<MESH>>> axis_container_;
	std::unordered_map<std::string, std::shared_ptr<cgogn::modeling::HandleDeformationTool<MESH>>> handle_container_;
	std::unordered_map<std::string, std::shared_ptr<cgogn::modeling::CageDeformationTool<MESH>>> cage_container_;
	std::unordered_map<std::string, std::shared_ptr<cgogn::modeling::GlobalCageDeformationTool<MESH>>>
		global_cage_container_;

	MultiToolsDeformation(const App& app)
		: ViewModule(app, "MultiToolsDeformation (" + std::string{mesh_traits<MESH>::name} + ")"), model_(nullptr),
		  selected_mesh_(nullptr), selected_cage_(nullptr), selected_handle_(nullptr), influence_cage_mode_(true)
	{
	}

	~MultiToolsDeformation()
	{
	}

	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.mesh_ = m;
		p.nb_max_tool_ = 100;
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m](MeshAttribute<Vec3>* attribute) {
					Parameters& p = parameters_[m];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6);
						p.update_selected_vertices_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));
		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<MeshVertex>>(
				m, [this, m](CellsSet<MESH, MeshVertex>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_vertices_set_ == set && p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
	}

	void set_model(MESH& m)
	{
		model_ = &m;
	}

	void set_vertex_position(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = float32(geometry::mean_edge_length(m, p.vertex_position_.get()) / 6); // 6 ???
			p.update_selected_vertices_vbo();
		}

		uint32 nbv_m = nb_cells<MeshVertex>(m);
		p.init_position_.resize(nbv_m);

		std::shared_ptr<MeshAttribute<uint32>> vertex_index = get_attribute<uint32, MeshVertex>(m, "vertex_index");

		parallel_foreach_cell(m, [&](MeshVertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(m, vertex_position, v);
			uint32 surface_point_idx = value<uint32>(m, vertex_index, v);

			p.init_position_[surface_point_idx] = surface_point;
			return true;
		});

		for (View* v : linked_views_)
			v->request_update();
	}

	template <typename FUNC>
	void foreach_cage(const FUNC& f)
	{
		static_assert(is_ith_func_parameter_same<FUNC, 0, MESH&>::value, "Wrong function parameter type");
		static_assert(is_ith_func_parameter_same<FUNC, 1, const std::string&>::value, "Wrong function parameter type");
		for (auto& [name, gcdt] : global_cage_container_)
			f(*(gcdt->global_cage_), name);

		for (auto& [name, cdt] : cage_container_)
			f(*(cdt->control_cage_), name);
	}

	template <typename FUNC>
	void foreach_local_cage(const FUNC& f)
	{
		static_assert(is_ith_func_parameter_same<FUNC, 0, MESH&>::value, "Wrong function parameter type");
		static_assert(is_ith_func_parameter_same<FUNC, 1, const std::string&>::value, "Wrong function parameter type");

		for (auto& [name, cdt] : cage_container_)
			f(*(cdt->control_cage_), name);
	}

	template <typename FUNC>
	void foreach_handle(const FUNC& f)
	{
		// static_assert(is_ith_func_parameter_same<FUNC, 0, GRAPH&>::value, "Wrong function parameter type");
		// static_assert(is_ith_func_parameter_same<FUNC, 1, const std::string&>::value, "Wrong function parameter
		// type");
		for (auto& [name, hdt] : handle_container_)
			f(*(hdt->control_handle_), name);
	}

	/////////////
	// SIGNALS //
	/////////////
	using cage_added = struct cage_added_ (*)(cgogn::modeling::CageDeformationTool<MESH>* c);
	using global_cage_added = struct global_cage_added_ (*)(cgogn::modeling::GlobalCageDeformationTool<MESH>* c);
	using handle_added = struct handle_added_ (*)(cgogn::modeling::HandleDeformationTool<MESH>* h);
	using axis_added = struct axis_added_ (*)(cgogn::modeling::AxisDeformationTool<MESH>* a);

private:
	// creation deformation tools
	MESH* generate_global_cage(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position, int i)
	{
		std::string cage_name = "global_cage";
		MESH* cage = mesh_provider_->add_mesh(cage_name);
		std::shared_ptr<MeshAttribute<Vec3>> cage_vertex_position =
			cgogn::add_attribute<Vec3, MeshVertex>(*cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*cage, cage_vertex_position);

		set_vertex_position(*cage, cage_vertex_position); 

		const auto [it, inserted] = global_cage_container_.emplace(
			cage_name, std::make_shared<cgogn::modeling::GlobalCageDeformationTool<MESH>>());
		cgogn::modeling::GlobalCageDeformationTool<MESH>* gcdt = it->second.get();

		if (inserted)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(m);
			Vec3 center = (md.bb_min_ + md.bb_max_) / Scalar(2);
			Vec3 bb_min = ((md.bb_min_ - center) * 1.1) + center;
			Vec3 bb_max = ((md.bb_max_ - center) * 1.1) + center;

			gcdt->create_global_cage(cage, cage_vertex_position.get(), bb_min, bb_max);

			mesh_provider_->emit_connectivity_changed(*cage);
			mesh_provider_->emit_attribute_changed(*cage, cage_vertex_position.get());

			ui::View* v1 = app_.current_view();

			surface_render_->set_vertex_position(*v1, *cage, cage_vertex_position);
			surface_render_->set_render_faces(*v1, *cage, false);

			surface_selection_->set_vertex_position(*cage, cage_vertex_position);

			surface_deformation_->set_vertex_position(*cage, cage_vertex_position);
		}

		return cage;
	}

	MESH* generate_local_cage(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position,
							  CellsSet<MESH, MeshVertex>* control_set)
	{
		int cage_number = cage_container_.size();
		std::string cage_name = "local_cage" + std::to_string(cage_number);
		MESH* l_cage = mesh_provider_->add_mesh(cage_name);

		std::shared_ptr<MeshAttribute<Vec3>> l_cage_vertex_position =
			cgogn::add_attribute<Vec3, MeshVertex>(*l_cage, "position");
		mesh_provider_->set_mesh_bb_vertex_position(*l_cage, l_cage_vertex_position);

		set_vertex_position(*l_cage, l_cage_vertex_position); 

		auto mesh_vertex_normal = get_attribute<Vec3, MeshVertex>(m, "normal");

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 local_min = {1000.0, 1000.0, 1000.0};
		Vec3 local_max = {0.0, 0.0, 0.0};

		Vec3 center = {0.0, 0.0, 0.0};
		Vec3 normal = {0.0, 0.0, 0.0};

		control_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			const Vec3& norm = value<Vec3>(m, mesh_vertex_normal, v);

			center += pos;
			normal += norm;

			for (size_t j = 0; j < 3; j++)
			{
				if (pos[j] < local_min[j])
				{
					local_min[j] = pos[j];
				}

				if (pos[j] > local_max[j])
				{
					local_max[j] = pos[j];
				}
			}
		});

		center /= control_set->size();
		normal /= control_set->size();

		Vec3 ray = normal;
		ray.normalize();

		const auto [it, inserted] =
			cage_container_.emplace(cage_name, std::make_shared<cgogn::modeling::CageDeformationTool<MESH>>());
		cgogn::modeling::CageDeformationTool<MESH>* cdt = it->second.get();

		if (inserted)
		{
			cdt->create_space_tool(l_cage, l_cage_vertex_position.get(), local_min, local_max, center, ray);

			mesh_provider_->emit_connectivity_changed(*l_cage);
			mesh_provider_->emit_attribute_changed(*l_cage, l_cage_vertex_position.get());

			View* v1 = app_.current_view();

			surface_render_->set_vertex_position(*v1, *l_cage, l_cage_vertex_position);
			surface_render_->set_render_faces(*v1, *l_cage, false);

			surface_selection_->set_vertex_position(*l_cage, l_cage_vertex_position);

			std::shared_ptr<MeshAttribute<Vec3>> l_cage_vertex_normal =
				cgogn::add_attribute<Vec3, MeshVertex>(*l_cage, "normal");

			surface_diff_pptes_->compute_normal(*l_cage, l_cage_vertex_position.get(), l_cage_vertex_normal.get());

			MESH* i_cage = mesh_provider_->add_mesh("i_cage" + cage_number);

			std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
				cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "position");

			mesh_provider_->set_mesh_bb_vertex_position(*i_cage, i_cage_vertex_position);

			Vec3 i_center = (local_min + local_max) / Scalar(2);
			Vec3 bb_min_ = ((local_min - i_center) * 2.f) + i_center;
			Vec3 bb_max_ = ((local_max - i_center) * 2.f) + i_center;

			cdt->set_center_control_cage(i_center);
			cdt->set_influence_cage(i_cage, i_cage_vertex_position.get(), bb_min_, bb_max_);

			mesh_provider_->emit_connectivity_changed(*i_cage);
			mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

			surface_render_->set_vertex_position(*v1, *i_cage, i_cage_vertex_position);
			surface_render_->set_render_faces(*v1, *i_cage, false);

			CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();
			cdt->influence_area_ = &i_set;

			//  cd.smoothing_factor_ = 0.3;

			// computeCagesAttenuationFactor(*l_cage);
			cdt->update_influence_area(m, vertex_position.get());

			mesh_provider_->emit_cells_set_changed(m, cdt->influence_area_);

			boost::synapse::emit<cage_added>(this, cdt);
		}

		return l_cage;
	}

	//
	// const std::string
	std::shared_ptr<cgogn::modeling::HandleDeformationTool<MESH>> generate_handle_tool(
		MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position, CellsSet<MESH, MeshVertex>* handle_set)
	{
		int handle_number = handle_container_.size();
		const std::string handle_name = "local_handle" + std::to_string(handle_number);
		GRAPH* handle = graph_provider_->add_mesh(handle_name);

		auto handle_vertex_position = add_attribute<Vec3, GraphVertex>(*handle, "position");
		auto handle_vertex_radius = add_attribute<Scalar, GraphVertex>(*handle, "radius");

		auto mesh_vertex_normal = get_attribute<Vec3, MeshVertex>(m, "normal");

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 local_min = {1000.f, 1000.f, 1000.f};
		Vec3 local_max = {0.f, 0.f, 0.f};

		Vec3 center = {0.0, 0.0, 0.0};
		Vec3 normal = {0.0, 0.0, 0.0};

		handle_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			const Vec3& norm = value<Vec3>(m, mesh_vertex_normal, v);

			center += pos;
			normal += norm;

			for (size_t j = 0; j < 3; j++)
			{
				if (pos[j] < local_min[j])
				{
					local_min[j] = pos[j];
				}

				if (pos[j] > local_max[j])
				{
					local_max[j] = pos[j];
				}
			}
		});

		center /= handle_set->size();
		normal /= handle_set->size();

		Vec3 ray = normal;
		ray.normalize();

		// const Vec3 handle_position = center;
		const Vec3 handle_position =
			center; //{center[0] + 0.1f * ray[0], center[1] + 0.1f * ray[1], center[2] + 0.1f * ray[2]};

		const Vec3 inner_handle_position = {center[0] - 2.f * ray[0], center[1] - 2.f * ray[1],
											center[2] - 2.f * ray[2]};

		MeshVertex closest_vertex;
		double min_dist = 1000000;
		handle_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);

			double dist = (pos - handle_position).squaredNorm();
			if (dist < min_dist)
			{
				min_dist = dist;
				closest_vertex = v;
			}
		});

		const auto [it, inserted] =
			handle_container_.emplace(handle_name, std::make_shared<cgogn::modeling::HandleDeformationTool<MESH>>());
		cgogn::modeling::HandleDeformationTool<MESH>* hdt = it->second.get();

		if (inserted)
		{

			Parameters& p = parameters_[&m];

			hdt->create_space_tool(handle, handle_vertex_position.get(), handle_vertex_radius.get(), handle_position,
								   inner_handle_position, ray, p.nb_tool_);

			hdt->set_handle_mesh_vertex(closest_vertex);

			graph_provider_->emit_connectivity_changed(*handle);
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_position.get());
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_radius.get());

			View* v1 = app_.current_view();
			graph_render_->set_vertex_position(*v1, *handle, handle_vertex_position);

			graph_selection_->set_vertex_position(*handle, handle_vertex_position);

			graph_deformation_->set_vertex_position(*handle, handle_vertex_position);

			graph_deformation_->set_displacement_normal(*handle, ray);

			auto object_geodesic = get_or_add_attribute<Scalar, MeshVertex>(*model_, "geodesic_distance");

			hdt->set_geodesic_distance(m, vertex_position);
			mesh_provider_->emit_attribute_changed(m, object_geodesic.get());

			boost::synapse::emit<handle_added>(this, hdt);

			p.nb_tool_++;
		}

		return handle_container_[handle_name];
		// return handle_name;
	}

	GRAPH* generate_axis(const MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position,
						 CellsSet<MESH, MeshVertex>* control_set)
	{
		int axis_number = axis_container_.size();
		std::string axis_name = "axis" + std::to_string(axis_number);
		GRAPH* axis = graph_provider_->add_mesh(axis_name);

		auto axis_vertex_position = add_attribute<Vec3, GraphVertex>(*axis, "position");
		auto axis_vertex_radius = add_attribute<Scalar, GraphVertex>(*axis, "radius");

		auto mesh_vertex_normal = get_attribute<Vec3, MeshVertex>(m, "normal");

		MeshData<MESH>& md = mesh_provider_->mesh_data(m);
		const Vec3& bb_min = md.bb_min_;
		const Vec3& bb_max = md.bb_max_;

		Vec3 local_min = {1000.f, 1000.f, 1000.f};
		Vec3 local_max = {0.f, 0.f, 0.f};

		Vec3 center = {0.0, 0.0, 0.0};
		Vec3 normal = {0.0, 0.0, 0.0};

		control_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			const Vec3& norm = value<Vec3>(m, mesh_vertex_normal, v);

			center += pos;
			normal += norm;

			for (size_t j = 0; j < 3; j++)
			{
				if (pos[j] < local_min[j])
				{
					local_min[j] = pos[j];
				}

				if (pos[j] > local_max[j])
				{
					local_max[j] = pos[j];
				}
			}
		});

		center /= control_set->size();
		normal /= control_set->size();

		Vec3 ray = normal;
		ray.normalize();

		Eigen::Matrix3f covariance_matrix;
		// inspired from https://gist.github.com/atandrau/847214/882418ab34737699a6b1394d3a28c66e2cc0856f
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				covariance_matrix(i, j) = 0.0;
				control_set->foreach_cell([&](MeshVertex v) {
					const Vec3& pos = value<Vec3>(m, vertex_position, v);
					covariance_matrix(i, j) += (center[i] - pos[i]) * (center[j] - pos[j]);
				});

				covariance_matrix(i, j) /= control_set->size() - 1;
			}

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> eigen_solver(
			covariance_matrix);
		Eigen::Matrix<float, 1, Eigen::Dynamic> eigen_values = eigen_solver.eigenvalues();
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> eigen_vectors = eigen_solver.eigenvectors();

		Eigen::Vector3f main_eigen_vector = cgogn::modeling::sort_eigen_vectors(eigen_values, eigen_vectors);
		main_eigen_vector.normalize();

		Vec3 main_direction = {main_eigen_vector[0], main_eigen_vector[1], main_eigen_vector[2]};

		// loop on control set to find extremum points and height
		double min_x = 1000, max_x = 0.0, max_normal = center.dot(ray);
		control_set->foreach_cell([&](MeshVertex v) {
			const Vec3& pos = value<Vec3>(m, vertex_position, v);

			double x_value = pos.dot(main_direction);
			double n_value = pos.dot(ray);

			if (x_value > max_x)
			{
				max_x = x_value;
			}

			if (x_value < min_x)
			{
				min_x = x_value;
			}

			if (n_value > max_normal)
			{
				max_normal = n_value;
			}
		});

		const Vec3 axis_center = center + 2.0 * (max_normal - center.dot(ray)) * ray;
		const Vec3 extrem_min = axis_center + (min_x - center.dot(main_direction)) * main_direction;
		const Vec3 extrem_max = axis_center + (max_x - center.dot(main_direction)) * main_direction;

		const double width = (extrem_max - extrem_min).norm() / 3.0;

		std::vector<Vec3> axis_vertices;
		axis_vertices.push_back(extrem_min);
		axis_vertices.push_back(axis_center);
		axis_vertices.push_back(extrem_max);

		const auto [it, inserted] =
			axis_container_.emplace(axis_name, std::make_shared<cgogn::modeling::AxisDeformationTool<MESH>>());
		cgogn::modeling::AxisDeformationTool<MESH>* adt = it->second.get();
		if (inserted)
		{

			adt->create_space_tool(axis, axis_vertex_position.get(), axis_vertex_radius.get(), axis_vertices);

			graph_provider_->emit_connectivity_changed(*axis);
			graph_provider_->emit_attribute_changed(*axis, axis_vertex_position.get());
			graph_provider_->emit_attribute_changed(*axis, axis_vertex_radius.get());

			View* v1 = app_.current_view();

			graph_render_->set_vertex_position(*v1, *axis, axis_vertex_position);

			graph_selection_->set_vertex_position(*axis, axis_vertex_position);

			MESH* i_cage = mesh_provider_->add_mesh("i_cage" + axis_number);

			std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
				cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "position");

			std::shared_ptr<MeshAttribute<Vec3>> i_cage_local_vertex_position =
				cgogn::add_attribute<Vec3, MeshVertex>(*i_cage, "local_position");

			std::shared_ptr<MeshAttribute<uint32>> i_cage_local_skeleton =
				cgogn::add_attribute<uint32, MeshVertex>(*i_cage, "local_skeleton");

			mesh_provider_->set_mesh_bb_vertex_position(*i_cage, i_cage_vertex_position);

			Vec3 i_center = (local_min + local_max) / Scalar(2);
			Vec3 bb_min_ = ((local_min - i_center) * 1.5f) + i_center;
			Vec3 bb_max_ = ((local_max - i_center) * 1.5f) + i_center;

			adt->set_influence_cage_axis(i_cage, i_cage_vertex_position.get(), i_cage_local_vertex_position.get(),
										 i_cage_local_skeleton.get(), bb_min_, bb_max_, main_direction, ray, width);

			mesh_provider_->emit_connectivity_changed(*i_cage);
			mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

			surface_render_->set_vertex_position(*v1, *i_cage, i_cage_vertex_position);
			surface_render_->set_render_faces(*v1, *i_cage, false);

			CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();
			adt->influence_area_ = &i_set;

			adt->update_influence_area(m, vertex_position.get());

			mesh_provider_->emit_cells_set_changed(m, adt->influence_area_);

			boost::synapse::emit<axis_added>(this, adt);
		}

		return axis;
	}

	/// BIND

	void bind_global_cage(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
						  MESH& global_cage, std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position,
						  std::string binding_type)
	{

		std::shared_ptr<cgogn::modeling::GlobalCageDeformationTool<MESH>> gcdt = global_cage_container_["global_cage"];

		if (binding_type == "MVC")
		{
			gcdt->bind_mvc(object, object_vertex_position.get());

			gcdt->cage_attribute_update_connection_ =
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					&global_cage, [&](MeshAttribute<Vec3>* attribute) {
						if (cage_vertex_position.get() == attribute)
						{

							std::shared_ptr<cgogn::modeling::GlobalCageDeformationTool<MESH>> current_gcdt =
								global_cage_container_["global_cage"];

							std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
								get_attribute<uint32, MeshVertex>(object, "vertex_index");

							std::shared_ptr<MeshAttribute<uint32>> cage_vertex_index =
								get_attribute<uint32, MeshVertex>(global_cage, "vertex_index");

							parallel_foreach_cell(object, [&](MeshVertex v) -> bool {
								uint32 vidx = value<uint32>(object, object_vertex_index, v);

								Vec3 new_pos_ = {0.0, 0.0, 0.0};

								foreach_cell(global_cage, [&](MeshVertex cv) -> bool {
									const Vec3& cage_point = value<Vec3>(global_cage, cage_vertex_position, cv);

									uint32 cage_point_idx = value<uint32>(global_cage, cage_vertex_index, cv);

									new_pos_ += current_gcdt->global_cage_coords_(vidx, cage_point_idx) * cage_point;

									return true;
								});

								value<Vec3>(object, object_vertex_position, v) = new_pos_;
								return true;
							});

							mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
						}
					});
		}

		if (binding_type == "Green")
		{
			gcdt->bind_green(object, object_vertex_position.get());

			gcdt->cage_attribute_update_connection_ =
				boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
					&global_cage, [&](MeshAttribute<Vec3>* attribute) {
						if (cage_vertex_position.get() == attribute)
						{
							std::shared_ptr<cgogn::modeling::GlobalCageDeformationTool<MESH>> current_gcdt =
								global_cage_container_["global_cage"];

							current_gcdt->update_green(object, object_vertex_position.get());

							mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
						}
					});
		}
	}

	void bind_influence_cage_mvc(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
								 MESH& control_cage, const std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position)

	{
		selected_cdt_->bind_mvc(object, object_vertex_position);
		selected_cdt_->set_up_attenuation(object, object_vertex_position);

		displayGammaColor(object);

		Parameters& p = parameters_[&object];

		selected_cdt_->cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				&control_cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cage_vertex_position.get() == attribute)
					{

						selected_cdt_->update_influence_cage_position();

						MESH* i_cage = selected_cdt_->influence_cage_;

						std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
							get_attribute<Vec3, MeshVertex>(*i_cage, "position");

						mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

						selected_cdt_->update_deformation_object_mvc(object, object_vertex_position, p.init_position_);

						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					}
				});
	}

	void bind_influence_cage_mvc_handle(MESH& object,
										const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
										GRAPH& control_handle,
										const std::shared_ptr<GraphAttribute<Vec3>>& handle_vertex_position)

	{
		selected_hdt_->bind_mvc(object, object_vertex_position);
		selected_hdt_->set_up_attenuation(object, object_vertex_position);

		displayGammaColor(object);

		Parameters& p = parameters_[&object];

		selected_hdt_->handle_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<GRAPH>::template attribute_changed_t<Vec3>>(
				&control_handle, [&](GraphAttribute<Vec3>* attribute) {
					if (handle_vertex_position.get() == attribute)
					{
						selected_hdt_->update_influence_cage_position();

						MESH* i_cage = selected_hdt_->influence_cage_;

						std::shared_ptr<MeshAttribute<Vec3>> i_cage_vertex_position =
							get_attribute<Vec3, MeshVertex>(*i_cage, "position");

						mesh_provider_->emit_attribute_changed(*i_cage, i_cage_vertex_position.get());

						selected_hdt_->update_deformation_object_mvc(object, object_vertex_position, p.init_position_);

						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					}
					else
					{
						std::cout << " no moving " << std::endl;
					}
				});
	}

	void bind_handle_influence_area(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
									CellsSet<MESH, MeshVertex>* influence_set, GRAPH& control_handle,
									const std::shared_ptr<GraphAttribute<Vec3>>& handle_vertex_position)
	{

		MeshData<MESH>& md = mesh_provider_->mesh_data(object);
		CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();
		// handle_container_[handle_name]->influence_area_ = &i_set;
		selected_hdt_->influence_area_ = &i_set;

		Parameters& p = parameters_[&object];

		// handle_container_[handle_name]->set_influence_area(object, object_vertex_position, influence_set);
		selected_hdt_->set_influence_area(object, object_vertex_position, influence_set);

		// mesh_provider_->emit_cells_set_changed(object, handle_container_[handle_name]->influence_area_);
		mesh_provider_->emit_cells_set_changed(object, selected_hdt_->influence_area_);

		displayGammaColor(*model_);

		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
			cgogn::get_attribute<uint32, MeshVertex>(object, "vertex_index");

		/*selected_hdt_->influence_area_->foreach_cell([&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(object, object_vertex_index, v);

			p.attenuation_matrix_(vidx, selected_hdt_->id_) = selected_hdt_->attenuation_[vidx];

			/*if (selected_hdt_->attenuation_[vidx] == 1.0){
				p.attenuation_matrix_(vidx, selected_hdt_->id_) = selected_hdt_->attenuation_[vidx];
			} else {
				Scalar rest = 1 - p.attenuation_matrix_.row(vidx).sum();

				if (rest > 0.0){
					p.attenuation_matrix_(vidx, selected_hdt_->id_) = std::min(rest, selected_hdt_->attenuation_[vidx]);

					std::cout << "set value" << p.attenuation_matrix_(vidx, selected_hdt_->id_) << std::endl;

					std::cout << "original value" << selected_hdt_->attenuation_[vidx] << std::endl;
				}

			//}

			/*return true;
		});*/

		selected_hdt_->handle_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<GRAPH>::template attribute_changed_t<Vec3>>(
				&control_handle, [&](GraphAttribute<Vec3>* attribute) {
					if (handle_vertex_position.get() == attribute)
					{

						const Vec3 new_deformation = selected_hdt_->get_handle_deformation();
						/*p.deformation_vector_(selected_hdt_->id_, 0) += new_deformation[0];
						p.deformation_vector_(selected_hdt_->id_, 1) += new_deformation[1];
						p.deformation_vector_(selected_hdt_->id_, 2) += new_deformation[2]; */

						std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(object, "vertex_index");

						selected_hdt_->influence_area_->foreach_cell([&](MeshVertex v) -> bool {
							uint32 vidx = value<uint32>(object, object_vertex_index, v);

							value<Vec3>(object, object_vertex_position, v) +=
								selected_hdt_->attenuation_[vidx] * new_deformation;
							/*p.attenuation_matrix_(vidx, selected_hdt_->id_)*new_deformation;  */

							return true;
						});

						// update_deformation_object(object);
						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					}
					else
					{
						std::cout << " no moving " << std::endl;
					}
				});
	}

	/*void updateDeformationObject(MESH& object)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			cgogn::get_attribute<uint32, MeshVertex>(object, "vertex_index");
		selected_hdt_->influence_area_->foreach_cell([&](MeshVertex v) -> bool {
			uint32 vidx = value<uint32>(object, object_vertex_index, v);


			return true;
		});


	}*/

	void displayGammaColor(MESH& object)
	{

		std::shared_ptr<MeshAttribute<Vec3>> gamma_color =
			cgogn::get_or_add_attribute<Vec3, MeshVertex>(object, "color_gamma");
		std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		Parameters& p = parameters_[&object];

		foreach_cell(object, [&](MeshVertex v) -> bool {
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			double gamma_value = 0.0;

			/*for (auto& [name, cdt] : cage_container_)
			{
				gamma_value += cdt->attenuation_(surface_point_idx);
			}*/

			for (auto& [name, hdt] : handle_container_)
			{
				gamma_value += hdt->attenuation_(surface_point_idx);
			}

			/*for (auto& [name, adt] : axis_container_)
			{
				gamma_value += adt->attenuation_(surface_point_idx);
			}*/

			Vec3 color;
			if (gamma_value <= 0.0)
				color = {0.0, 0.0, 0.0};
			else
			{
				gamma_value = (gamma_value > 1.0) ? 1.0 : gamma_value;
				color[0] = (1.0 + std::sin((gamma_value - 0.5) * M_PI)) / 2.0;
				color[1] = std::sin((gamma_value + 2.0) * M_PI);
				color[2] = (1.0 + std::sin((gamma_value + 0.5) * M_PI)) / 2.0;
			}

			value<Vec3>(object, gamma_color, v) = color;

			return true;
		});

		View* v1 = app_.current_view();
		surface_render_->set_vertex_color(*v1, object, gamma_color);
	}

protected:
	void init() override
	{

		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &MultiToolsDeformation<MESH, GRAPH>::init_mesh));
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

		graph_selection_ = static_cast<ui::GraphSelection<GRAPH>*>(
			app_.module("GraphSelection (" + std::string{mesh_traits<GRAPH>::name} + ")"));

		surface_deformation_ = static_cast<ui::SurfaceDeformation<MESH>*>(
			app_.module("SurfaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));

		graph_deformation_ = static_cast<ui::GraphDeformation<GRAPH>*>(
			app_.module("GraphDeformation (" + std::string{mesh_traits<GRAPH>::name} + ")"));
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		if (selected_mesh_ && view->shift_pressed())
		{

			Parameters& p = parameters_[selected_mesh_];

			p.object_update_ = true;
			if (p.vertex_position_)
			{
				rendering::GLVec3d near = view->unproject(x, y, 0.0);
				rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
				Vec3 A{near.x(), near.y(), near.z()};
				Vec3 B{far_d.x(), far_d.y(), far_d.z()};

				switch (p.selection_method_)
				{
				case SingleCell: {
					switch (p.selecting_cell_)
					{
					case VertexSelect:
						if (p.selected_vertices_set_)
						{
							std::vector<MeshVertex> picked;
							cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
							if (!picked.empty())
							{
								switch (button)
								{
								case 0:
									p.selected_vertices_set_->select(picked[0]);
									break;
								case 1:
									p.selected_vertices_set_->unselect(picked[0]);
									break;
								}
								mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
							}
						}
						break;
					}
					break;
				}
				case WithinSphere: {
					std::vector<MeshVertex> picked;
					cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
					if (!picked.empty())
					{

						CellCache<MESH> cache = geometry::within_sphere(*selected_mesh_, picked[0],
																		p.vertex_base_size_ * p.sphere_scale_factor_,
																		p.vertex_position_.get());

						switch (p.selecting_cell_)
						{
						case VertexSelect:
							if (p.selected_vertices_set_)
							{
								switch (button)
								{
								case 0:
									foreach_cell(cache, [&p](MeshVertex v) -> bool {
										p.selected_vertices_set_->select(v);
										return true;
									});

									break;
								case 1:
									foreach_cell(cache, [&p](MeshVertex v) -> bool {
										p.selected_vertices_set_->unselect(v);
										return true;
									});
									break;
								}
								if (!p.back_selection_)
								{
									mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
								}
							}
							break;
						}

						if (p.back_selection_)
						{

							CellCache<MESH> cache_back_ = geometry::within_sphere(
								*selected_mesh_, picked[1], p.vertex_base_size_ * p.sphere_scale_factor_,
								p.vertex_position_.get());

							if (p.selecting_cell_ == VertexSelect)
							{
								if (p.selected_vertices_set_)
								{
									switch (button)
									{
									case 0:
										foreach_cell(cache_back_, [&p](MeshVertex v) -> bool {
											p.selected_vertices_set_->select(v);
											return true;
										});

										break;

									case 1:
										foreach_cell(cache_back_, [&p](MeshVertex v) -> bool {
											p.selected_vertices_set_->unselect(v);
											return true;
										});
										break;
									}
									mesh_provider_->emit_cells_set_changed(*selected_mesh_, p.selected_vertices_set_);
								}
							}
						}
					}
					break;
				}
				}
			}
		}
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			const rendering::GLMat4& proj_matrix = view->projection_matrix();
			const rendering::GLMat4& view_matrix = view->modelview_matrix();

			if (p.selecting_cell_ == VertexSelect && p.selected_vertices_set_ && p.selected_vertices_set_->size() > 0 &&
				p.param_point_sprite_->attributes_initialized())
			{
				p.param_point_sprite_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
				p.param_point_sprite_->bind(proj_matrix, view_matrix);
				glDrawArrays(GL_POINTS, 0, p.selected_vertices_set_->size());
				p.param_point_sprite_->release();
			}
		}
	}

	void left_panel() override
	{
		imgui_mesh_selector(mesh_provider_, model_, "Object", [&](MESH& m) {
			model_ = &m;
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		bool need_update = false;

		selected_mesh_ = model_; 

		if (model_)
		{
			MeshData<MESH>& model_md = mesh_provider_->mesh_data(*model_);
			Parameters& model_p = parameters_[model_];

			auto object_vertex_position = get_attribute<Vec3, MeshVertex>(*model_, "position");

			set_vertex_position(*model_, object_vertex_position);

			if (model_p.vertex_position_)
			{
				ImGui::Separator();
				ImGui::Separator();
				ImGui::Text("Create tool");
				ImGui::Separator();
				ImGui::Text("Global");
				if (ImGui::Button("Generate global cage"))
				{
					generate_global_cage(*model_, model_p.vertex_position_, model_p.nb_cage);
				}

				if (global_cage_container_.size() == 1)
				{
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

					if (ImGui::Button("Bind global cage"))
					{
						std::shared_ptr<cgogn::modeling::GlobalCageDeformationTool<MESH>> gcdt =
							global_cage_container_["global_cage"];

						bind_global_cage(*model_, model_p.vertex_position_, *(gcdt->global_cage_),
										 gcdt->global_cage_vertex_position_, current_item);
					}
				}

				ImGui::Separator();
				ImGui::Text("Local");

				ImGui::RadioButton("Axis", reinterpret_cast<int*>(&selected_tool_), Axis);
				ImGui::SameLine();
				ImGui::RadioButton("Cage", reinterpret_cast<int*>(&selected_tool_), Cage);
				ImGui::SameLine();
				ImGui::RadioButton("Handle", reinterpret_cast<int*>(&selected_tool_), Handle);

				if (selected_tool_ == Cage)
				{

					ImGui::Separator();
					model_p.selection_method_ = WithinSphere;
					model_p.back_selection_ = true;

					ImGui::SliderFloat("Sphere radius", &(model_p.sphere_scale_factor_), 10.0f, 100.0f);

					if (ImGui::Button("Create control set##vertices_set"))
						model_md.template add_cells_set<MeshVertex>();
					imgui_combo_cells_set(model_md, model_p.selected_vertices_set_, "Sets", [&](CellsSet<MESH, MeshVertex>* cs) {
						model_p.selected_vertices_set_ = cs;
						model_p.update_selected_vertices_vbo();
						need_update = true;
					});

					if (model_p.selected_vertices_set_)
					{
						ImGui::Text("(nb elements: %d)", model_p.selected_vertices_set_->size());
						if (ImGui::Button("Clear##vertices_set"))
						{
							model_p.selected_vertices_set_->clear();
							mesh_provider_->emit_cells_set_changed(*model_, model_p.selected_vertices_set_);
						}
					}
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("color##vertices", model_p.param_point_sprite_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);

					if (need_update || model_p.object_update_)
					{
						for (View* v : linked_views_)
							v->request_update();
					}

					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					if (ImGui::Button("Accept control set##vertices_set"))
					{
						control_set = model_p.selected_vertices_set_;

						MESH* l_cage = generate_local_cage(*model_, model_p.vertex_position_, control_set);

						model_p.selected_vertices_set_ = nullptr;
					}

					/*imgui_combo_cells_set(md, control_set, "Control ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });*/

					bool newCage = false;
					model_p.object_update_ = false;

					ImGui::Separator();
					ImGui::Text("Binding");

					if (cage_container_.size() > 0)
					{

						MultiToolsDeformation<MESH, GRAPH>* space_deformation =
							static_cast<MultiToolsDeformation<MESH, GRAPH>*>(
								app_.module("MultiToolsDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));

						if (imgui_local_cage_selector(space_deformation, selected_cage_, "Local Cage",
													  [&](MESH& m) { selected_cage_ = &m; }))
						{

							const std::string cage_name = mesh_provider_->mesh_name(*selected_cage_);

							std::string prefix = "local_cage";

							if (cage_name.size() > 0 &&
								std::mismatch(prefix.begin(), prefix.end(), cage_name.begin()).first == prefix.end())
							{

								// inspired from https://github.com/ocornut/imgui/issues/1658
								const char* items[] = {"MVC", "Green"};
								std::string current_item = "MVC";
								ImGuiComboFlags flags = ImGuiComboFlags_NoArrowButton;

								ImGuiStyle& style = ImGui::GetStyle();
								float w = ImGui::CalcItemWidth();
								float spacing = style.ItemInnerSpacing.x;
								float button_sz = ImGui::GetFrameHeight();
								ImGui::PushItemWidth(w - spacing * 2.0f - button_sz * 2.0f);
								if (ImGui::BeginCombo("##custom combo", current_item.c_str(),
													  ImGuiComboFlags_NoArrowButton))
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

								std::shared_ptr<cgogn::modeling::CageDeformationTool<MESH>> cdt =
									cage_container_[cage_name];
								selected_cdt_ = cage_container_[cage_name];
								if (ImGui::Button("Bind object"))
								{
									// ImGui::Text("Binded");
									if (current_item == "MVC")
									{
										bind_influence_cage_mvc(*model_, model_p.vertex_position_, *selected_cage_,
																selected_cdt_->control_cage_vertex_position_);
									}

									else if (current_item == "Green")
									{
										/*if (!cd.local_def)
										{
											bind_object_green(*model_, model_p.vertex_position_,
										*selected_cage_, cd.cage_vertex_position_);
										}
										else
										{
											bind_local_green(*model_, model_p.vertex_position_,
										*selected_cage_, cd.cage_vertex_position_);
										}*/

										std::cout << "green" << std::endl;
									}
									else
									{
										std::cout << "not available yet" << std::endl;
									}
								}

								/*ImGui::Separator();
									float new_smoothing_factor = selected_cdt_->smoothing_factor_;
									ImGui::SliderFloat("Attenuation factor", &new_smoothing_factor, 0.0f, 1.0f);

									if (new_smoothing_factor != selected_cdt_->smoothing_factor_){
										selected_cdt_->smoothing_factor_ = new_smoothing_factor;
										selected_cdt_->update_attenuation(*model_);
									}*/
							}
						}
					}
				}

				if (selected_tool_ == Handle)
				{

					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					imgui_combo_cells_set(model_md, control_set, "Handle ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });

					if (control_set && control_set->size() > 0)
					{
						selected_hdt_ = generate_handle_tool(*model_, model_p.vertex_position_, control_set);
						// handle_name = generate_handle_tool(*model_, model_p.vertex_position_, control_set);
					}

					ImGui::Separator();
					ImGui::Text("Binding");

					CellsSet<MESH, MeshVertex>* influence_set = nullptr;

					imgui_combo_cells_set(model_md, influence_set, "Influence ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { influence_set = cs; });

					// A FIXER
					if (selected_hdt_ && influence_set)
					{

						// std::shared_ptr<cgogn::modeling::HandleDeformationTool<MESH>> current_handle =
						// handle_container_[handle_name];
						bind_handle_influence_area(*model_, model_p.vertex_position_, influence_set,
												   *(selected_hdt_->control_handle_),
												   selected_hdt_->control_handle_vertex_position_);

						influence_set = nullptr;

						surface_selection_->clear_selected_vertices_set(*model_);
					}

					ImGui::Separator();
					/*ImGui::Text("Binding");


					if (handle_container_.size() > 0)
					{
						MultiToolsDeformation<MESH, GRAPH>* space_deformation =
					static_cast<MultiToolsDeformation<MESH, GRAPH>*>( app_.module("MultiToolsDeformation (" +
					std::string{mesh_traits<MESH>::name} + ")")); if (!imgui_handle_selector(space_deformation,
					selected_handle_, "Handle",
												   [&](GRAPH& g) { selected_handle_ = &g; }))
						{
							//std::cout << " there is a bug" << std::endl;
						}
						else
						{

							const std::string handle_name = graph_provider_->mesh_name(*selected_handle_);

							std::string prefix = "local_handle";

							if (handle_name.size() > 0 &&
								 std::mismatch(prefix.begin(), prefix.end(), handle_name.begin()).first ==
					prefix.end())
							{
								const char* items[] = {"MVC", "Green"};
								std::string current_item = "MVC";
								ImGuiComboFlags flags = ImGuiComboFlags_NoArrowButton;

								ImGuiStyle& style = ImGui::GetStyle();
								float w = ImGui::CalcItemWidth();
								float spacing = style.ItemInnerSpacing.x;
								float button_sz = ImGui::GetFrameHeight();
								ImGui::PushItemWidth(w - spacing * 2.0f - button_sz * 2.0f);
								if (ImGui::BeginCombo("##custom combo", current_item.c_str(),
													  ImGuiComboFlags_NoArrowButton))
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

								selected_hdt_ = handle_container_[handle_name];
								if (ImGui::Button("Bind object"))
								{

									if (current_item == "MVC")
									{
										bind_influence_cage_mvc_handle(*model_, model_p.vertex_position_,
																	   *selected_handle_,
																	   selected_hdt_->control_handle_vertex_position_);
									}

									else
									{
										std::cout << "not available yet" << std::endl;
									}
								}
						}


						}
					}*/
				}

				if (selected_tool_ == Axis)
				{
					CellsSet<MESH, MeshVertex>* control_set = nullptr;

					imgui_combo_cells_set(model_md, control_set, "Control ",
										  [&](CellsSet<MESH, MeshVertex>* cs) { control_set = cs; });

					if (control_set && control_set->size() > 0)
					{

						generate_axis(*model_, model_p.vertex_position_, control_set);
					}
				}

				ImGui::Separator();
				ImGui::Separator();
				ImGui::Text("Deform Tool");

				ImGui::RadioButton("Deform Axis", reinterpret_cast<int*>(&deformed_tool_), Axis);
				ImGui::SameLine();
				ImGui::RadioButton("Deform Cage", reinterpret_cast<int*>(&deformed_tool_), Cage);
				ImGui::SameLine();
				ImGui::RadioButton("Deform Handle", reinterpret_cast<int*>(&deformed_tool_), Handle);

				if (deformed_tool_ == Cage)
				{
					if ((cage_container_.size() +  global_cage_container_.size() ) > 0 )
					{
						MultiToolsDeformation<MESH, GRAPH>* multi_tools_deformation =
							static_cast<MultiToolsDeformation<MESH, GRAPH>*>(
								app_.module("MultiToolsDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));

						imgui_cage_selector(multi_tools_deformation, selected_cage_, "Cage",
											[&](MESH& m) { selected_cage_ = &m; });

						if (selected_cage_)
						{

							selected_mesh_ = selected_cage_; 
							Parameters& cage_p = parameters_[selected_mesh_];

							const std::string cage_name = mesh_provider_->mesh_name(*selected_cage_);

							std::string prefix = "local_cage";
							if (cage_name.size() > 0)
							{
								if (std::mismatch(prefix.begin(), prefix.end(), cage_name.begin()).first ==
									prefix.end())
								{
									// local cage
									selected_cdt_ = cage_container_[cage_name];
								}
								else
								{
									selected_gcdt_ = global_cage_container_[cage_name];
								}

								MeshData<MESH>& cage_md = mesh_provider_->mesh_data(*selected_cage_);

								if (ImGui::Button("Choose cage vertices##vertices_set"))
									cage_md.template add_cells_set<MeshVertex>();

								imgui_combo_cells_set(cage_md, cage_p.selected_vertices_set_, "Cage D sets",
													  [&](CellsSet<MESH, MeshVertex>* cs) {
														  cage_p.selected_vertices_set_ = cs;
														  cage_p.update_selected_vertices_vbo();
														  need_update = true;
													  });

								if (cage_p.selected_vertices_set_)
								{
									ImGui::Text("(nb elements: %d)", cage_p.selected_vertices_set_->size());
									if (ImGui::Button("Clear cage selection##vertices_set"))
									{
										cage_p.selected_vertices_set_->clear();
										mesh_provider_->emit_cells_set_changed(*selected_cage_, cage_p.selected_vertices_set_);
									}
								}
								ImGui::TextUnformatted("Drawing parameters");
								need_update |=
									ImGui::ColorEdit3("color##vertices", cage_p.param_point_sprite_->color_.data(),
													  ImGuiColorEditFlags_NoInputs);

								if (need_update || cage_p.object_update_)
								{
									for (View* v : linked_views_)
										v->request_update();
								}

								CellsSet<MESH, MeshVertex>* control_set = nullptr;

								/*if (ImGui::Button("Accept control set##vertices_set"))
								{
									control_set = cage_p.selected_vertices_set_;

									MESH* l_cage =
										generate_local_cage(*selected_mesh_, cage_p.vertex_position_, control_set);

									cage_p.selected_vertices_set_ = nullptr;
								}*/
							}
						}

						/*if (ImGui::Button("Choose cage vertices##vertices_set"))
							md.template add_cells_set<MeshVertex>();
						imgui_combo_cells_set(md, cage_p.selected_vertices_set_, "Cage D sets",
											  [&](CellsSet<MESH, MeshVertex>* cs) {
												  cage_p.selected_vertices_set_ = cs;
												  cage_p.update_selected_vertices_vbo();
												  need_update = true;
											  });

						if (cage_p.selected_vertices_set_)
						{
							ImGui::Text("(nb elements: %d)", cage_p.selected_vertices_set_->size());
							if (ImGui::Button("Clear cage selection##vertices_set"))
							{
								cage_p.selected_vertices_set_->clear();
								mesh_provider_->emit_cells_set_changed(*selected_cage_, cage_p.selected_vertices_set_);
							}
						}
						ImGui::TextUnformatted("Drawing parameters");
						need_update |= ImGui::ColorEdit3("color##vertices", cage_p.param_point_sprite_->color_.data(),
														 ImGuiColorEditFlags_NoInputs);

						if (need_update || cage_p.object_update_)
						{
							for (View* v : linked_views_)
								v->request_update();
						}

						CellsSet<MESH, MeshVertex>* control_set = nullptr;

						if (ImGui::Button("Accept control set##vertices_set"))
						{
							control_set = cage_p.selected_vertices_set_;

							MESH* l_cage = generate_local_cage(*selected_mesh_, cage_p.vertex_position_, control_set);

							cage_p.selected_vertices_set_ = nullptr;
						}*/

						// const std::string cage_name = mesh_provider_->mesh_name(*selected_cage_);
					}
				}

				// mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
			}
		}
	}

private : 
MESH* model_;

MESH* selected_mesh_;
MESH* selected_cage_;
GRAPH* selected_handle_;

std::shared_ptr<cgogn::modeling::CageDeformationTool<MESH>> selected_cdt_;
std::shared_ptr<cgogn::modeling::GlobalCageDeformationTool<MESH>> selected_gcdt_;
std::shared_ptr<cgogn::modeling::HandleDeformationTool<MESH>> selected_hdt_;
std::unordered_map<const MESH*, Parameters> parameters_;

std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;

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
SelectionTool deformed_tool_;

bool influence_cage_mode_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MULTI_TOOLS_DEFORMATION_H_
