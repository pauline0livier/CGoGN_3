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

#include <cgogn/core/ui_modules/graph_provider.h>
#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/rendering/types.h>
#include <cgogn/rendering/ui_modules/graph_render.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/angle.h>

#include <GLFW/glfw3.h>

#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

#include <cgogn/modeling/types/axis_deformation_tool.h>
#include <cgogn/modeling/types/cage_deformation_tool.h>
#include <cgogn/modeling/types/global_cage_deformation_tool.h>
#include <cgogn/modeling/types/graph_parameters.h>
#include <cgogn/modeling/types/handle_deformation_tool.h>
#include <cgogn/modeling/types/mesh_parameters.h>

#include <boost/synapse/connect.hpp>

#include <iostream>

#include <algorithm>
#include <string>

namespace cgogn
{

namespace ui
{

enum SelectionTool
{
	Handle = 0,
	Axis,
	Cage
};

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

public:
	std::unordered_map<std::string, std::shared_ptr<modeling::HandleDeformationTool<MESH>>> handle_container_;

	MultiToolsDeformation(const App& app)
		: ViewModule(app, "MultiToolsDeformation (" + std::string{mesh_traits<MESH>::name} + ")"), model_(nullptr),
		  influence_set_(nullptr), selected_mesh_(nullptr), selected_graph_(nullptr), selected_handle_(nullptr),
		  new_tool_(false)//, //dragging_mesh_(false), dragging_handle_(false)
	{
	}

	~MultiToolsDeformation()
	{
	}

	void init_mesh(MESH* m)
	{
		
		parameters_[m] = new modeling::Parameters<MESH>();  
		modeling::Parameters<MESH>& p = *parameters_[m]; 
		p.mesh_ = m;

		mesh_connections_[m].push_back(
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				m, [this, m](MeshAttribute<Vec3>* attribute) {
					modeling::Parameters<MESH>& p = *parameters_[m];
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
					modeling::Parameters<MESH>& p = *parameters_[m];
					if (p.selected_vertices_set_ == set && p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));

		p.cells_set_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<MeshVertex>>(
				m, [this, m](CellsSet<MESH, MeshVertex>* set) {
					modeling::Parameters<MESH>& p = *parameters_[m];
					if (p.selected_vertices_set_ == set)
						std::cout << "ok " << std::endl;
					// p.solver_ready_ = false;
				});

		MeshData<MESH>& md = mesh_provider_->mesh_data(*m);
		md.template add_cells_set<MeshVertex>();
	}

	void init_graph(GRAPH* g)
	{
		graph_parameters_[g] = new modeling::GraphParameters<GRAPH>(); 
		modeling::GraphParameters<GRAPH>& p = *graph_parameters_[g];
		p.graph_ = g;

		graph_connections_[g].push_back(
			boost::synapse::connect<typename GraphProvider<GRAPH>::template attribute_changed_t<Vec3>>(
				g, [this, g](GraphAttribute<Vec3>* attribute) {
					modeling::GraphParameters<GRAPH>& p = *graph_parameters_[g];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = 6.0;
						p.update_selected_vertices_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));
		graph_connections_[g].push_back(
			boost::synapse::connect<typename GraphProvider<GRAPH>::template cells_set_changed<GraphVertex>>(
				g, [this, g](CellsSet<GRAPH, GraphVertex>* set) {
					modeling::GraphParameters<GRAPH>& p = *graph_parameters_[g];
					if (p.selected_vertices_set_ == set && p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));

		/*p.cells_set_connection_ =
			boost::synapse::connect<typename GraphProvider<GRAPH>::template cells_set_changed<GraphVertex>>(
				g, [this, g](CellsSet<GRAPH, GraphVertex>* set) {
					modeling::GraphParameters<GRAPH>& p = graph_parameters_[g];
					if (p.selected_vertices_set_ == set)
						std::cout << "ok " << std::endl;
					// p.solver_ready_ = false;
				});*/
		GraphData<GRAPH>& gd = graph_provider_->graph_data(*g);
		gd.template add_cells_set<GraphVertex>();
	}

	void set_model(MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position,
				   const std::shared_ptr<MeshAttribute<Vec3>>& vertex_normal)
	{
		model_ = &m;

		geometry::compute_normal<MeshVertex>(m, vertex_position.get(), vertex_normal.get());
		mesh_provider_->emit_attribute_changed(m, vertex_normal.get());

		set_vertex_position(m, vertex_position);
	}

	template <typename FUNC>
	void foreach_handle(const FUNC& f)
	{
		for (auto& [name, hdt] : handle_container_)
			f(*(hdt->control_handle_), name);
	}

	/////////////
	// SIGNALS //
	/////////////
	using handle_added = struct handle_added_ (*)(modeling::HandleDeformationTool<MESH>* h);

private:
	void set_vertex_position(MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position)
	{
		modeling::Parameters<MESH>& p = *parameters_[&m];
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

	void set_graph_vertex_position(const GRAPH& g, const std::shared_ptr<GraphAttribute<Vec3>>& vertex_position)
	{
		modeling::GraphParameters<GRAPH>& p = *graph_parameters_[&g];
		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = 5.0;
			p.update_selected_vertices_vbo();
		}

		for (View* v : linked_views_)
			v->request_update();
	}

	//
	// const std::string
	std::shared_ptr<modeling::HandleDeformationTool<MESH>> generate_handle_tool(
		MESH& m, const std::shared_ptr<MeshAttribute<Vec3>>& vertex_position, CellsSet<MESH, MeshVertex>* handle_set)
	{
		int handle_number = handle_container_.size();
		const std::string handle_name = "local_handle" + std::to_string(handle_number);
		GRAPH* handle = graph_provider_->add_graph(handle_name);

		auto handle_vertex_position = add_attribute<Vec3, GraphVertex>(*handle, "position");
		auto handle_vertex_radius = add_attribute<Scalar, GraphVertex>(*handle, "radius");

		set_graph_vertex_position(*handle, handle_vertex_position);

		auto mesh_vertex_normal = get_attribute<Vec3, MeshVertex>(m, "normal");

		const Vec3 handle_position = modeling::get_mean_value_attribute_from_set(m, vertex_position.get(), handle_set);

		Vec3 normal = modeling::get_mean_value_attribute_from_set(m, mesh_vertex_normal.get(), handle_set);
		normal.normalize();

		CMap2::Vertex closest_vertex =
			modeling::closest_vertex_in_set_from_value(m, vertex_position.get(), handle_set, handle_position);

		const auto [it, inserted] =
			handle_container_.emplace(handle_name, std::make_shared<modeling::HandleDeformationTool<MESH>>());
		modeling::HandleDeformationTool<MESH>* hdt = it->second.get();

		if (inserted)
		{

			modeling::Parameters<MESH>& p = *parameters_[&m];
			modeling::GraphParameters<GRAPH>& handle_p = *graph_parameters_[handle];

			hdt->create_space_tool(handle, handle_vertex_position.get(), handle_vertex_radius.get(), handle_position,
								   normal);

			hdt->set_handle_mesh_vertex(closest_vertex);

			handle_p.normal_ = normal;

			graph_provider_->emit_connectivity_changed(*handle);
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_position.get());
			graph_provider_->emit_attribute_changed(*handle, handle_vertex_radius.get());

			View* v1 = app_.current_view();
			graph_render_->set_vertex_position(*v1, *handle, handle_vertex_position);

			graph_render_->set_vertex_radius(*v1, *handle, handle_vertex_radius);

			auto object_geodesic = get_or_add_attribute<Scalar, MeshVertex>(*model_, "geodesic_distance");

			hdt->set_geodesic_distance(m, vertex_position);
			mesh_provider_->emit_attribute_changed(m, object_geodesic.get());

			boost::synapse::emit<handle_added>(this, hdt);
		}

		return handle_container_[handle_name];
		// return handle_name;
	}

	/// BIND

	void bind_local_handle(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position,
						   GRAPH& control_handle, const std::shared_ptr<GraphAttribute<Vec3>>& handle_vertex_position,
						   std::string binding_type)
	{

		modeling::Parameters<MESH>& p = *parameters_[&object];
		modeling::GraphParameters<GRAPH>& p_handle = *graph_parameters_[&control_handle];

		std::shared_ptr<modeling::HandleDeformationTool<MESH>> hdt = handle_container_[p_handle.name_];

		hdt->set_deformation_type(binding_type);

		MeshData<MESH>& md = mesh_provider_->mesh_data(object);
		CellsSet<MESH, MeshVertex>& i_set = md.template add_cells_set<MeshVertex>();

		hdt->influence_area_ = &i_set;

		hdt->set_influence_area(object, object_vertex_position, influence_set_);

		if (binding_type == "Spike")
		{
			hdt->set_attenuation_spike(object, object_vertex_position);
		}

		if (binding_type == "Round")
		{
			hdt->set_attenuation_round(object, object_vertex_position);
		}

		/*if (global_cage_container_.size() > 0)
		{
			for (auto& [name, gcdt] : global_cage_container_)
			{
				GraphVertex handle_vertex = hdt->get_handle_vertex();

				Eigen::VectorXf weights =
					gcdt->bind_mvc_handle(*(hdt->control_handle_), p_handle.vertex_position_, handle_vertex);

				hdt->global_cage_weights_ = weights;
			}
		}*/

		mesh_provider_->emit_cells_set_changed(object, hdt->influence_area_);

		hdt->handle_attribute_update_connection_ =
			boost::synapse::connect<typename GraphProvider<GRAPH>::template attribute_changed_t<Vec3>>(
				&control_handle, [&](GraphAttribute<Vec3>* attribute) {
					if (handle_vertex_position.get() == attribute)
					{
						if (deformed_tool_ == Handle)
						{
							std::shared_ptr<modeling::HandleDeformationTool<MESH>> current_hdt =
								handle_container_[p_handle.name_];

							// current_hdt-> update_deformation_object(object, object_vertex_position);
							std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
								get_attribute<uint32, MeshVertex>(object, "vertex_index");

							const Vec3 new_deformation = current_hdt->get_handle_deformation();

							/// graph_provider_->emit_attribute_changed(
							//*(current_hdt->control_handle_), current_hdt->control_handle_vertex_position_.get());

							current_hdt->influence_area_->foreach_cell([&](MeshVertex v) -> bool {
								uint32 vidx = value<uint32>(object, object_vertex_index, v);

								value<Vec3>(object, object_vertex_position, v) +=
									current_hdt->attenuation_[vidx] * new_deformation;

								return true;
							});

							mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
						}
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
			mesh_provider_, this, &MultiToolsDeformation<MESH, GRAPH>::init_mesh));

		graph_provider_ = static_cast<ui::GraphProvider<GRAPH>*>(
			app_.module("GraphProvider (" + std::string{mesh_traits<GRAPH>::name} + ")"));
		graph_provider_->foreach_graph([this](GRAPH& g, const std::string&) { init_graph(&g); });
		connections_.push_back(boost::synapse::connect<typename GraphProvider<GRAPH>::graph_added>(
			graph_provider_, this, &MultiToolsDeformation<MESH, GRAPH>::init_graph));

		surface_render_ = static_cast<ui::SurfaceRender<MESH>*>(
			app_.module("SurfaceRender (" + std::string{mesh_traits<MESH>::name} + ")"));

		graph_render_ = static_cast<ui::GraphRender<GRAPH>*>(
			app_.module("GraphRender (" + std::string{mesh_traits<GRAPH>::name} + ")"));
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		if (selected_mesh_ && view->shift_pressed())
		{

			modeling::Parameters<MESH>& p = *parameters_[selected_mesh_];

			p.select_vertices_set(view, button, x, y); 
			mesh_provider_->emit_cells_set_changed(*selected_mesh_,
													p.selected_vertices_set_);
		}

		if (selected_graph_ && view->control_pressed())
		{
			modeling::GraphParameters<GRAPH>& p = *graph_parameters_[selected_graph_];

			if (p.vertex_position_)
			{
				p.select_vertices_set(view, button, x, y); 
				graph_provider_->emit_cells_set_changed(*selected_graph_, p.selected_vertices_set_);
			}
		}
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_D)
		{
			if (selected_mesh_)
			{
				modeling::Parameters<MESH>& p = *parameters_[selected_mesh_];
				p.key_pressed_D_event(view); 
			}
		}

		if (key_code == GLFW_KEY_H)
		{ // handle
			if (selected_graph_)
			{
				modeling::GraphParameters<GRAPH>& p = *graph_parameters_[selected_graph_];
				p.key_pressed_H_event(view); 
			}
		}
	}

	void key_release_event(View* view, int32 key_code) override
	{
		unused_parameters(view);
		if (key_code == GLFW_KEY_D)
		{
			modeling::Parameters<MESH>& p = *parameters_[selected_mesh_];
			p.key_release_D_event(view); 
		}
		else if (key_code == GLFW_KEY_H || key_code == GLFW_KEY_Q || key_code == GLFW_KEY_A)
		{
			modeling::GraphParameters<GRAPH>& p = *graph_parameters_[selected_graph_];
			p.key_release_graph_event(view); 
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y) override
	{

		if (selected_mesh_)
			{
				modeling::Parameters<MESH>& p = *parameters_[selected_mesh_];
				p.mouse_displacement(view, x, y); 

				mesh_provider_->emit_attribute_changed(*selected_mesh_, p.vertex_position_.get());
			}

		if (selected_graph_)
			{
				modeling::GraphParameters<GRAPH>& p = *graph_parameters_[selected_graph_];
				p.mouse_displacement(view, x, y); 
				graph_provider_->emit_attribute_changed(*selected_graph_, p.vertex_position_.get());
			}
	}

	void draw(View* view) override
	{
		for (auto& [m, p] : parameters_)
		{
			p->local_draw(view); 
		}

		for (auto& [m, p] : graph_parameters_)
		{
			p->local_draw(view);
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

			modeling::Parameters<MESH>& model_p = *parameters_[model_];

			auto object_vertex_position = get_attribute<Vec3, MeshVertex>(*model_, "position");

			set_vertex_position(*model_, object_vertex_position);

			if (model_p.vertex_position_)
			{
				ImGui::Separator();
				ImGui::Separator();
				ImGui::Text("Create tool");
				ImGui::Separator();

				ImGui::Text("Local");
				ImGui::RadioButton("Handle", reinterpret_cast<int*>(&selected_tool_), Handle);

				if (selected_tool_ == Handle)
				{
					selected_mesh_ = model_;
					ImGui::Separator();

					imgui_combo_cells_set(model_md, model_p.selected_vertices_set_, "Set_for_handle",
										  [&](CellsSet<MESH, MeshVertex>* cs) {
											  model_p.selected_vertices_set_ = cs;
											  model_p.update_selected_vertices_vbo();
											  need_update = true;
										  });

					ImGui::RadioButton("Set_handle", reinterpret_cast<int*>(&model_p.selection_method_),
									   (int) modeling::SelectionMethod::SingleCell);
					ImGui::SameLine();
					ImGui::RadioButton("Set_influence_area", reinterpret_cast<int*>(&model_p.selection_method_),
									   (int) modeling::SelectionMethod::WithinSphere);

					if (model_p.selection_method_ == modeling::SelectionMethod::SingleCell)
					{
						if (model_p.selected_vertices_set_)
						{

							if (model_p.selected_vertices_set_->size() > 0)
							{
								CellsSet<MESH, MeshVertex>* control_set = model_p.selected_vertices_set_;

								generate_handle_tool(*model_, model_p.vertex_position_, control_set);

								new_tool_ = true;

								model_p.selected_vertices_set_->clear();
								mesh_provider_->emit_cells_set_changed(*model_, model_p.selected_vertices_set_);
							}
						}
					}

					if (model_p.selection_method_ == modeling::SelectionMethod::WithinSphere)
					{

						ImGui::SliderFloat("Sphere_radius", &(model_p.sphere_scale_factor_), 10.0f, 100.0f);

						if (model_p.selected_vertices_set_)
						{
							ImGui::Text("(nb elements: %d)", model_p.selected_vertices_set_->size());
							if (ImGui::Button("Clear handle area ##vertices_set"))
							{
								model_p.selected_vertices_set_->clear();
								mesh_provider_->emit_cells_set_changed(*model_, model_p.selected_vertices_set_);
							}
						}
						ImGui::TextUnformatted("Drawing parameters");
						need_update |=
							ImGui::ColorEdit3("color_handle##vertices", model_p.param_point_sprite_->color_.data(),
											  ImGuiColorEditFlags_NoInputs);

						if (need_update || model_p.object_update_)
						{
							for (View* v : linked_views_)
								v->request_update();
						}

						if (ImGui::Button("Accept handle influence area##vertices_set"))
						{
							influence_set_ = model_p.selected_vertices_set_;

							model_p.selected_vertices_set_ = nullptr;
						}

						model_p.object_update_ = false;
					}

					ImGui::Separator();
					ImGui::Text("Binding");

					if (new_tool_ && influence_set_)
					{

						MultiToolsDeformation<MESH, GRAPH>* space_deformation =
							static_cast<MultiToolsDeformation<MESH, GRAPH>*>(
								app_.module("MultiToolsDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));

						if (imgui_handle_selector(space_deformation, selected_handle_, "Local Handle",
												  [&](GRAPH& g) { selected_handle_ = &g; }))
						{

							const std::string handle_name = graph_provider_->graph_name(*selected_handle_);

							std::string prefix = "local_handle";

							if (handle_name.size() > 0 &&
								std::mismatch(prefix.begin(), prefix.end(), handle_name.begin()).first == prefix.end())
							{
								const char* items[] = {"Spike", "Round"};
								std::string current_item = "Spike";
								ImGuiComboFlags flags = ImGuiComboFlags_NoArrowButton;

								ImGuiStyle& style = ImGui::GetStyle();
								float w = ImGui::CalcItemWidth();
								float spacing = style.ItemInnerSpacing.x;
								float button_sz = ImGui::GetFrameHeight();
								ImGui::PushItemWidth(w - spacing * 2.0f - button_sz * 2.0f);
								if (ImGui::BeginCombo("##custom_combo", current_item.c_str(),
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

								if (ImGui::Button("Bind object"))
								{

									std::shared_ptr<modeling::HandleDeformationTool<MESH>> hdt =
										handle_container_[handle_name];

									modeling::GraphParameters<GRAPH>& local_p = *graph_parameters_[selected_handle_];

									local_p.name_ = handle_name;

									bind_local_handle(*model_, model_p.vertex_position_, *(hdt->control_handle_),
													  hdt->control_handle_vertex_position_, current_item);

									influence_set_ = nullptr;
									new_tool_ = false;
								}
							}
						}
					}
				}

				ImGui::Separator();
				ImGui::Separator();
				ImGui::Text("Deform Tool");

				ImGui::RadioButton("Deform Handle", reinterpret_cast<int*>(&deformed_tool_), Handle);

				if (deformed_tool_ == Handle)
				{

					if (handle_container_.size() > 0)
					{
						MultiToolsDeformation<MESH, GRAPH>* multi_tools_deformation =
							static_cast<MultiToolsDeformation<MESH, GRAPH>*>(
								app_.module("MultiToolsDeformation (" + std::string{mesh_traits<MESH>::name} + ")"));

						imgui_handle_selector(multi_tools_deformation, selected_handle_, "Handle",
											  [&](GRAPH& g) { selected_handle_ = &g; });

						if (selected_handle_)
						{

							selected_graph_ = selected_handle_;
							modeling::GraphParameters<GRAPH>& handle_p = *graph_parameters_[selected_graph_];

							const std::string handle_name = graph_provider_->graph_name(*selected_handle_);

							std::string prefix = "local_handle";
							if (handle_name.size() > 0)
							{
								selected_hdt_ = handle_container_[handle_name];

								GraphData<GRAPH>& handle_gd = graph_provider_->graph_data(*selected_handle_);

								//if (ImGui::Button("Choose handle vertices##vertices_set"))

									//handle_gd.template add_cells_set<GraphVertex>();

								imgui_combo_cells_set_graph(handle_gd, handle_p.selected_vertices_set_, "Handle_sets",
															[&](CellsSet<GRAPH, GraphVertex>* cs) {
																handle_p.selected_vertices_set_ = cs;
																handle_p.update_selected_vertices_vbo();
																need_update = true;
															});

								if (handle_p.selected_vertices_set_)
								{
									ImGui::Text("(nb elements: %d)", handle_p.selected_vertices_set_->size());
									if (ImGui::Button("Clear handle selection##vertices_set"))
									{
										handle_p.selected_vertices_set_->clear();
										graph_provider_->emit_cells_set_changed(*selected_handle_,
																				handle_p.selected_vertices_set_);
									}
								}
								ImGui::TextUnformatted("Drawing parameters");
								need_update |=
									ImGui::ColorEdit3("color##vertices", handle_p.param_point_sprite_->color_.data(),
													  ImGuiColorEditFlags_NoInputs);

								if (need_update || handle_p.object_update_)
								{
									for (View* v : linked_views_)
										v->request_update();

									/*std::shared_ptr<modeling::HandleDeformationTool<MESH>> current_hdt =
										handle_container_[handle_p.name_];

									graph_provider_->emit_attribute_changed(
										*(current_hdt->control_handle_),
										current_hdt->control_handle_vertex_position_.get());*/
								}
							}
						}
					}
				}

				// mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
			}
		}
	}

private:
	MESH* model_;

	MESH* selected_mesh_;
	GRAPH* selected_handle_;
	GRAPH* selected_graph_;

	std::shared_ptr<modeling::CageDeformationTool<MESH>> selected_cdt_;
	std::shared_ptr<modeling::GlobalCageDeformationTool<MESH>> selected_gcdt_;
	std::shared_ptr<modeling::HandleDeformationTool<MESH>> selected_hdt_;
	std::shared_ptr<modeling::AxisDeformationTool<MESH>> selected_adt_;

	std::unordered_map<const MESH*, modeling::Parameters<MESH>*> parameters_;
	std::unordered_map<const GRAPH*, modeling::GraphParameters<GRAPH>*> graph_parameters_;

	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;

	std::unordered_map<const GRAPH*, std::vector<std::shared_ptr<boost::synapse::connection>>> graph_connections_;

	MeshProvider<MESH>* mesh_provider_;
	GraphProvider<GRAPH>* graph_provider_;

	SurfaceRender<MESH>* surface_render_;
	GraphRender<GRAPH>* graph_render_;

	SelectionTool selected_tool_;
	SelectionTool deformed_tool_;

	bool axis_init_;

	CellsSet<MESH, MeshVertex>* influence_set_;
	bool new_tool_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MULTI_TOOLS_DEFORMATION_H_
