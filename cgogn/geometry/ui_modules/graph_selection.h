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

#ifndef CGOGN_MODULE_GRAPH_SELECTION_H_
#define CGOGN_MODULE_GRAPH_SELECTION_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/graph_provider.h>

#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/picking.h>
#include <cgogn/geometry/algos/selection.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/vbo_update.h>

#include <boost/synapse/connect.hpp>

#include <unordered_map>

#undef near
#undef far

namespace cgogn
{

namespace ui
{

template <typename GRAPH>
class GraphSelection : public ViewModule
{
	//static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceSelectionPO can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<GRAPH>::template Attribute<T>;

	using Vertex = typename mesh_traits<GRAPH>::Vertex;
	using Edge = typename mesh_traits<GRAPH>::Edge;

	enum SelectingCell
	{
		VertexSelect = 0,
	};

	enum SelectionMethod
	{
		SingleCell = 0,
	};

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), vertex_scale_factor_(1.0), sphere_scale_factor_(10.0),
			  selected_vertices_set_(nullptr),
			  selecting_cell_(VertexSelect), selection_method_(SingleCell)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_point_sprite_->set_vbos({&selected_vertices_vbo_});
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

	public:
		void update_selected_vertices_vbo()
		{
			if (selected_vertices_set_)
			{
				std::vector<Vec3> selected_vertices_position;
				selected_vertices_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell(
					[&](Vertex v) { selected_vertices_position.push_back(value<Vec3>(*graph_, vertex_position_, v)); });
				rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);
			}
		}

		GRAPH* graph_;
		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;

		rendering::VBO selected_vertices_vbo_;

		CellsSet<GRAPH, Vertex>* selected_vertices_set_;
		
		SelectingCell selecting_cell_;
		SelectionMethod selection_method_;
	};

public:
	GraphSelection(const App& app)
		: ViewModule(app, "GraphSelection (" + std::string{mesh_traits<GRAPH>::name} + ")"), selected_graph_(nullptr)
	{
	}

	~GraphSelection()
	{
	}

private:
	void init_graph(GRAPH* g)
	{
		Parameters& p = parameters_[g];
		p.graph_ = g;
		graph_connections_[g].push_back(
			boost::synapse::connect<typename GraphProvider<GRAPH>::template attribute_changed_t<Vec3>>(
				g, [this, g](Attribute<Vec3>* attribute) {
					Parameters& p = parameters_[g];
					if (p.vertex_position_.get() == attribute)
					{
						p.vertex_base_size_ = 5.0; //0.3; //float32(geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6);
						p.update_selected_vertices_vbo();
					}

					for (View* v : linked_views_)
						v->request_update();
				}));
		graph_connections_[g].push_back(
			boost::synapse::connect<typename GraphProvider<GRAPH>::template cells_set_changed<Vertex>>(
				g, [this, g](CellsSet<GRAPH, Vertex>* set) {
					Parameters& p = parameters_[g];
					if (p.selected_vertices_set_ == set && p.vertex_position_)
					{
						p.update_selected_vertices_vbo();
						for (View* v : linked_views_)
							v->request_update();
					}
				}));
	}

public:
	void set_vertex_position(const GRAPH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];

		p.vertex_position_ = vertex_position;
		if (p.vertex_position_)
		{
			p.vertex_base_size_ = 5.0; //0.3; //float32(geometry::mean_edge_length(m, p.vertex_position_.get()) / 6); // 6 ???
			p.update_selected_vertices_vbo(); 
		}

		for (View* v : linked_views_)
			v->request_update();
	}

	void set_selected_graph(GRAPH& m)
	{
		selected_graph_ = &m;
	}

	template <typename CELL>
	void set_selected_cells_set(const GRAPH& m, CellsSet<GRAPH, CELL>* cs)
	{
		Parameters& p = parameters_[&m];
		if constexpr (std::is_same_v<CELL, Vertex>)
		{
			p.selected_vertices_set_ = cs;
			p.update_selected_vertices_vbo();
		}

	}

protected:
	void init() override
	{
		graph_provider_ = static_cast<ui::GraphProvider<GRAPH>*>(
			app_.module("GraphProvider (" + std::string{mesh_traits<GRAPH>::name} + ")"));
		graph_provider_->foreach_graph([this](GRAPH& m, const std::string&) { init_graph(&m); });
		connections_.push_back(boost::synapse::connect<typename GraphProvider<GRAPH>::graph_added>(
			graph_provider_, this, &GraphSelection<GRAPH>::init_graph));
	}

	void mouse_press_event(View* view, int32 button, int32 x, int32 y) override
	{
		if (selected_graph_ && view->control_pressed())
		{
			Parameters& p = parameters_[selected_graph_];

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
							std::vector<Vertex> picked;
							cgogn::geometry::picking_sphere(*selected_graph_, p.vertex_position_.get(), 50, A, B, picked);
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
								graph_provider_->emit_cells_set_changed(*selected_graph_, p.selected_vertices_set_);
							} 
						}
						break;
					
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
		bool need_update = false;

		imgui_graph_selector(graph_provider_, selected_graph_, "Surface", [&](GRAPH& m) {
			selected_graph_ = &m;
			graph_provider_->graph_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_graph_)
		{
			// float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;
			Parameters& p = parameters_[selected_graph_];

			imgui_combo_attribute<Vertex, Vec3>(*selected_graph_, p.vertex_position_, "Position",
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													set_vertex_position(*selected_graph_, attribute);
												});

			if (p.vertex_position_)
			{
				ImGui::Separator();
				int* ptr_sel_cell = reinterpret_cast<int*>(&p.selecting_cell_);
				need_update |= ImGui::RadioButton("Vertex", ptr_sel_cell, VertexSelect);

				ImGui::RadioButton("Single", reinterpret_cast<int*>(&p.selection_method_), SingleCell);

				GraphData<GRAPH>& md = graph_provider_->graph_data(*selected_graph_);

				if (p.selecting_cell_ == VertexSelect)
				{
					if (ImGui::Button("Create set##vertices_set"))
						md.template add_cells_set<Vertex>();
					imgui_combo_cells_set_graph(md, p.selected_vertices_set_, "Sets", [&](CellsSet<GRAPH, Vertex>* cs) {
						p.selected_vertices_set_ = cs;
						p.update_selected_vertices_vbo();
						need_update = true;
					});
					if (p.selected_vertices_set_)
					{
						ImGui::Text("(nb elements: %d)", p.selected_vertices_set_->size());
						if (ImGui::Button("Clear##vertices_set"))
						{
							p.selected_vertices_set_->clear();
							graph_provider_->emit_cells_set_changed(*selected_graph_, p.selected_vertices_set_);
						}
					}
					ImGui::TextUnformatted("Drawing parameters");
					need_update |= ImGui::ColorEdit3("color##vertices", p.param_point_sprite_->color_.data(),
													 ImGuiColorEditFlags_NoInputs);
					need_update |= ImGui::SliderFloat("size##vertices", &(p.vertex_scale_factor_), 0.1f, 2.0f);
				}
			}
		}

		if (need_update)
			for (View* v : linked_views_)
				v->request_update();
	}

private:
	const GRAPH* selected_graph_;
	std::unordered_map<const GRAPH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	std::unordered_map<const GRAPH*, std::vector<std::shared_ptr<boost::synapse::connection>>> graph_connections_;
	GraphProvider<GRAPH>* graph_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_GRAPH_SELECTION_H_
