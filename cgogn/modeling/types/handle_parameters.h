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

#ifndef CGOGN_MODELING_TYPES_HANDLE_PARAMETERS_H_
#define CGOGN_MODELING_TYPES_HANDLE_PARAMETERS_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/vbo_update.h>

#include <cgogn/geometry/algos/picking.h>
#include <cgogn/geometry/algos/selection.h>

#include <cgogn/ui/view.h>
#include <algorithm>

namespace cgogn
{

namespace modeling
{

template <typename GRAPH>

/**
 * @Class graph parameters 
*/
class HandleParameters
{
	template <typename T>
	using Attribute = typename mesh_traits<GRAPH>::template Attribute<T>;
	using Vertex = typename mesh_traits<GRAPH>::Vertex;
	using Edge = typename mesh_traits<GRAPH>::Edge;

	using Vec3 = geometry::Vec3;

public:

	HandleParameters()
		: vertex_position_(nullptr), vertex_scale_factor_(1.0), 
			sphere_scale_factor_(10.0), selected_vertices_set_(nullptr), 
			object_update_(false), 
			selecting_cell_(SelectingCell::VertexSelect),
		  selection_method_(SelectionMethod::SingleCell), 
		  dragging_handle_(false)
	{
		param_point_sprite_ = 
			rendering::ShaderPointSprite::generate_param();
		param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
		param_point_sprite_->set_vbos({&selected_vertices_vbo_});
	}

	~HandleParameters()
	{
	}

	
	HandleParameters(HandleParameters&&) = default;
	HandleParameters& operator=(HandleParameters&&) = default;


	/// @brief update selected vertices 
	/// update vbo 
	void update_selected_vertices_vbo()
	{
		if (selected_vertices_set_)
		{
			std::vector<Vec3> selected_vertices_position;
			selected_vertices_position
				.reserve(selected_vertices_set_->size());

			selected_vertices_set_->foreach_cell(
				[&](Vertex v) { 
				selected_vertices_position
						.push_back(value<Vec3>(*graph_, vertex_position_, v)); 
				});
			rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);
		}
	}

	/**
	 * select vertices set 
	 * use picking sphere to raycast on a sphere
	*/

	/// @brief select vertices set 
	/// call picking sphere to raycast on a sphere 
	/// @param view current view
	/// @param button select (left) or unselect (right)
	/// @param x retrieved from mouse position
	/// @param y retrieved from mouse position
	void select_vertices_set(ui::View* view, const int32& button, 
							const int32& x, const int32& y)
	{
		if (selected_vertices_set_)
		{
			rendering::GLVec3d near_d = view->unproject(x, y, 0.0);
			rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
			Vec3 A{near_d.x(), near_d.y(), near_d.z()};
			Vec3 B{far_d.x(), far_d.y(), far_d.z()};

			std::vector<Vertex> picked;

			geometry::picking_sphere(*graph_, 
								vertex_position_.get(), 5.0, A, B, picked);
			if (!picked.empty())
			{
				switch (button)
				{
				case 0:
					selected_vertices_set_->select(picked[0]);
					break;
				case 1:
					selected_vertices_set_->unselect(picked[0]);
					break;
				}
			}
		}
	}
	
	/// @brief H key pressed
	/// initialize displacement of handle
	/// @param view current view 
	void key_pressed_H_event(ui::View* view)
	{
		if (vertex_position_ && selected_vertices_set_ && 
								selected_vertices_set_->size() > 0)
		{
			drag_z_ = 0.0;
			selected_vertices_set_->foreach_cell([&](Vertex v) {
				const Vec3& pos = value<Vec3>(*graph_, vertex_position_, v);
				rendering::GLVec4d vec(pos[0], pos[1], pos[2], 1.0);

				vec = 
					view->projection_matrix_d() * 
						view->modelview_matrix_d() * vec;
				vec /= vec[3];
				drag_z_ += (1.0 + vec[2]) / 2.0;
			});
			drag_z_ /= selected_vertices_set_->size();
			previous_drag_pos_ = 
				view->unproject(view->previous_mouse_x(), 
								view->previous_mouse_y(), drag_z_);
			dragging_handle_ = true;
		}
	}


	/// @brief h key released 
	/// reset dragging handle parameter
	/// @param view 
	void key_release_handle_event(ui::View* view)
	{
		if (dragging_handle_)
			dragging_handle_ = false;
	}


	/// @brief mouse displacement outcomes 
	/// call for the corresponding displacement functions 
	/// @param view current view
	/// @param x retrieved from mouse position
	/// @param y retrieved from mouse position
	void mouse_displacement(ui::View* view, const int32& x, const int32& y)
	{
		if (dragging_handle_)
		{
			mouse_handle_displacement(view, x, y); 
		}
	}

	
	/// @brief update rendering of selected vertices 
	/// @param view current view 
	void local_draw(ui::View* view)
	{
		const rendering::GLMat4& proj_matrix = view->projection_matrix();
		const rendering::GLMat4& view_matrix = view->modelview_matrix();

		if (selecting_cell_ == modeling::SelectingCell::VertexSelect && 
								selected_vertices_set_ &&
							selected_vertices_set_->size() > 0 && 
				param_point_sprite_->attributes_initialized())
		{
			param_point_sprite_->point_size_ = 
				vertex_base_size_ * vertex_scale_factor_;
			param_point_sprite_->bind(proj_matrix, view_matrix);
			glDrawArrays(GL_POINTS, 0, selected_vertices_set_->size());
			param_point_sprite_->release();
		}
	}

	GRAPH* graph_;
	std::string name_;

	std::shared_ptr<Attribute<Vec3>> vertex_position_;

	std::unique_ptr<rendering::ShaderPointSprite::Param> 
													param_point_sprite_;
	std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;

	float32 vertex_scale_factor_;
	float32 vertex_base_size_;
	float32 sphere_scale_factor_;

	rendering::VBO selected_vertices_vbo_;

	bool object_update_;

	std::vector<std::pair<Vertex, Vertex>> selected_depth_vertices_;

	ui::CellsSet<GRAPH, Vertex>* selected_vertices_set_;

	SelectingCell selecting_cell_;
	SelectionMethod selection_method_;

	std::shared_ptr<boost::synapse::connection> cells_set_connection_;

	Vec3 normal_;

	float64 drag_z_;
	rendering::GLVec3d previous_drag_pos_;

private: 
	bool dragging_handle_;

	/// @brief handle displacement caused by the mouse displacement
	/// handle translates along normal 
	/// @param view current view 
	/// @param x retrieved from mouse position
	/// @param y retrieved from mouse position
	void mouse_handle_displacement(ui::View* view, const int32& x, const int32& y)
	{
		rendering::GLVec3d drag_pos = view->unproject(x, y, drag_z_);
			Vec3 t = drag_pos - previous_drag_pos_;

			Vec3 t_bis = t.dot(normal_) * normal_;

			selected_vertices_set_->foreach_cell([&](Vertex v) { 
				value<Vec3>(*graph_, vertex_position_, v) += t_bis; });
			
			previous_drag_pos_ = drag_pos + t_bis;
	}



};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_TYPES_GRAPH_PARAMETERS_H_
