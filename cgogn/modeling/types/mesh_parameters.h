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

#ifndef CGOGN_MODELING_TYPES_MESH_PARAMETERS_H_
#define CGOGN_MODELING_TYPES_MESH_PARAMETERS_H_

#include <cgogn/core/types/cells_set.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

#include <cgogn/geometry/algos/picking.h>
#include <cgogn/geometry/algos/selection.h>

#include <cgogn/ui/view.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>

/**
 * @Class parameters 
*/
class Parameters
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;

	using Vec3 = geometry::Vec3;

public:
	Parameters()
		: vertex_position_(nullptr), 
		selection_method_(SelectionMethod::SingleCell),
		  selecting_cell_(SelectingCell::VertexSelect), 
		  selected_vertices_set_(nullptr), vertex_scale_factor_(1.0),
		  sphere_scale_factor_(10.0), object_update_(false), 
		  back_selection_(false), dragging_mesh_(false)
	{
		param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
		param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
		param_point_sprite_->set_vbos({&selected_vertices_vbo_});
	}

	~Parameters()
	{
	}

	/**
	 * update selected vertices vbo 
	 * selected vertices turns into red
	*/
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
				.push_back(value<Vec3>(*mesh_, vertex_position_, v)); });

			rendering::update_vbo(selected_vertices_position, 
									&selected_vertices_vbo_);
		}
	}

	/**
	 * select vertices set 
	 * use picking to raycast on object 
	*/
	void select_vertices_set(ui::View* view, const int32& button, 
							const int32& x, const int32& y)
	{
		if (vertex_position_)
		{
			rendering::GLVec3d near_d = view->unproject(x, y, 0.0);
			rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
			Vec3 A{near_d.x(), near_d.y(), near_d.z()};
			Vec3 B{far_d.x(), far_d.y(), far_d.z()};

			object_update_ = true;

			switch (selection_method_)
			{
				case modeling::SelectionMethod::SingleCell: 
				{
					switch (selecting_cell_)
					{
						case modeling::SelectingCell::VertexSelect:
							if (selected_vertices_set_)
							{
								std::vector<Vertex> picked;
								geometry::picking(*mesh_, vertex_position_.get(), 
												A, B, picked);
 
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

								if (back_selection_)
								{
									CellCache<MESH> cache_back_ = 
									geometry::within_sphere(*mesh_, 
												picked[1], 
									vertex_base_size_ * sphere_scale_factor_, 
									vertex_position_.get());
								selected_depth_vertices_
								.push_back(std::make_pair(picked[0], picked[1]));
								}
							}

						break;
					}
					break;
				}

				case modeling::SelectionMethod::WithinSphere: 
				{
					std::vector<Vertex> picked;
					geometry::picking(*mesh_, vertex_position_.get(), 
											A, B, picked);
					if (!picked.empty())
					{
						if (selection_for_handle_)
						{
							const Vec3& target_position = 
									value<Vec3>(*mesh_, vertex_position_, picked[0]);

							const Vec3& handle_position = 
								value<Vec3>(*mesh_, vertex_position_, handle_mesh_vertex_);
						
							//const float32 target_radius = (target_position - handle_position).norm() * 150.0; 
							//40.0; 
							const float32 target_radius = (target_position - handle_position).norm(); 
							//*1500.0; 

							CellCache<MESH> cache = geometry::within_sphere(
											*mesh_, handle_mesh_vertex_, 
									vertex_base_size_ * target_radius, 
											vertex_position_.get());
							

							switch (selecting_cell_)
							{
								case modeling::SelectingCell::VertexSelect:
									if (selected_vertices_set_)
									{
										switch (button)
										{
										case 0:
											foreach_cell(cache, [&](Vertex v) -> bool {
												selected_vertices_set_->select(v);
											return true;
											});
											break;

										case 1:
											foreach_cell(cache, [&](Vertex v) -> bool {
												selected_vertices_set_->unselect(v);
											return true;
											});
											break;
										}
									}
								break;
							}
		
						} else 
						{
							CellCache<MESH> cache = geometry::within_sphere(
											*mesh_, picked[0], 
									vertex_base_size_ * sphere_scale_factor_, 
											vertex_position_.get());

							switch (selecting_cell_)
							{
								case modeling::SelectingCell::VertexSelect:
									if (selected_vertices_set_)
									{
										switch (button)
										{
											case 0:
											foreach_cell(cache, [&](Vertex v) -> bool {
												selected_vertices_set_->select(v);
											return true;
											});
											break;

											case 1:
											foreach_cell(cache, [&](Vertex v) -> bool {
												selected_vertices_set_->unselect(v);
											return true;
											});
											break;
										}
									}
								break;
							}

							if (back_selection_)
							{
								if (picked.size() > 1)
								{
									CellCache<MESH> cache_back_ = 
										geometry::within_sphere(
										*mesh_, picked[1], 
									vertex_base_size_ * sphere_scale_factor_, 
										vertex_position_.get());

									if (selecting_cell_ == 
										modeling::SelectingCell::VertexSelect)
									{
										if (selected_vertices_set_)
										{
											switch (button)
											{
												case 0:
												foreach_cell(cache_back_,
													 [&](Vertex v) -> bool {
												selected_vertices_set_->select(v);
												return true;
												});
												break;

												case 1:
												foreach_cell(cache_back_, 
														[&](Vertex v) -> bool {
												selected_vertices_set_->unselect(v);
												return true;
												});
												break;
											}
										}
									}
								}
							}
						}
					}
					break;
				}
			}
		}
	}

	/**
	 * d key pressed 
	 * to displace vertices of mesh 
	*/
	void key_pressed_D_event(ui::View* view)
	{
		if (vertex_position_ && selected_vertices_set_ 
				&& selected_vertices_set_->size() > 0)
		{
			drag_z_ = 0.0;
			selected_vertices_set_->foreach_cell([&](Vertex v) {
				const Vec3& pos = 
					value<Vec3>(*mesh_, vertex_position_, v);
				rendering::GLVec4d vec(pos[0], pos[1], pos[2], 1.0);
				vec = 
					view->projection_matrix_d() * 
						view->modelview_matrix_d() * vec;
				vec /= vec[3];
				drag_z_ += (1.0 + vec[2]) / 2.0;
			});
			drag_z_ /= selected_vertices_set_->size();
			previous_drag_pos_ = view->unproject(view->previous_mouse_x(), view->previous_mouse_y(), drag_z_);
			dragging_mesh_ = true;
		}
	}

	/**
	 * release k
	 * reset dragging mesh parameter
	*/
	void key_release_D_event(ui::View* view)
	{
		if (dragging_mesh_)
			dragging_mesh_ = false;
	}

	/**
	 * mouse displacement 
	 * from the previous and current mouse position
	*/
	void mouse_displacement(ui::View* view, const int32& x, const int32& y)
	{
		if (dragging_mesh_)
		{
			rendering::GLVec3d drag_pos = view->unproject(x, y, drag_z_);
			Vec3 t = drag_pos - previous_drag_pos_;
			selected_vertices_set_->foreach_cell([&](Vertex v) { 
				value<Vec3>(*mesh_, vertex_position_, v) += t; });
			previous_drag_pos_ = drag_pos;
		}
	}

	/**
	 * local draw 
	*/
	void local_draw(ui::View* view)
	{
		const rendering::GLMat4& proj_matrix = view->projection_matrix();
		const rendering::GLMat4& view_matrix = view->modelview_matrix();

		if (selecting_cell_ == modeling::SelectingCell::VertexSelect && selected_vertices_set_ &&
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

	MESH* mesh_;
	std::string name_;

	std::shared_ptr<Attribute<Vec3>> vertex_position_;

	std::shared_ptr<boost::synapse::connection> cells_set_connection_;

	bool object_update_;

	bool back_selection_;

	CMap2::Vertex handle_mesh_vertex_; 
	bool selection_for_handle_; 

	std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;

	float32 vertex_scale_factor_;
	float32 vertex_base_size_;
	float32 sphere_scale_factor_;

	rendering::VBO selected_vertices_vbo_;

	std::vector<std::pair<Vertex, Vertex>> selected_depth_vertices_;

	ui::CellsSet<MESH, Vertex>* selected_vertices_set_;

	SelectingCell selecting_cell_;
	SelectionMethod selection_method_;

	bool dragging_mesh_;
	float64 drag_z_;
	rendering::GLVec3d previous_drag_pos_;
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_TYPES_MESH_PARAMETERS_H_
