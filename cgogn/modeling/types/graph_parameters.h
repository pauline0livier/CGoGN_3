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

#ifndef CGOGN_MODELING_TYPES_GRAPH_PARAMETERS_H_
#define CGOGN_MODELING_TYPES_GRAPH_PARAMETERS_H_

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
class GraphParameters
{
	template <typename T>
	using Attribute = typename mesh_traits<GRAPH>::template Attribute<T>;
	using Vertex = typename mesh_traits<GRAPH>::Vertex;
	using Edge = typename mesh_traits<GRAPH>::Edge;

	using Vec3 = geometry::Vec3;

public:

	GraphParameters()
		: vertex_position_(nullptr), vertex_scale_factor_(1.0), 
			sphere_scale_factor_(10.0), selected_vertices_set_(nullptr), 
			object_update_(false), 
			selecting_cell_(SelectingCell::VertexSelect),
		  selection_method_(SelectionMethod::SingleCell), 
		  dragging_handle_(false), dragging_axis_(false)
	{
		param_point_sprite_ = 
			rendering::ShaderPointSprite::generate_param();
		param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
		param_point_sprite_->set_vbos({&selected_vertices_vbo_});
	}

	~GraphParameters()
	{
	}

	
	GraphParameters(GraphParameters&&) = default;
	GraphParameters& operator=(GraphParameters&&) = default;


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

	/// @brief set number of handles in the axis 
	/// update the number of transformations 
	/// number of transformations = (number_of_handles -1) + 2
	/// number_of_handles -1 = number of bones 
	/// +2 for the virtual bones at the extremities 
	/// @param number_of_handles 
	void set_number_of_handles(const size_t& number_of_handles)
	{
		number_of_handles_ = number_of_handles; 

		transformations_.resize(number_of_handles+1);
		for (size_t t = 0; t < transformations_.size(); t++){
			transformations_[t].setIdentity();
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


	/// @brief A key pressed 
	/// initialize displacement of axis/partial skeleton  
	/// @param view current view
	void key_pressed_A_event(ui::View* view)
	{
		if (vertex_position_ && selected_vertices_set_ && 
									selected_vertices_set_->size() > 0)
		{
			if (selected_vertices_set_->size() == 1)
			{
				select_one_vertex(); 
			}
			else
			{
				select_multiple_vertices(); 
			}

			dragging_axis_ = true;
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


	/// @brief a or q key released 
	/// reset dragging axis parameter 
	/// reset the axis transformations 
	/// @param view 
	void key_release_axis_event(ui::View* view)
	{
		if (dragging_axis_){
			dragging_axis_ = false;

			for (size_t t = 0; t < transformations_.size(); t++){
				transformations_[t].setIdentity();
			}

			transformations_valid_indices_.clear(); 
			
		}
			
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

		if (dragging_axis_)
		{
			mouse_axis_displacement(view, x, y); 
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

	size_t number_of_handles_; 

	Vec3 rotation_center_;

	std::vector<rendering::Transfo3d> transformations_;
	std::vector<std::size_t> transformations_valid_indices_; 

	
	float64 drag_z_;
	rendering::GLVec3d previous_drag_pos_;

private: 
	bool dragging_handle_;
	bool dragging_axis_;

	/// @brief select one vertex
	/// Two cases: extremity vertex or not 
	/// TODO: non extremity case
	void select_one_vertex()
	{

		Vertex selected_vertex;  
		selected_vertices_set_->foreach_cell([&](Vertex v) {
			selected_vertex = v; 
		});

		std::vector<Edge>& edges = 
					(*graph_->vertex_incident_edges_)[selected_vertex.index_];

		Edge e0 = edges[0]; 
		std::pair<Vertex, Vertex> first_edge = 
								(*graph_->edge_incident_vertices_)[e0.index_]; 

		if (first_edge.first == selected_vertex){
			rotation_center_ = value<Vec3>(*graph_, vertex_position_, 
															first_edge.second);
		} else {
			rotation_center_ = value<Vec3>(*graph_, vertex_position_, 
															first_edge.first);
		} 

		std::cout << selected_vertex.index_ << std::endl; 
		if (selected_vertex.index_ == 0){
			std::cout << "case 0 index selected" << std::endl; 
			transformations_valid_indices_.push_back(1); 
		} else {
			if (selected_vertex.index_ == number_of_handles_ -1){
				std::cout << "case last index " << std::endl; 
				transformations_valid_indices_.push_back(selected_vertex.index_); 
			} else {
				std::cout << "TODO" << std::endl; 
			}
		}
	}
	

	/// @brief select multiple vertices
	/// Two cases: successive vertices or not 
	/// TODO: not successive vertices 
	void select_multiple_vertices()
	{

		std::vector<Vertex> selected_vertices; 
		selected_vertices_set_->foreach_cell([&](Vertex v) {
			selected_vertices.push_back(v); 
		});

		std::sort(selected_vertices.begin(), selected_vertices.end(), 
				[](Vertex a, Vertex b){return a.index_ < b.index_;});

		if (selected_vertices[0].index_ == 0){
			Vertex last_vertex = 
								selected_vertices[selected_vertices.size()-1]; 

			std::vector<Edge>& edges = 
						(*graph_->vertex_incident_edges_)[last_vertex.index_];
			
			Edge e0 = edges[0];
			Vertex previous_vertex = 
								selected_vertices[selected_vertices.size()-2];

			std::pair<Vertex, Vertex> first_edge = 
								(*graph_->edge_incident_vertices_)[e0.index_]; 
			
			if (first_edge.first == previous_vertex || 
				first_edge.second == previous_vertex)
			{
				Edge e1 = edges[1]; 
				std::pair<Vertex, Vertex> second_edge = 
								(*graph_->edge_incident_vertices_)[e1.index_]; 
				
				if (second_edge.first == last_vertex){
					rotation_center_ = 
						value<Vec3>(*graph_, vertex_position_, second_edge.second);
				} else {
					rotation_center_ = 
						value<Vec3>(*graph_, vertex_position_, second_edge.first);
				}
				
			} else {
				if (first_edge.first == last_vertex){
					rotation_center_ = 
						value<Vec3>(*graph_, vertex_position_, first_edge.second);
				} else {
					rotation_center_ = 
						value<Vec3>(*graph_, vertex_position_, first_edge.first);
				}
			}

			for (std::size_t i = 0; i < selected_vertices.size() -1; i++){
				transformations_valid_indices_.push_back(i+1); 
			}
		} 
		else if (selected_vertices[selected_vertices.size()-1].index_ == 
														number_of_handles_ -1)
		{
			Vertex first_vertex = selected_vertices[0]; 

			std::vector<Edge>& edges = 
						(*graph_->vertex_incident_edges_)[first_vertex.index_];
			
			Edge e0 = edges[0];
			Vertex next_vertex = selected_vertices[1];

			std::pair<Vertex, Vertex> first_edge = 
								(*graph_->edge_incident_vertices_)[e0.index_]; 
			
			if (first_edge.first == next_vertex || 
				first_edge.second == next_vertex)
			{
				Edge e1 = edges[1]; 
				std::pair<Vertex, Vertex> second_edge = 
								(*graph_->edge_incident_vertices_)[e1.index_]; 
				
				if (second_edge.first == first_vertex){
					rotation_center_ = 
						value<Vec3>(*graph_, vertex_position_, second_edge.second);
				} else {
					rotation_center_ = 
						value<Vec3>(*graph_, vertex_position_, second_edge.first);
				}
				
			} else {
				if (first_edge.first == first_vertex){
					rotation_center_ = 
						value<Vec3>(*graph_, vertex_position_, first_edge.second);
				} else {
					rotation_center_ = 
						value<Vec3>(*graph_, vertex_position_, first_edge.first);
				}
			}

			for (std::size_t i = 0; i < selected_vertices.size() -1 ; i++){
				transformations_valid_indices_.push_back(number_of_handles_ - (i+1)); 
			}
		} else {
			std::cout << "TODO" << std::endl; 
		}
		
	}

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

	/// @brief axis displacement caused by mouse displacement
	/// axis rotates around defined center and axis of rotation
	/// @param view current view
	/// @param x retrieved from mouse position
	/// @param y retrieved from mouse position
	void mouse_axis_displacement(ui::View* view, const int32& x, const int32& y)
	{
		float64 dx = float64(x - view->previous_mouse_x());
			float64 dy = float64(y - view->previous_mouse_y());

			if (std::abs(dx) + std::abs(dy) > 0.0)
			{
				rendering::GLVec3d axis(dy, dx, 0.0);
				float64 spinning_speed = axis.norm();

				axis /= spinning_speed;
				spinning_speed *= 0.005;

				float64 sign = 1.0; 
				if (dx < 0.0)
				{
					sign = -1.0;
				}
				
				rendering::Transfo3d inv_camera = 
									view->camera().frame_.inverse();
				rendering::Transfo3d 
				sm(Eigen::AngleAxisd(sign * 1.0 * spinning_speed, normal_)); // 2.0
				rendering::Transfo3d 
					rot((inv_camera * sm * view->camera().frame_).linear());

				rendering::Transfo3d M =
					Eigen::Translation3d(rotation_center_) * rot 
						* Eigen::Translation3d(-rotation_center_);

				selected_vertices_set_->foreach_cell([&](Vertex v) {
					Vec3& axis_vertex_position = value<Vec3>(*graph_, 
															vertex_position_, v);
						axis_vertex_position = M * axis_vertex_position;
					}); 

				for (size_t t = 0; t < transformations_valid_indices_.size(); t++){ 
					transformations_[transformations_valid_indices_[t]] = M;
				}
			}
				
	}



};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_TYPES_GRAPH_PARAMETERS_H_
