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

#ifndef CGOGN_MODELING_TYPES_AXIS_PARAMETERS_H_
#define CGOGN_MODELING_TYPES_AXIS_PARAMETERS_H_

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
class AxisParameters
{
	template <typename T>
	using Attribute = typename mesh_traits<GRAPH>::template Attribute<T>;
	using Vertex = typename mesh_traits<GRAPH>::Vertex;
	using Edge = typename mesh_traits<GRAPH>::Edge;

	using Vec3 = geometry::Vec3;

	struct skeletonNode
	{
		bool root; 
		bool selected_node;
		Vertex axis_vertex; 
		rendering::Transfo3d local_transformation; 
	}; 

public:

	AxisParameters()
		: vertex_position_(nullptr), vertex_scale_factor_(1.0), 
			sphere_scale_factor_(10.0), selected_vertices_set_(nullptr), 
			object_update_(false), 
			selecting_cell_(SelectingCell::VertexSelect),
		  selection_method_(SelectionMethod::SingleCell), dragging_axis_(false)
	{
		param_point_sprite_ = 
			rendering::ShaderPointSprite::generate_param();
		param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
		param_point_sprite_->set_vbos({&selected_vertices_vbo_});
	}

	~AxisParameters()
	{
	}

	
	AxisParameters(AxisParameters&&) = default;
	AxisParameters& operator=(AxisParameters&&) = default;


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


	/// @brief number of transformations = (number of handles -1) + 2
	/// +2 for the virtual bones at the extremities 
	/// @param number_of_handles 
	void set_skeleton_data(const size_t& number_of_handles, 
							const std::vector<Vertex> axis_skeleton)
	{
		std::cout << "set skeleton  begin " << std::endl; 
		number_of_handles_ = number_of_handles; 
		size_t number_of_transformations = number_of_handles_ +1; 
		
		skeleton_data_.resize(number_of_transformations); 

		skeletonNode rootNode; 
		rootNode.selected_node = false; 
		rootNode.local_transformation.setIdentity();
		skeleton_data_[0] = rootNode; 

		for (std::size_t t = 1; t < number_of_transformations -1; t++)
		{
			skeletonNode currentNode; 
			currentNode.axis_vertex = axis_skeleton[t-1];   
			currentNode.selected_node = false; 
			currentNode.local_transformation.setIdentity(); 
			skeleton_data_[t] = currentNode; 
		}

		skeletonNode lastNode; 
		lastNode.selected_node = true; 
		lastNode.local_transformation.setIdentity(); 
		skeleton_data_[number_of_transformations-1] = lastNode; 

	}

	/// @brief reset transformations of axis
	void reset_transformations()
	{
		for (size_t t = 0; t < skeleton_data_.size(); t++)
		{
			skeleton_data_[t].local_transformation.setIdentity();
			skeleton_data_[t].selected_node = false;
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


	/// @brief a or q key released 
	/// reset dragging axis parameter 
	/// reset the axis transformations 
	/// @param view 
	void key_release_axis_event(ui::View* view)
	{
		if (dragging_axis_){
			dragging_axis_ = false;

			reset_transformations(); 
			
		}
			
	}

	/// @brief mouse displacement outcomes 
	/// call for the corresponding displacement functions 
	/// @param view current view
	/// @param x retrieved from mouse position
	/// @param y retrieved from mouse position
	void mouse_displacement(ui::View* view, const int32& x, const int32& y)
	{
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

	std::vector<skeletonNode> skeleton_data_;  

private: 
	bool dragging_axis_;

	/// @brief select one vertex
	/// Two cases: extremity vertex or not 
	/// TODO: non extremity case
	void select_one_vertex()
	{
		uint32 vertex_index; 
		selected_vertices_set_->foreach_cell([&](Vertex v) {
			vertex_index = v.index_; 
			skeleton_data_[vertex_index].selected_node = true;
		});

		std::cout << "single begin " << std::endl; 
		if (vertex_index == 0)
		{
			Vertex next_vertex = skeleton_data_[1].axis_vertex; 

			rotation_center_ =  value<Vec3>(*graph_, vertex_position_, next_vertex);
		} else {
			Vertex previous_vertex = skeleton_data_[vertex_index -1].axis_vertex; 

			rotation_center_ =  value<Vec3>(*graph_, vertex_position_, previous_vertex);
		}
		std::cout << "single end " << std::endl; 

		/*Vertex selected_vertex;  
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
 
		skeleton_data_[selected_vertex.index].selected_node = true;*/  
	}
	

	/// @brief select multiple vertices
	/// Two cases: successive vertices or not 
	/// TODO: not successive vertices 
	void select_multiple_vertices()
	{

		/*std::vector<Vertex> selected_vertices; 
		selected_vertices_set_->foreach_cell([&](Vertex v) {
			selected_vertices.push_back(v); 
		});*/ 
		selected_vertices_set_->foreach_cell([&](Vertex v) {
			skeleton_data_[v.index_].selected_node = true; 
		});

		std::cout << "valid here" << std::endl; 
		for (std::size_t n = 0; n < skeleton_data_.size() -1; n++)
		{
			std::cout << skeleton_data_[n].selected_node << std::endl; 
			std::cout << skeleton_data_[n+1].selected_node << std::endl; 
			if ((!skeleton_data_[n].selected_node) && skeleton_data_[n+1].selected_node)
			{
				rotation_center_ = value<Vec3>(*graph_, vertex_position_, skeleton_data_[n].axis_vertex);
				break; 
			}
		}

		std::cout << "multiple end " << std::endl; 

		/*std::sort(selected_vertices.begin(), selected_vertices.end(), 
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

			for (std::size_t i = 0; i < selected_vertices.size() +1; i++){  
				transformations_valid_indices_.push_back(number_of_handles_ -i); 
			}

		} else {
			std::cout << "TODO" << std::endl; 
		}*/
		
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
				sm(Eigen::AngleAxisd(sign * 1.0 * spinning_speed, normal_)); 

				rendering::Transfo3d 
					rot((inv_camera * sm * view->camera().frame_).linear());

				rendering::Transfo3d M =
					Eigen::Translation3d(rotation_center_) * rot 
						* Eigen::Translation3d(-rotation_center_);

				std::cout << "check seg fault " << std::endl; 
				std::size_t n = 1; 
				while ( n < skeleton_data_.size() && !skeleton_data_[n].selected_node){
					n += 1; 
				}
 
				if (n == skeleton_data_.size()){
					std::cout << "no selection " << std::endl; 
				} else {
					skeletonNode first_selected_node = skeleton_data_[n]; 
					skeleton_data_[n].local_transformation = M; 
			
					Vertex first_vertex = skeleton_data_[n].axis_vertex; 
					Vec3& first_vertex_position = value<Vec3>(*graph_, 
															vertex_position_, first_vertex);

					value<Vec3>(*graph_, vertex_position_, first_vertex) = M*first_vertex_position; 
 
					if (n < skeleton_data_.size() -2)
					{


						for (std::size_t c = n+1; c < skeleton_data_.size() -1; c++){
							if (!skeleton_data_[c].selected_node)
							{
								//rendering::Transfo3d inverse_M = M.inverse();  

								/*Vertex local_vertex = skeleton_data_[c].axis_vertex; 
								Vec3& local_vertex_position = value<Vec3>(*graph_, 
															vertex_position_, local_vertex);

								value<Vec3>(*graph_, vertex_position_, local_vertex) = skeleton_data_[c].inverse_M*local_vertex_position;*/
								break; 
							} else {

								Vertex local_vertex = skeleton_data_[c].axis_vertex; 
								Vec3& local_vertex_position = value<Vec3>(*graph_, 
															vertex_position_, local_vertex);

								skeleton_data_[c].local_transformation = M; 

								value<Vec3>(*graph_, vertex_position_, local_vertex) = M*local_vertex_position;

							}
						}
					}

				}

				
				/*selected_vertices_set_->foreach_cell([&](Vertex v) {
					Vec3& axis_vertex_position = value<Vec3>(*graph_, 
															vertex_position_, v);
						axis_vertex_position = M * axis_vertex_position;


					}); */

				/*for (size_t t = 0; t < transformations_valid_indices_.size(); t++){ 
					//transformations_[transformations_valid_indices_[t]] = M;
					transformations_[transformations_valid_indices_[t]] = rot;

				}*/

				/*transformations_[1] = rot; 
				transformations_[2] = rot; 
				transformations_[3] = rot;*/

			}
				
	}



};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_TYPES_AXIS_PARAMETERS_H_
