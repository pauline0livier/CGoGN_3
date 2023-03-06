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

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/modeling/algos/deformation/deformation_definitions.h>
#include <cgogn/core/types/cells_set.h>


namespace cgogn
{

namespace modeling
{

template <typename MESH>

class Parameters
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;

	using Vec3 = geometry::Vec3;

public:

	Parameters()
			: vertex_position_(nullptr), selection_method_(SelectionMethod::SingleCell), selecting_cell_(SelectingCell::VertexSelect),
			  selected_vertices_set_(nullptr), vertex_scale_factor_(1.0), sphere_scale_factor_(10.0),
			  object_update_(false), back_selection_(false)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 0.65f);
			param_point_sprite_->set_vbos({&selected_vertices_vbo_});
		}

		~Parameters()
		{
		}


	void update_selected_vertices_vbo()
		{

			if (selected_vertices_set_)
			{
				std::vector<Vec3> selected_vertices_position;
				selected_vertices_position.reserve(selected_vertices_set_->size());
				selected_vertices_set_->foreach_cell([&](Vertex v) {
					selected_vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_, v));
				});
				rendering::update_vbo(selected_vertices_position, &selected_vertices_vbo_);
			}
		}

		//CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		MESH* mesh_;
		std::string name_;

		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		std::vector<Vec3> init_position_;

		std::shared_ptr<boost::synapse::connection> cells_set_connection_;

		bool object_update_;

		bool back_selection_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;

		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
		float32 sphere_scale_factor_;

		rendering::VBO selected_vertices_vbo_;

		std::vector<std::pair<Vertex, Vertex>> selected_depth_vertices_; 

		ui::CellsSet<MESH, Vertex>* selected_vertices_set_;

		SelectingCell selecting_cell_;
		SelectionMethod selection_method_;
	
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_MODELING_TYPES_MESH_PARAMETERS_H_
