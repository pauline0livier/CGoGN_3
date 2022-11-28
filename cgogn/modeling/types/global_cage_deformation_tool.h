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

#ifndef CGOGN_MODELING_GLOBAL_CAGE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_GLOBAL_CAGE_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cmap/cmap2.h>

#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/imgui_helpers.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <boost/synapse/connect.hpp>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class GlobalCageDeformationTool
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

public:
	MESH* global_cage_;
	std::shared_ptr<Attribute<Vec3>> global_cage_vertex_position_;

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_cage_coords_;
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_cage_normal_coords_;

	GlobalCageDeformationTool() : global_cage_vertex_position_(nullptr)
	{
	}

	~GlobalCageDeformationTool()
	{
	}

	void create_global_cage(MESH* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max)
	{
		global_cage_ = m;
		cgogn::modeling::create_bounding_box(*m, vertex_position, bb_min, bb_max);

		global_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		std::shared_ptr<Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, Vertex>(*global_cage_, "vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*global_cage_, vertex_index.get());

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, Vertex>(*global_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*global_cage_, marked_vertices.get());
	}

	void bind_mvc(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position){

		uint32 nbv_object = nb_cells<Vertex>(object);

		uint32 nbv_cage = nb_cells<Vertex>(*global_cage_);

		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, Vertex>(object, "vertex_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*global_cage_, "vertex_index");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked =
			get_attribute<bool, Vertex>(*global_cage_, "marked_vertices");

		global_cage_coords_.resize(nbv_object, nbv_cage); 

		parallel_foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			std::cout << surface_point_idx << std::endl; 
			DartMarker dm(*global_cage_);

			float sumMVC = 0.0;
			for (Dart d = global_cage_->begin(), end = global_cage_->end(); d != end; d = global_cage_->next(d))
			{
				std::cout << "dart session" << std::endl; 
				Vertex cage_vertex = CMap2::Vertex(d);

				bool vc_marked = value<bool>(*global_cage_, cage_vertex_marked, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{

					const Vec3& cage_point = value<Vec3>(*global_cage_, global_cage_vertex_position_, cage_vertex);
					uint32 cage_point_idx = value<uint32>(*global_cage_, cage_vertex_index, cage_vertex);

					float mvc_value = cgogn::modeling::compute_mvc(surface_point, d, *global_cage_, cage_point, global_cage_vertex_position_.get());

					global_cage_coords_(surface_point_idx, cage_point_idx) = mvc_value;

					dm.mark(d);

					value<bool>(*global_cage_, cage_vertex_marked, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			float sum_lambda = 0.0; 

			parallel_foreach_cell(*global_cage_, [&](Vertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(*global_cage_, cage_vertex_index, vc);

				global_cage_coords_(surface_point_idx, cage_point_idx2) = global_cage_coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				sum_lambda += global_cage_coords_(surface_point_idx, cage_point_idx2);

				value<bool>(*global_cage_, cage_vertex_marked, vc) = false;

				return true;
			});

			std::cout << sum_lambda << std::endl; 

			return true;
		});
	}

	void bind_green(MESH& object, CMap2::Attribute<Vec3>* object_vertex_position){

	}


private:
	

};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
