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

#ifndef CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_

#include <cgogn/core/types/cells_set.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>
#include <cgogn/modeling/algos/deformation/deformation_utils.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class SpaceDeformationTool
{
  public:
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

	MESH* influence_cage_;
	cgogn::ui::CellsSet<MESH, Vertex>* influence_area_;

	float smoothing_factor_;

	Eigen::VectorXd attenuation_;

	SpaceDeformationTool()
		: influence_cage_(nullptr), influence_area_(nullptr), influence_cage_vertex_position_(nullptr),
		  smoothing_factor_(0.5f)
	{
	}

	virtual ~SpaceDeformationTool()
	{
	}

	virtual void create_space_tool()
	{
	}

	void set_influence_cage(MESH* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max)
	{
		influence_cage_ = m;
		cgogn::modeling::create_box(*m, vertex_position, bb_min, bb_max);

		influence_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		std::shared_ptr<Attribute<uint32>> position_indices =
			cgogn::add_attribute<uint32, Vertex>(*influence_cage_, "position_indices");
		cgogn::modeling::set_attribute_position_indices(*influence_cage_, position_indices.get());

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, Vertex>(*influence_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*influence_cage_, marked_vertices.get());
	}

	void update_influence_area(const MESH& m, const CMap2::Attribute<Vec3>* vertex_position)
	{
		foreach_cell(m, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(m, vertex_position, v); 

			bool inside_cage = local_mvc_pt_surface(surface_point);

			if (inside_cage)
			{
				influence_area_->select(v);
			}

			return true;
		});
	}

	void bind_mvc(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_position_indices =
			get_attribute<uint32, Vertex>(object, "position_indices");

		std::shared_ptr<Attribute<uint32>> cage_position_indices =
			get_attribute<uint32, Vertex>(*influence_cage_, "position_indices");

		std::shared_ptr<Attribute<bool>> cage_marked_vertices =
			get_attribute<bool, Vertex>(*influence_cage_, "marked_vertices");

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*influence_cage_);

		coords_.resize(nbv_object, nbv_cage);
		coords_.setZero();

		influence_area_->foreach_cell([&](Vertex v) { 
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_position_indices, v);

			DartMarker dm(*influence_cage_);
			float sumMVC = 0.0;

			for (Dart d = influence_cage_->begin(), end = influence_cage_->end(); d != end;
				 d = influence_cage_->next(d))
			{
				Vertex cage_vertex = CMap2::Vertex(d);
				bool vc_marked = value<bool>(*influence_cage_, cage_marked_vertices, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{
					const Vec3& cage_point =
						value<Vec3>(*influence_cage_, influence_cage_vertex_position_, cage_vertex);
					uint32 cage_point_idx = value<uint32>(*influence_cage_, cage_position_indices, cage_vertex);

					float mvc_value = compute_mvc(surface_point, d, *influence_cage_, cage_point,
												  influence_cage_vertex_position_.get());

					coords_(surface_point_idx, cage_point_idx) = mvc_value;

					dm.mark(d);

					value<bool>(*influence_cage_, cage_marked_vertices, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			parallel_foreach_cell(*influence_cage_, [&](Vertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(*influence_cage_, cage_position_indices, vc);

				coords_(surface_point_idx, cage_point_idx2) = coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				value<bool>(*influence_cage_, cage_marked_vertices, vc) = false;

				return true;
			});
		});

		attenuation_.resize(nbv_object);
		attenuation_.setZero();

		compute_attenuation(object);
	}

	void update_deformation_object_mvc(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position, const std::vector<Vec3>& init_position)
	{
		std::cout << "update deformation object" << std::endl; 
		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			cgogn::get_attribute<uint32, Vertex>(object, "position_indices");

		std::shared_ptr<Attribute<uint32>> influence_cage_vertex_index =
			cgogn::get_attribute<uint32, Vertex>(*influence_cage_, "position_indices");

		influence_area_->foreach_cell([&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(object, object_vertex_index, v);

			float current_attenuation = attenuation_(vidx);

			Vec3 new_pos_ = {0.0, 0.0, 0.0};

			foreach_cell(*influence_cage_, [&](Vertex cv) -> bool {
				const Vec3& i_cage_point = value<Vec3>(*influence_cage_, influence_cage_vertex_position_, cv);
				uint32 i_cage_point_idx = value<uint32>(*influence_cage_, influence_cage_vertex_index, cv);

				new_pos_ += coords_(vidx, i_cage_point_idx) * i_cage_point;

				return true;
			});

			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			new_pos_ = new_pos_ * current_attenuation;

			Vec3 current_pos = (1.0f - current_attenuation) * init_position[vidx];

			value<Vec3>(object, object_vertex_position, v) = new_pos_ + current_pos;

			return true;
		});
	}

protected:
	std::shared_ptr<Attribute<Vec3>> influence_cage_vertex_position_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> coords_;

	Eigen::Matrix<Vec2, Eigen::Dynamic, Eigen::Dynamic> n_coords_;

private:
	void compute_attenuation(MESH& object)
	{
		const float h = smoothing_factor_;
		assert((h <= 1.0f && h >= 0.0f) || !"Cage's attenuation factor must be computed!");

		std::shared_ptr<Attribute<uint32>> cage_face_indices =
			add_attribute<uint32, Face>(*influence_cage_, "face_indices");

		cgogn::modeling::set_attribute_face_indices(*influence_cage_, cage_face_indices.get());

		std::shared_ptr<Attribute<uint32>> object_position_indices =
			get_attribute<uint32, Vertex>(object, "position_indices");

		uint32 nbf_cage = 2 * nb_cells<Face>(*influence_cage_);
		uint32 nbv_cage = nb_cells<Vertex>(*influence_cage_);

		influence_area_->foreach_cell([&](Vertex v) {
			uint32 surface_point_idx = value<uint32>(object, object_position_indices, v);

			float i_dist = cage_influence_distance(surface_point_idx, nbf_cage, nbv_cage);

			attenuation_(surface_point_idx) =
				(i_dist < h) ? (0.5f * ((float)sin(M_PI * ((i_dist / h) - 0.5f))) + 0.5f) : 1.0f;
		});
	}

	float cage_influence_distance(const int object_point_idx, const int& nbf_cage, const int& nbv_cage)
	{

		std::shared_ptr<Attribute<uint32>> cage_position_index =
			get_attribute<uint32, Vertex>(*influence_cage_, "position_indices");

		std::shared_ptr<Attribute<uint32>> cage_face_index =
			get_attribute<uint32, Face>(*influence_cage_, "face_indices");

		float r = 1.0;
		foreach_cell(*influence_cage_, [&](Face fc) -> bool {
			uint32 cage_face_idx = value<uint32>(*influence_cage_, cage_face_index, fc);

			std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*influence_cage_, fc);

			// triangle 1
			const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3], face_vertices_[0]};
			const std::vector<uint32> t1_index = {value<uint32>(*influence_cage_, cage_position_index, triangle1[0]),
												  value<uint32>(*influence_cage_, cage_position_index, triangle1[1]),
												  value<uint32>(*influence_cage_, cage_position_index, triangle1[2])};

			const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2], face_vertices_[3]};
			const std::vector<uint32> t2_index = {value<uint32>(*influence_cage_, cage_position_index, triangle2[0]),
												  value<uint32>(*influence_cage_, cage_position_index, triangle2[1]),
												  value<uint32>(*influence_cage_, cage_position_index, triangle2[2])};

			r *= (1.f - (coords_(object_point_idx, t1_index[0]) + coords_(object_point_idx, t1_index[1]) +
						 coords_(object_point_idx, t1_index[2])));

			r *= (1.f - (coords_(object_point_idx, t2_index[0]) + coords_(object_point_idx, t2_index[1]) +
						 coords_(object_point_idx, t2_index[2])));

			return true;
		});

		r /= std::pow((1.0 - (3.0 / nbv_cage)), nbf_cage);

		if (r > 1.0f)
		{
			r = 1.0f;
		}

		return r;
	}

	bool local_mvc_pt_surface(Vec3 pt)
	{
		std::shared_ptr<Attribute<uint32>> i_position_indices =
			get_attribute<uint32, Vertex>(*influence_cage_, "position_indices");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked =
			get_attribute<bool, Vertex>(*influence_cage_, "marked_vertices");

		DartMarker dm(*influence_cage_);

		bool checked = true;
		for (Dart d = influence_cage_->begin(), end = influence_cage_->end(); d != end; d = influence_cage_->next(d))
		{
			Vertex cage_vertex = CMap2::Vertex(d);
			bool vc_marked = value<bool>(*influence_cage_, cage_vertex_marked, cage_vertex);

			if (!dm.is_marked(d) && !vc_marked)
			{

				const Vec3& cage_point =
					value<Vec3>(*influence_cage_,influence_cage_vertex_position_, cage_vertex);

				float mvc_value = cgogn::modeling::compute_mvc(pt, d, *influence_cage_, cage_point, influence_cage_vertex_position_.get());

				dm.mark(d);

				value<bool>(*influence_cage_, cage_vertex_marked, cage_vertex) = true;

				if (mvc_value < 0)
				{
					checked = false;
					break;
				}
			}
		}

		parallel_foreach_cell(*influence_cage_, [&](Vertex vc) -> bool {
			value<bool>(*influence_cage_, cage_vertex_marked, vc) = false;
			return true;
		});

		return checked;
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_SPACE_DEFORMATION_TOOL_H_
