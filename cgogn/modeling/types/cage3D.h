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

#ifndef CGOGN_MODELING_CAGE_TOOL_H_
#define CGOGN_MODELING_CAGE_TOOL_H_

#include <cgogn/core/types/cmap/cmap2.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class CageTool : public CMap2
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

public:
	Cage3D()
	{
		cage_ = nullptr;
		vertex_position_ = nullptr;
		normal_position_ = nullptr;
	}

	create_cage(CMap2* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max)
	{
		cage_ = m;
		vertex_position = vertex_position;

		create_box(bb_min, bb_max);

		create_attribute_position_indices();

		create_attribute_marked_vertex();
	}

	create_attribute_face_indices()
	{
		std::shared_ptr<Attribute<uint32>> face_index = add_attribute<uint32, Face>(cage_, "face_indices");
		uint32 nb_faces = 0;
		foreach_cell(cage_, [&](Face f) -> bool {
			value<uint32>(cage_, face_index, f) = nb_faces++;
			return true;
		});
	}

	get_vertex_position()
	{
		return vertex_position_;
	}

	get_normal_position()
	{
		return normal_position_;
	}

private:
	MESH* cage_;
	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	std::shared_ptr<Attribute<Vec3>> normal_position_;

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;

	std::string name;

	void create_box(const Vec3& bb_min, const Vec3& bb_max)
	{
		CMap2::Volume v = add_prism(cage_, 4);
		Dart f1 = v.dart;
		Dart f2 = phi<2, 1, 1, 2>(m, f1);
		std::vector<CMap2::Vertex> vertices = {CMap2::Vertex(f1),
											   CMap2::Vertex(phi1(m, f1)),
											   CMap2::Vertex(phi<1, 1>(m, f1)),
											   CMap2::Vertex(phi_1(m, f1)),
											   CMap2::Vertex(f2),
											   CMap2::Vertex(phi1(m, f2)),
											   CMap2::Vertex(phi<1, 1>(m, f2)),
											   CMap2::Vertex(phi_1(m, f2))};

		Vec3 center = (bb_min + bb_max) / Scalar(2);
		Vec3 bb_min_ = ((bb_min - center) * 1.5) + center;
		Vec3 bb_max_ = ((bb_max - center) * 1.5) + center;

		value<Vec3>(m, vertex_position, vertices[0]) = bb_min_;
		value<Vec3>(m, vertex_position, vertices[1]) = {bb_min_[0], bb_max_[1], bb_min_[2]};
		value<Vec3>(m, vertex_position, vertices[2]) = {bb_max_[0], bb_max_[1], bb_min_[2]};
		value<Vec3>(m, vertex_position, vertices[3]) = {bb_max_[0], bb_min_[1], bb_min_[2]};

		value<Vec3>(m, vertex_position, vertices[4]) = {bb_min_[0], bb_max_[1], bb_max_[2]};
		value<Vec3>(m, vertex_position, vertices[5]) = {bb_min_[0], bb_min_[1], bb_max_[2]};
		value<Vec3>(m, vertex_position, vertices[6]) = {bb_max_[0], bb_min_[1], bb_max_[2]};
		value<Vec3>(m, vertex_position, vertices[7]) = {bb_max_[0], bb_max_[1], bb_max_[2]};
	}

	create_attribute_position_indices()
	{
		std::shared_ptr<Attribute<uint32>> vertex_index = cgogn::add_attribute<uint32, Vertex>(cage_, "pos_indices");

		uint32 nb_vertices_ = 0;
		foreach_cell(*cage_, [&](Vertex v) -> bool {
			value<uint32>(*cage_, vertex_index, v) = nb_vertices_++;
			return true;
		});
	}

	create_attribute_marked_vertex()
	{
		std::shared_ptr<Attribute<bool>> vertex_marked = cgogn::add_attribute<bool, Vertex>(*cage_, "marked_vertex");

		parallel_foreach_cell(*cage_, [&](Vertex v) -> bool {
			value<bool>(*cage_, vertex_marked, v) = false;
			return true;
		});
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_H_
