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

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/modeling/algos/deformation/creation_space_tool.h>

namespace cgogn
{

namespace modeling
{

std::vector<CMap2::Vertex> create_bounding_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min,
											   const Vec3& bb_max)
{
	CMap2::Volume v = add_prism(m, 4);
	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi_1(m, f2))};

	value<Vec3>(m, vertex_position, vertices[0]) = bb_min;
	value<Vec3>(m, vertex_position, vertices[1]) = {bb_min[0], bb_max[1], bb_min[2]};
	value<Vec3>(m, vertex_position, vertices[2]) = {bb_max[0], bb_max[1], bb_min[2]};
	value<Vec3>(m, vertex_position, vertices[3]) = {bb_max[0], bb_min[1], bb_min[2]};

	value<Vec3>(m, vertex_position, vertices[4]) = {bb_min[0], bb_max[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[5]) = {bb_min[0], bb_min[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[6]) = {bb_max[0], bb_min[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[7]) = {bb_max[0], bb_max[1], bb_max[2]};

	return vertices;
}

void update_bounding_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, const std::vector<CMap2::Vertex>& vertices,
						 const Vec3& bb_min, const Vec3& bb_max)
{

	value<Vec3>(m, vertex_position, vertices[0]) = bb_min;
	value<Vec3>(m, vertex_position, vertices[1]) = {bb_min[0], bb_max[1], bb_min[2]};
	value<Vec3>(m, vertex_position, vertices[2]) = {bb_max[0], bb_max[1], bb_min[2]};
	value<Vec3>(m, vertex_position, vertices[3]) = {bb_max[0], bb_min[1], bb_min[2]};

	value<Vec3>(m, vertex_position, vertices[4]) = {bb_min[0], bb_max[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[5]) = {bb_min[0], bb_min[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[6]) = {bb_max[0], bb_min[1], bb_max[2]};
	value<Vec3>(m, vertex_position, vertices[7]) = {bb_max[0], bb_max[1], bb_max[2]};
}

void create_cage_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max,
					 const Vec3& center, const Vec3& normal)
{
	CMap2::Volume v = add_prism(m, 4);
	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {
		CMap2::Vertex(f1), CMap2::Vertex(phi1(m, f1)), CMap2::Vertex(phi<1, 1>(m, f1)), CMap2::Vertex(phi_1(m, f1)),
		CMap2::Vertex(f2), CMap2::Vertex(phi1(m, f2)), CMap2::Vertex(phi<1, 1>(m, f2)), CMap2::Vertex(phi_1(m, f2))};

	// bbmin belongs to the plane
	float d_min = -(normal.dot(bb_min));

	// find alpha such that center + alpha*handle_normal belongs to plane
	float alpha_min = -d_min - (normal.dot(center));

	Eigen::Vector3d center_min_plane = center + alpha_min * normal;

	Vec3 local_frame_x = bb_min - center_min_plane;
	local_frame_x.normalize();

	Vec3 local_frame_y = normal.cross(local_frame_x);
	local_frame_y.normalize();

	Eigen::Matrix3d local_frame, frame_inverse;

	local_frame.row(0) = local_frame_x;
	local_frame.row(1) = local_frame_y;
	local_frame.row(2) = normal;
	frame_inverse = local_frame.inverse();

	const Vec3 local_bb_min = local_frame * (bb_min - center);
	const Vec3 local_bb_max = local_frame * (bb_max - center);

	const double radius = local_bb_min.norm();

	double min_n, max_n;
	if (local_bb_min[2] > local_bb_max[2])
	{
		min_n = local_bb_max[2];
		max_n = local_bb_min[2];
	}
	else
	{
		min_n = local_bb_min[2];
		max_n = local_bb_max[2];
	}

	Eigen::Vector3d local_vertex0 = {radius * std::cos(0), radius * std::sin(0), max_n};
	Eigen::Vector3d local_vertex1 = {radius * std::cos(M_PI / 2), radius * std::sin(M_PI / 2), max_n};
	Eigen::Vector3d local_vertex2 = {radius * std::cos(M_PI), radius * std::sin(M_PI), max_n};
	Eigen::Vector3d local_vertex3 = {radius * std::cos(3 * M_PI / 2), radius * std::sin(3 * M_PI / 2), max_n};

	Eigen::Vector3d local_vertex4 = {radius * std::cos(-3 * M_PI / 2), radius * std::sin(-3 * M_PI / 2), min_n};
	Eigen::Vector3d local_vertex5 = {radius * std::cos(0), radius * std::sin(0), min_n};
	Eigen::Vector3d local_vertex6 = {radius * std::cos(-M_PI / 2), radius * std::sin(-M_PI / 2), min_n};
	Eigen::Vector3d local_vertex7 = {radius * std::cos(-M_PI), radius * std::sin(-M_PI), min_n};

	value<Vec3>(m, vertex_position, vertices[0]) = (frame_inverse * local_vertex0) + center;
	value<Vec3>(m, vertex_position, vertices[1]) = (frame_inverse * local_vertex1) + center;
	value<Vec3>(m, vertex_position, vertices[2]) = frame_inverse * local_vertex2 + center;
	value<Vec3>(m, vertex_position, vertices[3]) = frame_inverse * local_vertex3 + center;
	value<Vec3>(m, vertex_position, vertices[4]) = frame_inverse * local_vertex4 + center;
	value<Vec3>(m, vertex_position, vertices[5]) = frame_inverse * local_vertex5 + center;
	value<Vec3>(m, vertex_position, vertices[6]) = frame_inverse * local_vertex6 + center;
	value<Vec3>(m, vertex_position, vertices[7]) = frame_inverse * local_vertex7 + center;
}

void create_handle_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, CMap2::Attribute<Vec3>* local_vertex_position,
					   const Vec3& handle_position, const double& radius, const Eigen::Matrix3d& frame_inverse,
					   const float& local_max_depth, const float& local_min_depth)
{

	CMap2::Volume v = add_prism(m, 6);

	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {CMap2::Vertex(f1),
										   CMap2::Vertex(phi1(m, f1)),
										   CMap2::Vertex(phi<1, 1>(m, f1)),
										   CMap2::Vertex(phi<1, 1, 1>(m, f1)),
										   CMap2::Vertex(phi<1, 1, 1, 1>(m, f1)),
										   CMap2::Vertex(phi_1(m, f1)),
										   CMap2::Vertex(f2),
										   CMap2::Vertex(phi1(m, f2)),
										   CMap2::Vertex(phi<1, 1>(m, f2)),
										   CMap2::Vertex(phi<1, 1, 1>(m, f2)),
										   CMap2::Vertex(phi<1, 1, 1, 1>(m, f2)),
										   CMap2::Vertex(phi_1(m, f2))};

	Eigen::Vector3d local_vertex0 = {radius * std::cos(0), radius * std::sin(0), local_max_depth};

	Eigen::Vector3d local_vertex1 = {radius * std::cos(M_PI / 3), radius * std::sin(M_PI / 3), local_max_depth};

	Eigen::Vector3d local_vertex2 = {radius * std::cos(2 * M_PI / 3), radius * std::sin(2 * M_PI / 3), local_max_depth};

	Eigen::Vector3d local_vertex3 = {radius * std::cos(M_PI), radius * std::sin(M_PI), local_max_depth};

	Eigen::Vector3d local_vertex4 = {radius * std::cos(4 * M_PI / 3), radius * std::sin(4 * M_PI / 3), local_max_depth};

	Eigen::Vector3d local_vertex5 = {radius * std::cos(5 * M_PI / 3), radius * std::sin(5 * M_PI / 3), local_max_depth};

	value<Vec3>(m, local_vertex_position, vertices[0]) = local_vertex0;
	value<Vec3>(m, local_vertex_position, vertices[1]) = local_vertex1;
	value<Vec3>(m, local_vertex_position, vertices[2]) = local_vertex2;
	value<Vec3>(m, local_vertex_position, vertices[3]) = local_vertex3;
	value<Vec3>(m, local_vertex_position, vertices[4]) = local_vertex4;
	value<Vec3>(m, local_vertex_position, vertices[5]) = local_vertex5;

	value<Vec3>(m, vertex_position, vertices[0]) = (frame_inverse * local_vertex0) + handle_position;
	value<Vec3>(m, vertex_position, vertices[1]) = (frame_inverse * local_vertex1) + handle_position;
	value<Vec3>(m, vertex_position, vertices[2]) = frame_inverse * local_vertex2 + handle_position;
	value<Vec3>(m, vertex_position, vertices[3]) = frame_inverse * local_vertex3 + handle_position;
	value<Vec3>(m, vertex_position, vertices[4]) = frame_inverse * local_vertex4 + handle_position;
	value<Vec3>(m, vertex_position, vertices[5]) = frame_inverse * local_vertex5 + handle_position;

	Eigen::Vector3d local_vertex6 = {radius * std::cos(-5 * M_PI / 3), radius * std::sin(-5 * M_PI / 3),
									 local_min_depth};
	Eigen::Vector3d local_vertex7 = {radius * std::cos(0), radius * std::sin(0), local_min_depth};
	Eigen::Vector3d local_vertex8 = {radius * std::cos(-M_PI / 3), radius * std::sin(-M_PI / 3), local_min_depth};
	Eigen::Vector3d local_vertex9 = {radius * std::cos(-2 * M_PI / 3), radius * std::sin(-2 * M_PI / 3),
									 local_min_depth};
	Eigen::Vector3d local_vertex10 = {radius * std::cos(-M_PI), radius * std::sin(-M_PI), local_min_depth};
	Eigen::Vector3d local_vertex11 = {radius * std::cos(-4 * M_PI / 3), radius * std::sin(-4 * M_PI / 3),
									  local_min_depth};

	value<Vec3>(m, local_vertex_position, vertices[6]) = local_vertex6;
	value<Vec3>(m, local_vertex_position, vertices[7]) = local_vertex7;
	value<Vec3>(m, local_vertex_position, vertices[8]) = local_vertex8;
	value<Vec3>(m, local_vertex_position, vertices[9]) = local_vertex9;
	value<Vec3>(m, local_vertex_position, vertices[10]) = local_vertex10;
	value<Vec3>(m, local_vertex_position, vertices[11]) = local_vertex11;

	value<Vec3>(m, vertex_position, vertices[6]) = frame_inverse * local_vertex6 + handle_position;
	value<Vec3>(m, vertex_position, vertices[7]) = frame_inverse * local_vertex7 + handle_position;
	value<Vec3>(m, vertex_position, vertices[8]) = frame_inverse * local_vertex8 + handle_position;
	value<Vec3>(m, vertex_position, vertices[9]) = frame_inverse * local_vertex9 + handle_position;
	value<Vec3>(m, vertex_position, vertices[10]) = frame_inverse * local_vertex10 + handle_position;
	value<Vec3>(m, vertex_position, vertices[11]) = frame_inverse * local_vertex11 + handle_position;
}

void create_axis_box(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, CMap2::Attribute<Vec3>* local_vertex_position,
					 Graph::Attribute<uint32>* skeleton_vertex, const std::vector<Vec3>& vertex_coords,
					 const std::vector<Vec3>& local_vertex_coords)
{

	int volume_size = vertex_coords.size() / 2;
	CMap2::Volume v = add_prism(m, volume_size);

	Dart f1 = v.dart;
	Dart f2 = phi<2, 1, 1, 2>(m, f1);
	std::vector<CMap2::Vertex> vertices = {CMap2::Vertex(f1),
										   CMap2::Vertex(phi1(m, f1)),
										   CMap2::Vertex(phi<1, 1>(m, f1)),
										   CMap2::Vertex(phi<1, 1, 1>(m, f1)),
										   CMap2::Vertex(phi<1, 1, 1, 1>(m, f1)),
										   CMap2::Vertex(phi_1(m, f1)),
										   CMap2::Vertex(f2),
										   CMap2::Vertex(phi1(m, f2)),
										   CMap2::Vertex(phi<1, 1>(m, f2)),
										   CMap2::Vertex(phi<1, 1, 1>(m, f2)),
										   CMap2::Vertex(phi<1, 1, 1, 1>(m, f2)),
										   CMap2::Vertex(phi_1(m, f2))};

	// top
	value<uint32>(m, skeleton_vertex, vertices[0]) = 0;
	value<uint32>(m, skeleton_vertex, vertices[1]) = 0;
	value<uint32>(m, skeleton_vertex, vertices[2]) = 1;
	value<uint32>(m, skeleton_vertex, vertices[3]) = 2;
	value<uint32>(m, skeleton_vertex, vertices[4]) = 2;
	value<uint32>(m, skeleton_vertex, vertices[5]) = 1;

	value<Vec3>(m, local_vertex_position, vertices[0]) = local_vertex_coords[4 * 0 + 0];
	value<Vec3>(m, local_vertex_position, vertices[1]) = local_vertex_coords[4 * 0 + 1];
	value<Vec3>(m, local_vertex_position, vertices[2]) = local_vertex_coords[4 * 1 + 1];
	value<Vec3>(m, local_vertex_position, vertices[3]) = local_vertex_coords[4 * 2 + 1];
	value<Vec3>(m, local_vertex_position, vertices[4]) = local_vertex_coords[4 * 2 + 0];
	value<Vec3>(m, local_vertex_position, vertices[5]) = local_vertex_coords[4 * 1 + 0];

	value<Vec3>(m, vertex_position, vertices[0]) = vertex_coords[0];
	value<Vec3>(m, vertex_position, vertices[1]) = vertex_coords[1];
	value<Vec3>(m, vertex_position, vertices[2]) = vertex_coords[5];
	value<Vec3>(m, vertex_position, vertices[3]) = vertex_coords[9];
	value<Vec3>(m, vertex_position, vertices[4]) = vertex_coords[8];
	value<Vec3>(m, vertex_position, vertices[5]) = vertex_coords[4];

	// bottom
	value<uint32>(m, skeleton_vertex, vertices[6]) = 0;
	value<uint32>(m, skeleton_vertex, vertices[7]) = 0;
	value<uint32>(m, skeleton_vertex, vertices[8]) = 1;
	value<uint32>(m, skeleton_vertex, vertices[9]) = 2;
	value<uint32>(m, skeleton_vertex, vertices[10]) = 2;
	value<uint32>(m, skeleton_vertex, vertices[11]) = 1;

	value<Vec3>(m, local_vertex_position, vertices[6]) = local_vertex_coords[4 * 0 + 3];
	value<Vec3>(m, local_vertex_position, vertices[7]) = local_vertex_coords[4 * 0 + 2];
	value<Vec3>(m, local_vertex_position, vertices[8]) = local_vertex_coords[4 * 1 + 2];
	value<Vec3>(m, local_vertex_position, vertices[9]) = local_vertex_coords[4 * 2 + 2];
	value<Vec3>(m, local_vertex_position, vertices[10]) = local_vertex_coords[4 * 2 + 3];
	value<Vec3>(m, local_vertex_position, vertices[11]) = local_vertex_coords[4 * 1 + 3];

	value<Vec3>(m, vertex_position, vertices[6]) = vertex_coords[3];
	value<Vec3>(m, vertex_position, vertices[7]) = vertex_coords[2];
	value<Vec3>(m, vertex_position, vertices[8]) = vertex_coords[6];
	value<Vec3>(m, vertex_position, vertices[9]) = vertex_coords[10];
	value<Vec3>(m, vertex_position, vertices[10]) = vertex_coords[11];
	value<Vec3>(m, vertex_position, vertices[11]) = vertex_coords[7];
}

void set_attribute_vertex_index(CMap2& cage, CMap2::Attribute<uint32>* position_indices)
{
	uint32 nb_vertices = 0;

	foreach_cell(cage, [&](CMap2::Vertex v) -> bool {
		value<uint32>(cage, position_indices, v) = nb_vertices++;
		return true;
	});
}

void set_attribute_vertex_index_graph(IncidenceGraph& graph, IncidenceGraph::Attribute<uint32>* position_indices)
{
	uint32 nb_vertices = 0;

	foreach_cell(graph, [&](IncidenceGraph::Vertex v) -> bool {
		value<uint32>(graph, position_indices, v) = nb_vertices++;
		return true;
	});
}

void set_attribute_marked_vertices(CMap2& cage, CMap2::Attribute<bool>* marked_vertices)
{

	parallel_foreach_cell(cage, [&](CMap2::Vertex v) -> bool {
		value<bool>(cage, marked_vertices, v) = false;
		return true;
	});
}

void set_attribute_face_indices(CMap2& cage, CMap2::Attribute<uint32>* face_indices)
{
	uint32 nb_faces = 0;
	foreach_cell(cage, [&](CMap2::Face f) -> bool {
		value<uint32>(cage, face_indices, f) = nb_faces++;
		return true;
	});
}

void set_attribute_face_normal(CMap2& cage, CMap2::Attribute<Vec3>* vertex_position,
							   CMap2::Attribute<Vec3>* face_normal)
{
	foreach_cell(cage, [&](CMap2::Face fc) -> bool {
		std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(cage, fc);

		// triangle 1

		const std::vector<CMap2::Vertex> triangle1 = {face_vertices_[1], face_vertices_[3], face_vertices_[0]};
		const std::vector<Vec3> t1_values = {value<Vec3>(cage, vertex_position, triangle1[0]),
											 value<Vec3>(cage, vertex_position, triangle1[1]),
											 value<Vec3>(cage, vertex_position, triangle1[2])};

		//Vec3 t1_normal = (cgogn::geometry::normal(t1_values[0], t1_values[1], t1_values[2])).normalized();

		// triangle 2
		/*const std::vector<CMap2::Vertex> triangle2 = {face_vertices_[1], face_vertices_[2], face_vertices_[3]};
		const std::vector<Vec3> t2_values = {value<Vec3>(cage, vertex_position, triangle2[0]),
											 value<Vec3>(cage, vertex_position, triangle2[1]),
											 value<Vec3>(cage, vertex_position, triangle2[2])};*/

		//Vec3 t2_normal = (cgogn::geometry::normal(t2_values[0], t2_values[1], t2_values[2])).normalized();

		value<Vec3>(cage, face_normal, fc) = (cgogn::geometry::normal(t1_values[0], t1_values[1], t1_values[2])).normalized(); 

		return true;
	} );
}

Graph::Vertex create_handle(Graph& g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius,
							const Vec3& center1)
{
	Graph::Vertex nv = add_vertex(g);

	value<Vec3>(g, vertex_position, nv) = center1;
	value<Scalar>(g, vertex_radius, nv) = Scalar(5);

	return nv;
}

std::vector<Graph::Vertex> create_axis(Graph& g, Graph::Attribute<Vec3>* vertex_position,
									   Graph::Attribute<Scalar>* vertex_radius,
									   const std::vector<Vec3>& vertices_positions)
{
	std::vector<Graph::Vertex> list_vertex;
	Graph::Vertex nv = add_vertex(g);

	list_vertex.push_back(nv);

	value<Vec3>(g, vertex_position, nv) = vertices_positions[0];
	value<Scalar>(g, vertex_radius, nv) = Scalar(5);

	Graph::Vertex lastVertex = nv;

	for (unsigned i = 1; i < vertices_positions.size(); i++)
	{
		Graph::Vertex nv1 = add_vertex(g);
		list_vertex.push_back(nv1);

		value<Vec3>(g, vertex_position, nv1) = vertices_positions[i];
		value<Scalar>(g, vertex_radius, nv1) = Scalar(5);

		connect_vertices(g, lastVertex, nv1);

		lastVertex = nv1;
	}

	return list_vertex;
}

void set_graph_attribute_position_indices(Graph& g, Graph::Attribute<uint32>* position_indices)
{
	uint32 nb_vertices = 0;

	foreach_cell(g, [&](Graph::Vertex v) -> bool {
		value<uint32>(g, position_indices, v) = nb_vertices++;
		return true;
	});
}

} // namespace modeling

} // namespace cgogn